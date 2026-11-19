"""Request body size limit middleware (pure ASGI)."""

from __future__ import annotations

from collections.abc import Awaitable, Callable, MutableMapping
from typing import Any

Scope = MutableMapping[str, Any]
Message = MutableMapping[str, Any]
Receive = Callable[[], Awaitable[Message]]
Send = Callable[[Message], Awaitable[None]]
ASGIApp = Callable[[Scope, Receive, Send], Awaitable[None]]

DEFAULT_MAX_REQUEST_BYTES = 1024 * 1024  # 1 MiB
# The OAuth endpoints take a handful of short form/JSON fields; nothing legitimate is large.
AUTH_MAX_REQUEST_BYTES = 64 * 1024
AUTH_PATHS = frozenset({"/register", "/token", "/revoke", "/authorize"})


class BodySizeLimitMiddleware:
    """Reject HTTP requests whose body exceeds a byte limit with ``413``.

    Two enforcement points, so neither a lying nor a missing ``Content-Length`` gets past:

    * **Declared length.** A ``Content-Length`` over the limit is rejected before the app
      runs (no body is read). A malformed or conflicting ``Content-Length`` is a ``400``.
    * **Streamed length.** ``receive`` is wrapped and counts body bytes as they arrive
      (chunked uploads have no ``Content-Length``). On overflow the middleware answers ``413``
      itself if the app has not started a response, then hands the app an
      ``http.disconnect`` so it stops reading; anything the app sends afterwards is dropped.

    The OAuth endpoints (``/register``, ``/token``, ``/revoke``, ``/authorize``) get a fixed
    64 KiB cap; every other path uses ``max_bytes``.
    """

    def __init__(
        self,
        app: ASGIApp,
        max_bytes: int = DEFAULT_MAX_REQUEST_BYTES,
        auth_max_bytes: int = AUTH_MAX_REQUEST_BYTES,
        auth_paths: frozenset[str] = AUTH_PATHS,
    ) -> None:
        self.app = app
        self.max_bytes = max_bytes
        self.auth_max_bytes = auth_max_bytes
        self.auth_paths = auth_paths

    def limit_for(self, path: str) -> int:
        normalized = path.rstrip("/") or "/"
        if normalized in self.auth_paths:
            return min(self.auth_max_bytes, self.max_bytes)
        return self.max_bytes

    @staticmethod
    async def _respond(send: Send, status: int, text: str) -> None:
        body = text.encode()
        await send(
            {
                "type": "http.response.start",
                "status": status,
                "headers": [
                    (b"content-type", b"text/plain; charset=utf-8"),
                    (b"content-length", str(len(body)).encode()),
                    (b"connection", b"close"),
                ],
            }
        )
        await send({"type": "http.response.body", "body": body})

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return

        limit = self.limit_for(scope.get("path", ""))

        declared: set[bytes] = {
            value.strip()
            for name, value in scope.get("headers", [])
            if name.lower() == b"content-length"
        }
        if declared:
            if len(declared) != 1 or not next(iter(declared)).isdigit():
                await self._respond(send, 400, "Invalid Content-Length")
                return
            if int(next(iter(declared))) > limit:
                await self._respond(send, 413, "Request body too large")
                return

        received = 0
        response_started = False
        rejected = False
        answered = False

        async def limited_receive() -> Message:
            nonlocal received, rejected, answered
            if rejected:
                return {"type": "http.disconnect"}
            message = await receive()
            if message["type"] == "http.request":
                received += len(message.get("body", b""))
                if received > limit:
                    rejected = True
                    if not response_started:
                        answered = True
                        await self._respond(send, 413, "Request body too large")
                    return {"type": "http.disconnect"}
            return message

        async def guarded_send(message: Message) -> None:
            nonlocal response_started
            if answered:
                return  # we already sent the 413; the app's response is discarded
            if message["type"] == "http.response.start":
                response_started = True
            await send(message)

        try:
            await self.app(scope, limited_receive, guarded_send)
        except Exception:
            # Once we have answered 413, the app's reaction to the synthetic disconnect
            # (e.g. Starlette's ClientDisconnect from request.body()) is expected noise.
            if not answered:
                raise

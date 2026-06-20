---
name: provider-oci
description: How Oracle Cloud (OCI) works in a Greenlight setup — the `target: oci` runtime for stateful MCP servers (BAMCP), the Always-Free idle-reclaim trap solved by converting to PAYG (manual, not code), health-check keepalive, and API-key signing. Use when wiring/debugging an oci-target tool or the idle-reclaim problem.
---

# provider-oci

OCI is the `target: oci` runtime for **stateful** services that don't fit serverless —
the canonical case is **BAMCP** (a stateful MCP server). The deploy adapter for oci is a
follow-up (not built); today oci tools deploy via their own repo's CI, and the wrapper
manages DNS (a Cloudflare CNAME / Tunnel) + keepalive.

## The idle-reclaim trap — fixed manually, NOT by code

OCI **Always Free** reclaims idle compute. The harness does **not** prevent this. The fix is
**manual**: convert the tenancy to **Pay-As-You-Go** (+ set a low billing alarm so it stays
effectively free). The keepalive worker only **health-checks** the service and **nags** via
`doctor` / the alert sink — it cannot stop reclaim. Never imply otherwise.

See `docs/oci-payg-runbook.md` for the PAYG conversion + billing-alarm steps.

## Auth — API-key signing (not a bearer token)

OCI uses **request signing** with an API key (`~/.oci/config` + a PEM key), not a bearer
token — so there's no cheap `verify()` curl. Setup is the runbook; Greenlight treats the
credential as presence-only.

## Terraform

The `tool` module manages the subdomain DNS (CNAME / Tunnel hostname) for the oci service.
The compute itself is out of the V1 declarative scope (the oci deploy adapter is deferred).

## MCP / verify
An oci-hosted tool is typically an **MCP server** — verify it with `mode: mcp` (initialize →
`tools/list` → call a tool & assert shape → assert auth rejection), and connect at
`<name>.<domain>/mcp` (streamable HTTP).

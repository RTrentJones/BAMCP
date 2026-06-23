---
name: provider-oci
description: How Oracle Cloud (OCI) works in a Greenlight setup — the `target: oci` runtime for stateful MCP servers (BAMCP) on a free-tier Ampere A1 Container Instance, the provider-agnostic build-via-GitHub→GHCR model, Greenlight-owned compute + tunnel Terraform, the OCI token CLI, deploy = restart, and staying on the free tier (no PAYG; recover-on-alert). Use when wiring/debugging an oci-target tool.
---

# provider-oci

OCI is the `target: oci` runtime for **stateful** services that don't fit serverless — the
canonical case is **BAMCP** (a stateful MCP server). The split: **the tool is provider-agnostic
and just builds a container via GitHub; Greenlight owns the OCI infra** (compute + tunnel + DNS),
configured in the wrapper.

## Free tier — A1 Container Instance + GHCR

The Always-Free path is an **OCI Container Instance** on **Ampere A1** (the A1 allotment — 2 OCPU
/ 12 GB as of 2026-06-15 — is shared across VM / Bare-Metal / Container-Instances). No VM to
provision, no cloud-init, no SSH. The image comes from **GHCR** (free); **OCI's own registry
(OCIR) is paid** — that was the trap in the old BAMCP pipeline. The container instance runs the
tool container + a **cloudflared sidecar** (shared netns → localhost), exposed at `<name>.<domain>`.

## Provider-agnostic tool → GHCR

The tool repo (BAMCP) has ONE job: a GitHub Actions workflow that **builds + pushes the container
to GHCR**. No OCI, no SSH, no deploy logic — portable to any provider's infra.

## Greenlight OCI infra (Terraform, in the wrapper)

`greenlight add/adopt` emits, per oci tool:
- **`oci-container-instance`** module — the container instance (tool image from GHCR + cloudflared
  sidecar with the tunnel token), `CI.Standard.A1.Flex` within the free allotment, restart ALWAYS.
- **`tunnel`** module — cloudflared tunnel + ingress `<name>.<domain> → http://localhost:8000` + token.
- **`tool`** DNS module — CNAME → the tunnel.
The `oci` provider (auth below) is added to `infra/main.tf`.

## OCI token CLI

`greenlight secrets gather <tool> --repo <o/r>` pushes the OCI creds straight to GitHub secrets
(hidden prompts, no disk/logs). **The only manual OCI inputs are the API-key auth values** —
`TF_VAR_OCI_TENANCY_OCID`, `TF_VAR_OCI_USER_OCID`, `TF_VAR_OCI_FINGERPRINT`, `TF_VAR_OCI_PRIVATE_KEY`
(PEM), `TF_VAR_OCI_REGION`. `TF_VAR_OCI_COMPARTMENT_ID` is **optional** (blank → the tenancy/root
compartment). Auth is API-key request signing — no bearer, so no fetch-verify. The container
instance OCID is **not** a manual input — the deploy workflow resolves it at deploy time from OCI
by the instance's display name (= the tool name), so it's abstracted from the developer.

**The VCN, subnet, and availability domain are NOT manual** — they're Terraform: the `oci-network`
module creates the VCN + a public (egress-only) subnet, and the container-instance module looks the
AD up via an `oci_identity_availability_domains` data source. So the bootstrap is just "create one
API key" — Terraform can't create the credential it needs to authenticate, but it owns everything
after that. (Out-of-A1-capacity in one AD? set `availability_domain` on the instance module to pin
another — the only time you touch it.)

**Shortcut — feed the API-key config preview directly.** After *Add API key*, OCI shows a
"Configuration file preview" (the `[DEFAULT]` block) and you download the `.pem`. Pass both:
`greenlight secrets gather <tool> --repo <o/r> --oci-config ~/path/config [--oci-key ~/path/key.pem]`
— it auto-fills the 5 auth secrets (incl. the multi-line PEM, read from the file so it's never
pasted) and only prompts for the 3 placement values + the option-B deploy PATs.

## Deploy = restart (re-pull)

`greenlight deploy <tool>` (oci) just runs `oci container-instances container-instance restart
--container-instance-id <OCID>` — the instance re-pulls the latest GHCR image. The tool's CI
builds; an event trigger (the chosen deploy option) fires the restart. The adapter does NOT build.

## Idle-reclaim — stay free, recover on alert

OCI Always-Free can reclaim idle compute. We **stay on the free tier** and accept that: the
instance runs restart-policy ALWAYS, keepalive health-checks it, and if it's ever stopped/reclaimed
the alert fires and a **re-apply / redeploy** restores it. **PAYG is NOT used** — it's an optional
last resort (see `docs/oci-payg-runbook.md`) only if reclaim ever becomes a recurring problem.

## Verify
The tool is typically an **MCP server**: verify with `mode: mcp`, connect at `<name>.<domain>/mcp`
(FastMCP's `streamable_http_app()` serves `/mcp` by default — the convention). If auth gates
`initialize`, supply a token or use an `api`-mode 401 check. Keepalive health-checks `target: oci`.

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

`greenlight add`/`init` gather the OCI creds into `.greenlight/secrets.env` (+ GH secrets):
**provider auth** `TF_VAR_oci_tenancy_ocid`, `TF_VAR_oci_user_ocid`, `TF_VAR_oci_fingerprint`,
`TF_VAR_oci_private_key` (PEM), `TF_VAR_oci_region`; **placement** `TF_VAR_oci_compartment_id`,
`TF_VAR_oci_availability_domain`, `TF_VAR_oci_subnet_id`; and `OCI_CONTAINER_INSTANCE_OCID`
(the Terraform output, for deploy). Auth is API-key request signing — no bearer, so no fetch-verify.

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

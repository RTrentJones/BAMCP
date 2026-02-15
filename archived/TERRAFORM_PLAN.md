# Plan: Migrate BAMCP Infrastructure to Terraform + GitHub Actions

## Context

BAMCP's OCI infrastructure (container instance with BAMCP + cloudflared sidecar, OCIR repo) is currently managed by interactive shell scripts (`setup-cloudflared.sh`, `setup-oci-secrets.sh`). This makes infrastructure hard to reproduce, audit, or change safely. Moving to Terraform + GitHub Actions CI gives us declarative config, plan-before-apply reviews, and automated deployments.

## Architecture

**Two separate workflows, two separate concerns:**

```
Code change (v* tag / manual)          Infra change (push to infra/)
         │                                        │
    deploy.yml                               infra.yml
         │                                        │
  Build ARM64 image                     terraform plan (on PR)
  Push to OCIR                          terraform apply (on merge)
  Restart container instance            Recreates instance if needed
```

- **deploy.yml** (high frequency): Build image, push `:latest`, restart instance. No Terraform. Fast, no downtime.
- **infra.yml** (low frequency): `terraform plan` on PRs touching `infra/`, `terraform apply` on merge to main. Only runs when infrastructure config changes (env vars, shape, new containers).

This separation matters because OCI Container Instances are **immutable** — Terraform changes force destroy/recreate, but `restart` just pulls the latest image with no downtime.

**State backend**: OCI Object Storage with S3-compatible API. State stored remotely so GitHub Actions can access it.

## File Structure

```
infra/
  main.tf              — Provider, backend, data sources, all resources
  variables.tf         — Input variable declarations
  outputs.tf           — Instance OCID, OCIR repo path
  versions.tf          — Required providers + Terraform version
  terraform.tfvars.example — Documented example (no secrets, committed)

.github/workflows/
  infra.yml            — NEW: Terraform plan on PR, apply on merge
  deploy.yml           — MODIFIED: Add concurrency group
  ci.yml               — Unchanged
  release.yml          — Unchanged

scripts/
  bootstrap-state.sh   — NEW: One-time OCI Object Storage bucket + S3 creds setup
  debug-container.sh   — KEPT: Operational debugging (unchanged)
  setup-cloudflared.sh — ARCHIVED: Move to archived/scripts/
  setup-oci-secrets.sh — ARCHIVED: Move to archived/scripts/

.gitignore             — MODIFIED: Add Terraform patterns
CLAUDE.md              — MODIFIED: Document infra/ and Terraform workflow
```

Flat structure (no modules) — only 2 resources, modules would be overengineering.

## Terraform Resources

### `infra/versions.tf`

```hcl
terraform {
  required_version = ">= 1.5.0"
  required_providers {
    oci = {
      source  = "oracle/oci"
      version = "~> 6.0"
    }
  }
  backend "s3" {
    # Configured via -backend-config at init time
    # See scripts/bootstrap-state.sh
  }
}
```

### `infra/main.tf`

```hcl
provider "oci" {
  tenancy_ocid = var.tenancy_ocid
  user_ocid    = var.user_ocid
  fingerprint  = var.fingerprint
  private_key  = var.private_key
  region       = var.region
}

data "oci_identity_availability_domains" "ads" {
  compartment_id = var.tenancy_ocid
}

data "oci_objectstorage_namespace" "ns" {
  compartment_id = var.tenancy_ocid
}

# --- OCIR Repository ---
resource "oci_artifacts_container_repository" "bamcp" {
  compartment_id = var.compartment_ocid
  display_name   = "bamcp"
  is_public      = false
  is_immutable   = false
}

# --- Container Instance (BAMCP + cloudflared sidecar) ---
resource "oci_container_instances_container_instance" "bamcp" {
  compartment_id           = var.compartment_ocid
  availability_domain      = data.oci_identity_availability_domains.ads.availability_domains[var.availability_domain_index].name
  display_name             = "bamcp"
  shape                    = "CI.Standard.A1.Flex"
  container_restart_policy = "ALWAYS"

  shape_config {
    ocpus         = 1
    memory_in_gbs = 2
  }

  containers {
    display_name          = "bamcp"
    image_url             = "${var.region}.ocir.io/${data.oci_objectstorage_namespace.ns.namespace}/bamcp:latest"
    environment_variables = var.bamcp_env_vars
  }

  containers {
    display_name = "cloudflared"
    image_url    = "docker.io/cloudflare/cloudflared:latest"
    arguments    = ["tunnel", "run", "--token", var.cloudflare_tunnel_token]
  }

  vnics {
    subnet_id             = var.subnet_ocid
    is_public_ip_assigned = false
  }

  image_pull_secrets {
    registry_endpoint = "${var.region}.ocir.io"
    secret_type       = "BASIC"
    username          = base64encode(var.ocir_username)
    password          = base64encode(var.ocir_auth_token)
  }
}
```

### `infra/variables.tf`

All secrets passed via `TF_VAR_*` environment variables in GitHub Actions (mapped from existing GitHub Secrets):

| Variable | Source Secret | Sensitive |
|----------|-------------|-----------|
| `tenancy_ocid` | `OCI_CLI_TENANCY` | no |
| `user_ocid` | `OCI_CLI_USER` | no |
| `fingerprint` | `OCI_CLI_FINGERPRINT` | no |
| `private_key` | `OCI_CLI_KEY_CONTENT` | yes |
| `region` | `OCI_CLI_REGION` | no |
| `compartment_ocid` | `OCI_COMPARTMENT_OCID` | no |
| `subnet_ocid` | `OCI_SUBNET_OCID` (new) | no |
| `ocir_username` | `OCIR_USERNAME` (new) | yes |
| `ocir_auth_token` | `OCI_AUTH_TOKEN` | yes |
| `cloudflare_tunnel_token` | `CLOUDFLARE_TUNNEL_TOKEN` | yes |
| `availability_domain_index` | default `0` | no |
| `bamcp_env_vars` | defaults in variables.tf | no |

The `bamcp_env_vars` map has defaults matching current config (BAMCP_TRANSPORT=streamable-http, etc.), so it doesn't need a secret.

### `infra/outputs.tf`

```hcl
output "container_instance_id" {
  value = oci_container_instances_container_instance.bamcp.id
}
output "ocir_repo_path" {
  value = oci_artifacts_container_repository.bamcp.display_name
}
```

## GitHub Actions Workflows

### `infra.yml` (new)

```
Triggers: push to main (paths: infra/**), PR to main (paths: infra/**), manual dispatch
```

- **Plan job** (runs on all triggers):
  1. Checkout, setup Terraform 1.9
  2. `terraform init` with `-backend-config` flags (S3 endpoint, bucket, credentials from secrets)
  3. `terraform validate`
  4. `terraform plan -no-color -out=tfplan`
  5. On PRs: post plan output as PR comment

- **Apply job** (runs only on push to main, requires plan success):
  1. Same init
  2. `terraform apply -auto-approve`
  3. Print new container instance OCID in job summary (for manual secret update if changed)

- **Environment gate**: Apply job uses `environment: production` requiring manual approval

### `deploy.yml` (modified)

Only change: add concurrency group to prevent race conditions with infra workflow:

```yaml
concurrency:
  group: bamcp-oci
  cancel-in-progress: false
```

### New GitHub Secrets Needed

| Secret | Purpose | How to Get |
|--------|---------|-----------|
| `OCI_SUBNET_OCID` | Private subnet for container instance | OCI Console or `oci network subnet list` |
| `OCIR_USERNAME` | OCIR pull credentials (`namespace/email`) | Already known from setup-cloudflared.sh runs |
| `OCI_OS_NAMESPACE` | Object Storage namespace for backend URL | `oci os ns get` |
| `OCI_S3_ACCESS_KEY` | S3-compatible access key for TF state | OCI Console > Identity > Customer Secret Keys |
| `OCI_S3_SECRET_KEY` | S3-compatible secret key for TF state | Same (shown once at creation) |

## Bootstrap (One-Time, Manual)

### `scripts/bootstrap-state.sh`

1. Create OCI Object Storage bucket `bamcp-tfstate` with versioning enabled
2. Print instructions to create Customer Secret Keys (S3-compatible credentials)
3. Print `gh secret set` commands for the new secrets

### Import Existing Resources

After bootstrap + `terraform init`:

```bash
cd infra
terraform import oci_artifacts_container_repository.bamcp <OCIR_REPO_OCID>
terraform import oci_container_instances_container_instance.bamcp <INSTANCE_OCID>
terraform plan  # Verify minimal drift
```

First `terraform apply` after import may force a recreate (image_pull_secrets encoding diff). Schedule during a maintenance window — takes ~2 minutes, cloudflared reconnects automatically.

## Implementation Steps

1. Create `infra/versions.tf`, `main.tf`, `variables.tf`, `outputs.tf`, `terraform.tfvars.example`
2. Create `scripts/bootstrap-state.sh`
3. Create `.github/workflows/infra.yml`
4. Modify `.github/workflows/deploy.yml` (add concurrency group)
5. Add Terraform patterns to `.gitignore`
6. Move `scripts/setup-cloudflared.sh` and `scripts/setup-oci-secrets.sh` to `archived/scripts/`
7. Update `CLAUDE.md` deployment section

## Verification

1. Run `bootstrap-state.sh` to create state bucket and S3 credentials
2. Push new secrets to GitHub (`OCI_SUBNET_OCID`, `OCIR_USERNAME`, `OCI_OS_NAMESPACE`, `OCI_S3_ACCESS_KEY`, `OCI_S3_SECRET_KEY`)
3. `cd infra && terraform init -backend-config=backend.hcl` (local test)
4. `terraform import` existing resources
5. `terraform plan` — should show minimal or no changes
6. Push to a branch, open PR — verify `infra.yml` posts plan as PR comment
7. Merge — verify apply succeeds, container instance is still running
8. Tag a new version — verify `deploy.yml` still builds/pushes/restarts correctly
9. `curl -I https://bamcp.rtrentjones.dev/mcp` — still works

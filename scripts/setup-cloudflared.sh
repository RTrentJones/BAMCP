#!/usr/bin/env bash
# ============================================================================
# Setup Cloudflare Tunnel sidecar for BAMCP OCI Container Instance.
#
# Recreates the container instance with two containers:
#   1. BAMCP (from OCIR)
#   2. cloudflared (Cloudflare Tunnel connector)
#
# Prerequisites:
#   - OCI CLI configured (`oci setup config`)
#   - GitHub CLI authenticated (`gh auth login`)
#   - Cloudflare Tunnel created (Zero Trust -> Networks -> Tunnels)
#   - BAMCP image pushed to OCIR
#
# Usage:
#   ./scripts/setup-cloudflared.sh [--profile PROFILE] [--dry-run]
# ============================================================================

set -euo pipefail

PROFILE="DEFAULT"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --profile) PROFILE="$2"; shift 2 ;;
        --dry-run)  DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: $0 [--profile PROFILE] [--dry-run]"
            echo "  --profile   OCI config profile (default: DEFAULT)"
            echo "  --dry-run   Print config without creating resources"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------- dependency checks ----------

for cmd in oci gh; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: '$cmd' is required but not found in PATH." >&2
        exit 1
    fi
done

# ---------- read OCI config ----------

OCI_CONFIG="${OCI_CONFIG_FILE:-$HOME/.oci/config}"
if [[ ! -f "$OCI_CONFIG" ]]; then
    echo "Error: OCI config not found at $OCI_CONFIG" >&2
    echo "Run 'oci setup config' first." >&2
    exit 1
fi

read_config() {
    local key="$1"
    awk -v profile="[$PROFILE]" -v key="$key" '
        $0 == profile { found=1; next }
        /^\[/ { found=0 }
        found && $0 ~ "^"key"[[:space:]]*=" {
            sub(/^[^=]*=[[:space:]]*/, "")
            print
            exit
        }
    ' "$OCI_CONFIG"
}

TENANCY_OCID=$(read_config "tenancy")
REGION=$(read_config "region")
KEY_FILE=$(read_config "key_file")
KEY_FILE="${KEY_FILE/#\~/$HOME}"

for var in TENANCY_OCID REGION; do
    if [[ -z "${!var}" ]]; then
        echo "Error: Could not read '$var' from [$PROFILE] in $OCI_CONFIG" >&2
        exit 1
    fi
done

# ---------- handle encrypted key passphrase ----------

PASSPHRASE_FILE=""
cleanup() {
    if [[ -n "$PASSPHRASE_FILE" && -f "$PASSPHRASE_FILE" ]]; then
        rm -f "$PASSPHRASE_FILE"
    fi
    if [[ -n "${CONTAINER_JSON_FILE:-}" && -f "${CONTAINER_JSON_FILE:-}" ]]; then
        rm -f "$CONTAINER_JSON_FILE"
    fi
    if [[ -n "${REGISTRY_CREDS_FILE:-}" && -f "${REGISTRY_CREDS_FILE:-}" ]]; then
        rm -f "$REGISTRY_CREDS_FILE"
    fi
}
trap cleanup EXIT

OCI_PASS_ARGS=()
if [[ -f "$KEY_FILE" ]] && head -1 "$KEY_FILE" | grep -q "ENCRYPTED"; then
    echo "Private key is encrypted."
    read -rsp "Enter passphrase: " PASSPHRASE
    echo

    PASSPHRASE_FILE=$(mktemp)
    chmod 600 "$PASSPHRASE_FILE"
    cat > "$PASSPHRASE_FILE" <<OCIEOF
[$PROFILE]
user=$(read_config "user")
tenancy=$TENANCY_OCID
fingerprint=$(read_config "fingerprint")
key_file=$KEY_FILE
region=$REGION
pass_phrase=$PASSPHRASE
OCIEOF
    OCI_PASS_ARGS=(--config-file "$PASSPHRASE_FILE")
    unset PASSPHRASE
fi

run_oci() {
    oci "${OCI_PASS_ARGS[@]}" "$@"
}

# ---------- verify OCI CLI auth ----------

echo "Verifying OCI CLI authentication..."
if ! run_oci iam region list --query "data[0].name" --raw-output &>/dev/null; then
    echo "Error: OCI CLI authentication failed." >&2
    exit 1
fi
echo "  OK"
echo

# ---------- prompt for tunnel token and public URL ----------

read -rsp "Paste Cloudflare Tunnel token (input hidden): " TUNNEL_TOKEN
echo
if [[ -z "$TUNNEL_TOKEN" ]]; then
    echo "Error: Tunnel token cannot be empty." >&2
    exit 1
fi

DEFAULT_PUBLIC_URL="https://bamcp.rtrentjones.dev"
read -rp "Public URL [$DEFAULT_PUBLIC_URL]: " PUBLIC_URL_INPUT
BAMCP_PUBLIC_URL="${PUBLIC_URL_INPUT:-$DEFAULT_PUBLIC_URL}"
# Strip trailing slash
BAMCP_PUBLIC_URL="${BAMCP_PUBLIC_URL%/}"
echo

# ---------- select compartment ----------

echo "Fetching compartments..."
COMPARTMENT_NAMES=("(root/tenancy)")
COMPARTMENT_IDS=("$TENANCY_OCID")

while IFS=$'\t' read -r name id; do
    [[ -z "$name" ]] && continue
    COMPARTMENT_NAMES+=("$name")
    COMPARTMENT_IDS+=("$id")
done < <(run_oci iam compartment list \
    --compartment-id "$TENANCY_OCID" \
    --compartment-id-in-subtree true \
    --lifecycle-state ACTIVE \
    --query 'data[].{a:name, b:id}' \
    --output table 2>/dev/null | tail -n +3 | grep -v '^+' | sed 's/  *| */\t/g; s/^|//; s/|$//')

echo "Available compartments:"
for i in "${!COMPARTMENT_NAMES[@]}"; do
    echo "  [$i] ${COMPARTMENT_NAMES[$i]}"
done
echo
read -rp "Select compartment [0-$((${#COMPARTMENT_NAMES[@]} - 1))]: " COMP_IDX
COMPARTMENT_OCID="${COMPARTMENT_IDS[$COMP_IDX]}"
echo "  Selected: ${COMPARTMENT_NAMES[$COMP_IDX]}"
echo

# ---------- select availability domain ----------

echo "Fetching availability domains..."
AD_NAMES=()
while IFS= read -r name; do
    [[ -z "$name" ]] && continue
    AD_NAMES+=("$name")
done < <(run_oci iam availability-domain list \
    --compartment-id "$TENANCY_OCID" \
    --query 'data[].name' \
    --raw-output 2>/dev/null | tr -d '[]"' | tr ',' '\n' | sed 's/^ *//')

if [[ ${#AD_NAMES[@]} -eq 1 ]]; then
    AD_NAME="${AD_NAMES[0]}"
    echo "  Using: $AD_NAME"
else
    echo "Available ADs:"
    for i in "${!AD_NAMES[@]}"; do
        echo "  [$i] ${AD_NAMES[$i]}"
    done
    echo
    read -rp "Select AD [0-$((${#AD_NAMES[@]} - 1))]: " AD_IDX
    AD_NAME="${AD_NAMES[$AD_IDX]}"
    echo "  Selected: $AD_NAME"
fi
echo

# ---------- select subnet ----------

echo "Fetching subnets..."
SUBNET_NAMES=()
SUBNET_IDS=()

while IFS=$'\t' read -r name id; do
    [[ -z "$name" ]] && continue
    SUBNET_NAMES+=("$name")
    SUBNET_IDS+=("$id")
done < <(run_oci network subnet list \
    --compartment-id "$COMPARTMENT_OCID" \
    --lifecycle-state AVAILABLE \
    --query 'data[].{a:"display-name", b:id}' \
    --output table 2>/dev/null | tail -n +3 | grep -v '^+' | sed 's/  *| */\t/g; s/^|//; s/|$//')

if [[ ${#SUBNET_NAMES[@]} -eq 0 ]]; then
    echo "Error: No subnets found in compartment." >&2
    exit 1
fi

echo "Available subnets:"
for i in "${!SUBNET_NAMES[@]}"; do
    echo "  [$i] ${SUBNET_NAMES[$i]}"
done
echo
read -rp "Select subnet [0-$((${#SUBNET_NAMES[@]} - 1))]: " SUB_IDX
SUBNET_OCID="${SUBNET_IDS[$SUB_IDX]}"
echo "  Selected: ${SUBNET_NAMES[$SUB_IDX]}"
echo

# ---------- detect OCIR image URL ----------

NAMESPACE=$(run_oci os ns get --query 'data' --raw-output 2>/dev/null)
DEFAULT_IMAGE="${REGION}.ocir.io/${NAMESPACE}/bamcp:latest"
echo "OCIR image URL:"
echo "  Default: $DEFAULT_IMAGE"
read -rp "  Use this? [Y/n] or paste custom URL: " IMAGE_INPUT

if [[ -z "$IMAGE_INPUT" || "$IMAGE_INPUT" == [yY] ]]; then
    OCIR_IMAGE="$DEFAULT_IMAGE"
else
    OCIR_IMAGE="$IMAGE_INPUT"
fi
echo "  Using: $OCIR_IMAGE"
echo

# ---------- OCIR registry credentials ----------

OCIR_USER=$(read_config "user")
OCIR_SERVER="${REGION}.ocir.io"
OCIR_USERNAME="${NAMESPACE}/$(run_oci iam user get \
    --user-id "$OCIR_USER" \
    --query 'data.name' \
    --raw-output 2>/dev/null)"

echo "OCIR pull credentials:"
echo "  Registry: $OCIR_SERVER"
echo "  Username: $OCIR_USERNAME"
read -rsp "  Auth token (same one used for Docker login): " OCIR_AUTH_TOKEN
echo
echo

# ---------- check for existing instance ----------

echo "Checking for existing 'bamcp' container instance..."
EXISTING_ID=$(run_oci container-instances container-instance list \
    --compartment-id "$COMPARTMENT_OCID" \
    --display-name "bamcp" \
    --lifecycle-state ACTIVE \
    --query 'data[0].id' \
    --raw-output 2>/dev/null || true)

if [[ -n "$EXISTING_ID" && "$EXISTING_ID" != "None" && "$EXISTING_ID" != "null" ]]; then
    echo "  Found existing instance: ${EXISTING_ID:0:40}..."
    read -rp "  Delete and recreate? [y/N] " CONFIRM_DELETE
    if [[ "$CONFIRM_DELETE" != [yY] ]]; then
        echo "Aborted."
        exit 0
    fi

    if ! $DRY_RUN; then
        echo "  Deleting..."
        run_oci container-instances container-instance delete \
            --container-instance-id "$EXISTING_ID" \
            --force 2>/dev/null

        echo "  Waiting for deletion..."
        for i in $(seq 1 30); do
            STATE=$(run_oci container-instances container-instance get \
                --container-instance-id "$EXISTING_ID" \
                --query 'data."lifecycle-state"' \
                --raw-output 2>/dev/null || echo "DELETED")
            if [[ "$STATE" == "DELETED" || "$STATE" == "null" ]]; then
                echo "  Deleted."
                break
            fi
            sleep 5
        done
    else
        echo "  [dry-run] Would delete instance."
    fi
else
    echo "  No existing instance found."
fi
echo

# ---------- build container instance JSON ----------

CONTAINER_JSON_FILE=$(mktemp)
chmod 600 "$CONTAINER_JSON_FILE"
cat > "$CONTAINER_JSON_FILE" <<JSONEOF
[
  {
    "displayName": "bamcp",
    "imageUrl": "$OCIR_IMAGE",
    "environmentVariables": {
      "BAMCP_TRANSPORT": "streamable-http",
      "BAMCP_HOST": "0.0.0.0",
      "BAMCP_PORT": "8000",
      "BAMCP_AUTH_ENABLED": "true",
      "BAMCP_ALLOW_REMOTE_FILES": "true",
      "BAMCP_RATE_LIMIT": "60",
      "BAMCP_MAX_READS": "10000",
      "BAMCP_DEFAULT_WINDOW": "500",
      "BAMCP_MIN_VAF": "0.1",
      "BAMCP_MIN_DEPTH": "10",
      "BAMCP_MIN_MAPQ": "0"
    }
  },
  {
    "displayName": "cloudflared",
    "imageUrl": "docker.io/cloudflare/cloudflared:latest",
    "arguments": ["tunnel", "run", "--token", "$TUNNEL_TOKEN"]
  }
]
JSONEOF

VNICS_JSON="[{\"subnetId\": \"$SUBNET_OCID\", \"isPublicIpAssigned\": false}]"

OCIR_USERNAME_B64=$(printf '%s' "$OCIR_USERNAME" | base64 -w0)
OCIR_AUTH_TOKEN_B64=$(printf '%s' "$OCIR_AUTH_TOKEN" | base64 -w0)

REGISTRY_CREDS_FILE=$(mktemp)
chmod 600 "$REGISTRY_CREDS_FILE"
cat > "$REGISTRY_CREDS_FILE" <<REGEOF
[
  {
    "secretType": "BASIC",
    "registryEndpoint": "$OCIR_SERVER",
    "username": "$OCIR_USERNAME_B64",
    "password": "$OCIR_AUTH_TOKEN_B64"
  }
]
REGEOF

# ---------- summary ----------

echo "============================================"
echo "Container instance to create:"
echo "============================================"
echo "  Display name:    bamcp"
echo "  Shape:           CI.Standard.A1.Flex (1 OCPU, 2 GB)"
echo "  AD:              $AD_NAME"
echo "  Subnet:          ${SUBNET_NAMES[$SUB_IDX]}"
echo "  Public IP:       no"
echo "  Restart policy:  ALWAYS"
echo "  Containers:"
echo "    [1] bamcp         — $OCIR_IMAGE"
echo "    [2] cloudflared   — docker.io/cloudflare/cloudflared:latest"
echo "  Tunnel token:    ${TUNNEL_TOKEN:0:20}..."
echo "============================================"
echo

if $DRY_RUN; then
    echo "[dry-run] Would create the above. Exiting."
    exit 0
fi

read -rp "Create this container instance? [y/N] " CONFIRM
if [[ "$CONFIRM" != [yY] ]]; then
    echo "Aborted."
    exit 0
fi

# ---------- create container instance ----------

echo
echo "Creating container instance..."
CREATE_OUTPUT=$(run_oci container-instances container-instance create \
    --compartment-id "$COMPARTMENT_OCID" \
    --availability-domain "$AD_NAME" \
    --shape "CI.Standard.A1.Flex" \
    --shape-config '{"ocpus": 1, "memoryInGBs": 2}' \
    --display-name "bamcp" \
    --container-restart-policy "ALWAYS" \
    --containers "file://$CONTAINER_JSON_FILE" \
    --vnics "$VNICS_JSON" \
    --image-pull-secrets "file://$REGISTRY_CREDS_FILE" \
    2>&1) || true

echo "$CREATE_OUTPUT"

NEW_INSTANCE_ID=$(echo "$CREATE_OUTPUT" | grep -o 'ocid1\.computecontainerinstance\.[^ "]*' | head -1)

if [[ -z "$NEW_INSTANCE_ID" ]]; then
    echo "Error: Failed to create container instance." >&2
    exit 1
fi

echo "  Instance: ${NEW_INSTANCE_ID:0:50}..."

# ---------- wait for ACTIVE state ----------

echo "  Waiting for ACTIVE state..."
for i in $(seq 1 60); do
    STATE=$(run_oci container-instances container-instance get \
        --container-instance-id "$NEW_INSTANCE_ID" \
        --query 'data."lifecycle-state"' \
        --raw-output 2>/dev/null || echo "UNKNOWN")
    if [[ "$STATE" == "ACTIVE" ]]; then
        echo "  Instance is ACTIVE."
        break
    elif [[ "$STATE" == "FAILED" ]]; then
        echo "Error: Instance entered FAILED state." >&2
        exit 1
    fi
    sleep 5
done
echo

# ---------- update GitHub secrets ----------

echo "Updating GitHub secrets..."
echo "$TUNNEL_TOKEN" | gh secret set CLOUDFLARE_TUNNEL_TOKEN
echo "  CLOUDFLARE_TUNNEL_TOKEN"

echo "$NEW_INSTANCE_ID" | gh secret set OCI_CONTAINER_INSTANCE_OCID
echo "  OCI_CONTAINER_INSTANCE_OCID"

echo "$BAMCP_PUBLIC_URL" | gh secret set BAMCP_PUBLIC_URL
echo "  BAMCP_PUBLIC_URL"
echo

# ---------- done ----------

echo "============================================"
echo "Setup complete!"
echo "============================================"
echo
echo "Next steps:"
echo "  1. Check Cloudflare dashboard — tunnel should show 'Healthy'"
echo "  2. Test:  curl -I ${BAMCP_PUBLIC_URL}/mcp"
echo "  3. Configure Claude Desktop:"
echo "     {\"mcpServers\": {\"bamcp\": {\"url\": \"${BAMCP_PUBLIC_URL}/mcp\"}}}"
echo
echo "Instance OCID: $NEW_INSTANCE_ID"

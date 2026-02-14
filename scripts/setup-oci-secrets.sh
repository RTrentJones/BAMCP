#!/usr/bin/env bash
# ============================================================================
# Setup OCI secrets for GitHub Actions deploy workflow.
#
# Reads from ~/.oci/config, interactively selects compartment and container
# instance, then pushes all required secrets via `gh secret set`.
#
# Prerequisites:
#   - OCI CLI configured (`oci setup config`)
#   - GitHub CLI authenticated (`gh auth login`)
#
# Usage:
#   ./scripts/setup-oci-secrets.sh [--profile PROFILE] [--dry-run]
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
            echo "  --dry-run   Print secrets without pushing to GitHub"
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
    # Parse INI-style config: find [PROFILE] section, then key=value
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

USER_OCID=$(read_config "user")
TENANCY_OCID=$(read_config "tenancy")
FINGERPRINT=$(read_config "fingerprint")
KEY_FILE=$(read_config "key_file")
REGION=$(read_config "region")

# Expand ~ in key_file path
KEY_FILE="${KEY_FILE/#\~/$HOME}"

for var in USER_OCID TENANCY_OCID FINGERPRINT KEY_FILE REGION; do
    if [[ -z "${!var}" ]]; then
        echo "Error: Could not read '$var' from [$PROFILE] in $OCI_CONFIG" >&2
        exit 1
    fi
done

if [[ ! -f "$KEY_FILE" ]]; then
    echo "Error: Key file not found: $KEY_FILE" >&2
    exit 1
fi

KEY_BYTES=$(wc -c < "$KEY_FILE")

echo "OCI config [$PROFILE]:"
echo "  User:        ${USER_OCID:0:40}..."
echo "  Tenancy:     ${TENANCY_OCID:0:40}..."
echo "  Fingerprint: $FINGERPRINT"
echo "  Region:      $REGION"
echo "  Key file:    $KEY_FILE ($KEY_BYTES bytes)"
echo

# ---------- handle encrypted key passphrase ----------

PASSPHRASE_FILE=""
cleanup() {
    if [[ -n "$PASSPHRASE_FILE" && -f "$PASSPHRASE_FILE" ]]; then
        rm -f "$PASSPHRASE_FILE"
    fi
}
trap cleanup EXIT

OCI_PASS_ARGS=()
if head -1 "$KEY_FILE" | grep -q "ENCRYPTED"; then
    echo "Private key is encrypted."
    read -rsp "Enter passphrase (asked once, used for all OCI calls): " PASSPHRASE
    echo
    echo

    # Write passphrase to temp file readable only by us, pass via --config-file
    # with a temporary profile that includes pass_phrase
    PASSPHRASE_FILE=$(mktemp)
    chmod 600 "$PASSPHRASE_FILE"
    cat > "$PASSPHRASE_FILE" <<OCIEOF
[$PROFILE]
user=$USER_OCID
tenancy=$TENANCY_OCID
fingerprint=$FINGERPRINT
key_file=$KEY_FILE
region=$REGION
pass_phrase=$PASSPHRASE
OCIEOF
    OCI_PASS_ARGS=(--config-file "$PASSPHRASE_FILE")
    # Clear passphrase from memory
    unset PASSPHRASE
fi

# Wrapper: runs oci cli with passphrase config if needed
run_oci() {
    oci "${OCI_PASS_ARGS[@]}" "$@"
}

# ---------- verify OCI CLI auth works ----------

echo "Verifying OCI CLI authentication..."
if ! run_oci iam region list --query "data[0].name" --raw-output &>/dev/null; then
    echo "Error: OCI CLI authentication failed." >&2
    echo "Check your config, key file, and passphrase." >&2
    exit 1
fi
echo "  Authentication OK."
echo

# ---------- select compartment ----------

echo "Fetching compartments..."

# Use --output tsv to avoid jq dependency. Query returns "name\tid" per line.
COMPARTMENT_LINES=()
COMPARTMENT_IDS=()
COMPARTMENT_NAMES=()

# Always include the root (tenancy) compartment as option 0
COMPARTMENT_NAMES+=("(root/tenancy)")
COMPARTMENT_IDS+=("$TENANCY_OCID")

# Fetch sub-compartments via OCI CLI with TSV output
while IFS=$'\t' read -r name id; do
    [[ -z "$name" ]] && continue
    COMPARTMENT_NAMES+=("$name")
    COMPARTMENT_IDS+=("$id")
done < <(run_oci iam compartment list \
    --compartment-id "$TENANCY_OCID" \
    --compartment-id-in-subtree true \
    --lifecycle-state ACTIVE \
    --query 'data[].{a:name, b:id}' \
    --output table 2>/dev/null | tail -n +3 | sed 's/  *| */\t/g; s/^|//; s/|$//')
# Note: --output table gives header+separator+rows. tail -n +3 skips header.
# sed converts table pipes to tabs for clean parsing.

COMPARTMENT_COUNT=${#COMPARTMENT_NAMES[@]}

if [[ "$COMPARTMENT_COUNT" -eq 0 ]]; then
    echo "Error: No compartments found." >&2
    exit 1
fi

echo "Available compartments:"
for i in "${!COMPARTMENT_NAMES[@]}"; do
    echo "  [$i] ${COMPARTMENT_NAMES[$i]}"
done
echo

read -rp "Select compartment [0-$((COMPARTMENT_COUNT - 1))]: " COMP_IDX

if [[ "$COMP_IDX" -lt 0 || "$COMP_IDX" -ge "$COMPARTMENT_COUNT" ]] 2>/dev/null; then
    echo "Error: Invalid selection." >&2
    exit 1
fi

COMPARTMENT_OCID="${COMPARTMENT_IDS[$COMP_IDX]}"
COMPARTMENT_NAME="${COMPARTMENT_NAMES[$COMP_IDX]}"
echo "  Selected: $COMPARTMENT_NAME"
echo

# ---------- select container instance ----------

echo "Fetching container instances in '$COMPARTMENT_NAME'..."

INSTANCE_NAMES=()
INSTANCE_IDS=()
INSTANCE_STATES=()

while IFS=$'\t' read -r name id state; do
    [[ -z "$name" ]] && continue
    INSTANCE_NAMES+=("$name")
    INSTANCE_IDS+=("$id")
    INSTANCE_STATES+=("$state")
done < <(run_oci container-instances container-instance list \
    --compartment-id "$COMPARTMENT_OCID" \
    --query 'data[].{a:"display-name", b:id, c:"lifecycle-state"}' \
    --output table 2>/dev/null | tail -n +3 | sed 's/  *| */\t/g; s/^|//; s/|$//' || true)

INSTANCE_COUNT=${#INSTANCE_NAMES[@]}

CONTAINER_INSTANCE_OCID=""
if [[ "$INSTANCE_COUNT" -eq 0 ]]; then
    echo "  No container instances found in this compartment."
    read -rp "Paste container instance OCID manually (or leave blank to skip): " CONTAINER_INSTANCE_OCID
else
    echo "Available container instances:"
    for i in "${!INSTANCE_NAMES[@]}"; do
        echo "  [$i] ${INSTANCE_NAMES[$i]} (${INSTANCE_STATES[$i]})"
    done
    echo

    read -rp "Select container instance [0-$((INSTANCE_COUNT - 1))]: " INST_IDX
    CONTAINER_INSTANCE_OCID="${INSTANCE_IDS[$INST_IDX]}"
    echo "  Selected: ${INSTANCE_NAMES[$INST_IDX]}"
fi
echo

# ---------- auth token ----------

echo "An OCIR auth token is needed for Docker registry login."
echo "If you don't have one, create it at:"
echo "  OCI Console -> Identity -> Users -> Your User -> Auth Tokens -> Generate Token"
echo
read -rsp "Paste your auth token (input hidden): " AUTH_TOKEN
echo
echo

# ---------- summary ----------

echo "============================================"
echo "Secrets to push:"
echo "============================================"
echo "  OCI_CLI_USER              = ${USER_OCID:0:40}..."
echo "  OCI_CLI_TENANCY           = ${TENANCY_OCID:0:40}..."
echo "  OCI_CLI_FINGERPRINT       = $FINGERPRINT"
echo "  OCI_CLI_KEY_CONTENT       = ($KEY_BYTES bytes from $KEY_FILE)"
echo "  OCI_CLI_REGION            = $REGION"
echo "  OCI_COMPARTMENT_OCID      = ${COMPARTMENT_OCID:0:40}..."
if [[ -n "$CONTAINER_INSTANCE_OCID" ]]; then
    echo "  OCI_CONTAINER_INSTANCE_OCID = ${CONTAINER_INSTANCE_OCID:0:40}..."
else
    echo "  OCI_CONTAINER_INSTANCE_OCID = (skipped)"
fi
echo "  OCI_AUTH_TOKEN            = (${#AUTH_TOKEN} chars)"
echo "============================================"
echo

if $DRY_RUN; then
    echo "[dry-run] Would push the above secrets. Exiting."
    exit 0
fi

read -rp "Push these secrets to GitHub? [y/N] " CONFIRM
if [[ "$CONFIRM" != [yY] ]]; then
    echo "Aborted."
    exit 0
fi

# ---------- push secrets ----------

echo
echo "Pushing secrets..."

push_secret() {
    local name="$1" value="$2"
    if [[ -z "$value" ]]; then
        echo "  Skipping $name (empty)"
        return
    fi
    echo "$value" | gh secret set "$name"
    echo "  $name"
}

push_secret "OCI_CLI_USER"        "$USER_OCID"
push_secret "OCI_CLI_TENANCY"     "$TENANCY_OCID"
push_secret "OCI_CLI_FINGERPRINT" "$FINGERPRINT"

# Pipe key file directly to gh — avoids holding key material in a shell variable
gh secret set "OCI_CLI_KEY_CONTENT" < "$KEY_FILE"
echo "  OCI_CLI_KEY_CONTENT"

push_secret "OCI_CLI_REGION"      "$REGION"
push_secret "OCI_COMPARTMENT_OCID" "$COMPARTMENT_OCID"
push_secret "OCI_CONTAINER_INSTANCE_OCID" "$CONTAINER_INSTANCE_OCID"
push_secret "OCI_AUTH_TOKEN"      "$AUTH_TOKEN"

echo
echo "Done. Verify with: gh secret list"

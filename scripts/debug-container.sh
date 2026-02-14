#!/usr/bin/env bash
# ============================================================================
# Debug BAMCP container instance — list containers and fetch logs.
#
# Usage:
#   ./scripts/debug-container.sh [--profile PROFILE]
# ============================================================================

set -euo pipefail

PROFILE="DEFAULT"

while [[ $# -gt 0 ]]; do
    case $1 in
        --profile) PROFILE="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--profile PROFILE]"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------- read OCI config ----------

OCI_CONFIG="${OCI_CONFIG_FILE:-$HOME/.oci/config}"

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

# ---------- handle encrypted key ----------

PASSPHRASE_FILE=""
cleanup() {
    [[ -n "$PASSPHRASE_FILE" && -f "$PASSPHRASE_FILE" ]] && rm -f "$PASSPHRASE_FILE"
}
trap cleanup EXIT

OCI_PASS_ARGS=()
if [[ -f "$KEY_FILE" ]] && head -1 "$KEY_FILE" | grep -q "ENCRYPTED"; then
    read -rsp "Enter passphrase: " PASSPHRASE
    echo
    PASSPHRASE_FILE=$(mktemp)
    chmod 600 "$PASSPHRASE_FILE"
    cat > "$PASSPHRASE_FILE" <<EOF
[$PROFILE]
user=$(read_config "user")
tenancy=$TENANCY_OCID
fingerprint=$(read_config "fingerprint")
key_file=$KEY_FILE
region=$REGION
pass_phrase=$PASSPHRASE
EOF
    OCI_PASS_ARGS=(--config-file "$PASSPHRASE_FILE")
    unset PASSPHRASE
fi

run_oci() {
    oci "${OCI_PASS_ARGS[@]}" "$@"
}

# ---------- find bamcp instance ----------

echo "Finding 'bamcp' container instance..."
INSTANCE_OCID=$(run_oci container-instances container-instance list \
    --compartment-id "$TENANCY_OCID" \
    --display-name "bamcp" \
    --query 'data[0].id' \
    --raw-output 2>/dev/null)

if [[ -z "$INSTANCE_OCID" || "$INSTANCE_OCID" == "None" || "$INSTANCE_OCID" == "null" ]]; then
    echo "Error: No 'bamcp' container instance found." >&2
    exit 1
fi

echo "  Instance: $INSTANCE_OCID"
echo

# ---------- get instance details ----------

echo "Instance status:"
run_oci container-instances container-instance get \
    --container-instance-id "$INSTANCE_OCID" \
    --query 'data.{"state":"lifecycle-state","created":"time-created","shape":shape,"ad":"availability-domain"}' \
    --output table 2>/dev/null
echo

# ---------- list containers ----------

echo "Containers:"
run_oci container-instances container-instance get \
    --container-instance-id "$INSTANCE_OCID" \
    --query 'data.containers[].{"name":"display-name","id":id,"state":"lifecycle-state"}' \
    --output table 2>/dev/null
echo

# ---------- get container IDs ----------

CONTAINER_IDS=()
CONTAINER_NAMES=()

while IFS=$'\t' read -r name id; do
    [[ -z "$name" ]] && continue
    CONTAINER_NAMES+=("$name")
    CONTAINER_IDS+=("$id")
done < <(run_oci container-instances container-instance get \
    --container-instance-id "$INSTANCE_OCID" \
    --query 'data.containers[].{"a":"display-name","b":id}' \
    --output table 2>/dev/null | tail -n +3 | sed 's/  *| */\t/g; s/^|//; s/|$//')

# ---------- fetch logs for each container ----------

for i in "${!CONTAINER_NAMES[@]}"; do
    echo "============================================"
    echo "Logs: ${CONTAINER_NAMES[$i]}"
    echo "============================================"
    run_oci container-instances container retrieve-logs \
        --container-id "${CONTAINER_IDS[$i]}" --file - 2>/dev/null | head -100 || echo "  (no logs available)"
    echo
done

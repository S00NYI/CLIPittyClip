#!/bin/bash

# CLIPittyClip_wizard.sh
# Triage entry point for the CLIPittyClip Suite wizard.
# Opens the tool-selection menu (CLIPittyClip / PREPittyPrep / MAPittyMap / PEAKittyPeak),
# walks the user through the chosen tool's wizard, then launches it.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/wizard.sh"

# Help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo ""
    echo "Usage: $0"
    echo ""
    echo "CLIPittyClip Suite Wizard — interactive triage that asks which tool you"
    echo "want to run, walks you through its options, and launches it."
    echo ""
    echo "If you already know which tool you want, you can also call its own"
    echo "wizard directly: CLIPittyClip.sh -w  ·  PREPittyPrep.sh -w  ·"
    echo "MAPittyMap.sh -w  ·  PEAKittyPeak.sh -w"
    echo ""
    exit 0
fi

# Make the four pipeline scripts findable regardless of cwd
export PATH="$SCRIPT_DIR:$PATH"

# Run triage. Dispatcher exports CLIPITTY_EQUIV_CMD on success.
run_wizard_dispatcher || exit $?

if [[ -z "$CLIPITTY_EQUIV_CMD" ]]; then
    echo "[ERROR] Wizard finished without building an equivalent command." >&2
    exit 1
fi

echo ""
echo "─── Launching ─────────────────────────────────────────────────"
echo "$CLIPITTY_EQUIV_CMD"
echo "──────────────────────────────────────────────────────────────"
echo ""
eval "$CLIPITTY_EQUIV_CMD"

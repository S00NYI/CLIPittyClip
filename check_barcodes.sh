#!/bin/bash
# check_barcodes.sh
# Check if a set of barcodes remains unique with N mismatches allowed.
# Usage: check_barcodes.sh -f <barcode_file> -m <mismatches>

# Default values
MISMATCHES=1

function show_usage {
    echo "Usage: $0 -f <barcode_file> [-m <mismatches>]"
    echo "  -f : Path to barcode file (Format: Name<tab>Sequence)"
    echo "  -m : Number of mismatches allowed (Default: 1)"
    exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--file) BARCODE_FILE="$2"; shift 2 ;;
        -m|--mismatches) MISMATCHES="$2"; shift 2 ;;
        -h|--help) show_usage ;;
        *) echo "Unknown option: $1"; show_usage ;;
    esac
done

# Validation
if [[ -z "$BARCODE_FILE" ]]; then
    echo "[ERROR] Missing barcode file."
    show_usage
fi

if [[ ! -f "$BARCODE_FILE" ]]; then
    echo "[ERROR] File not found: $BARCODE_FILE"
    exit 1
fi

echo "  > Checking barcodes in '$BARCODE_FILE' with $MISMATCHES allowed mismatches..."

# Check Uniqueness with AWK
# Recursively generate all N-mismatch variants and check for overlap
COLLISIONS=$(awk -v m="$MISMATCHES" '
    # Recursive function to generate N-substituted variants
    function gen(seq, n, start_pos,    i, new_seq) {
        print seq
        if (n == 0) return
        for (i = start_pos; i <= length(seq); i++) {
            new_seq = substr(seq, 1, i-1) "N" substr(seq, i+1)
            gen(new_seq, n-1, i+1)
        }
    }
    {
        # Generate variants for current barcode (column 2)
        gen($2, m, 1)
    }
' "$BARCODE_FILE" | sort | uniq -d | head -n 5)

if [[ -n "$COLLISIONS" ]]; then
    echo "  > [FAIL] Barcode collisions detected!"
    echo "  > The following sequences could belong to multiple barcodes:"
    echo "$COLLISIONS"
    exit 1
else
    echo "  > All barcodes are unique with $MISMATCHES mismatches."
    exit 0
fi

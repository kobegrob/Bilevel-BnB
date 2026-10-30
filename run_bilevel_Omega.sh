#!/usr/bin/env bash
set -u

EXEC="/home/kobe/Documents/BnB/Bilevel-BnB/bilevel_Omega"
INPUT_DIR="${1:-.}"
OUTPUT_CSV="${2:-results.csv}"

# Check whether the executable exists
if [ ! -x "$EXEC" ]; then
    echo "Error: $EXEC not found or not executable."
    exit 1
fi

# Check whether the input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input folder '$INPUT_DIR' does not exist."
    exit 1
fi

# Create CSV file with header
echo "problem_name,solution_time_seconds,omega_value" > "$OUTPUT_CSV"

shopt -s nullglob

# Run all .knp instances
for file in "$INPUT_DIR"/*.knp; do

    base="$(basename "$file" .knp)"

    echo "Running $base..."

    tmp_out="$(mktemp)"

    # The program expects the problem name without .knp
    # Run from the input directory so the instance can be found
    (
        cd "$INPUT_DIR" || exit 1
        printf '%s\n' "$base" | "$EXEC"
    ) > "$tmp_out" 2>&1

    # Extract omega(1, b1, b2, 0)
    value="$(grep -E 'glb = ' "$tmp_out" | tail -n1 | sed -E 's/.*glb = *([-0-9]+).*/\1/')"

    # Extract computation time
    time_sec="$(grep -E 'Computation time:' "$tmp_out" | tail -n1 | sed -E 's/.*Computation time: *([0-9.]+) seconds.*/\1/')"

    # Fallback if parsing fails
    value="${value:-NA}"
    time_sec="${time_sec:-NA}"

    # Save results
    echo "$base,$time_sec,$value" >> "$OUTPUT_CSV"

    rm -f "$tmp_out"

done

echo "Done. Results written to $OUTPUT_CSV"
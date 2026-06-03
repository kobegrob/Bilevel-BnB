#!/usr/bin/env bash
set -u

EXEC="/home/kobe/Documents/BnB/Bilevel-BnB/bileveldensity"
INPUT_DIR="${1:-.}"
OUTPUT_CSV="${2:-results.csv}"

if [ ! -x "$EXEC" ]; then
  echo "Error: $EXEC not found or not executable."
  exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: input folder '$INPUT_DIR' does not exist."
  exit 1
fi

echo "problem_name,density" > "$OUTPUT_CSV"

shopt -s nullglob
for file in "$INPUT_DIR"/*.knp; do
  base="$(basename "$file" .knp)"

  echo "Running $base..."

  tmp_out="$(mktemp)"

  (
    cd "$INPUT_DIR" || exit 1
    printf '%s\n' "$base" | "$EXEC"
  ) > "$tmp_out" 2>&1

  density="$(grep -E 'Omega density = ' "$tmp_out" | tail -n1 | sed -E 's/.*Omega density = *([0-9.]+).*/\1/')"
  density="${density:-NA}"

  echo "$base,$density" >> "$OUTPUT_CSV"

  rm -f "$tmp_out"
done

echo "Done. Results written to $OUTPUT_CSV"
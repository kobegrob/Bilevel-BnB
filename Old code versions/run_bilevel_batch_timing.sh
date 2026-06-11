#!/usr/bin/env bash
set -u

EXEC="/home/kobe/Documents/BnB/Bilevel-BnB/bilevel_timing"
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

echo "problem_name,solution_time_seconds,greedyupper_seconds,solveP1_seconds,omega_seconds,solution_value,gap_percent,nodes_visited" > "$OUTPUT_CSV"

shopt -s nullglob
for file in "$INPUT_DIR"/*.knp; do
  base="$(basename "$file" .knp)"

  echo "Running $base..."

  tmp_out="$(mktemp)"

  (
    cd "$INPUT_DIR" || exit 1
    printf '%s\n' "$base" | "$EXEC"
  ) > "$tmp_out" 2>&1

  value="$(grep -E 'Value = ' "$tmp_out" | tail -n1 | sed -E 's/.*Value = *([-0-9]+).*/\1/')"

  time_sec="$(grep -E 'Computation time:' "$tmp_out" | tail -n1 | sed -E 's/.*Computation time: *([0-9.]+) seconds.*/\1/')"

  greedy_time="$(sed -nE 's/.*greedyupper\(\) time:[[:space:]]*([0-9.]+).*/\1/p' "$tmp_out" | tail -n1)"

  solvep1_time="$(sed -nE 's/.*solveP1:[[:space:]]*([0-9.]+).*/\1/p' "$tmp_out" | tail -n1)"

  omega_time="$(sed -nE 's/.*omega:[[:space:]]*([0-9.]+).*/\1/p' "$tmp_out" | tail -n1)"

  gap="$(grep -E 'Gap = ' "$tmp_out" | tail -n1 | sed -E 's/.*Gap = *([0-9.]+)%.*/\1/')"

  nodes="$(grep -E 'Nodes Visited:' "$tmp_out" | tail -n1 | sed -E 's/.*Nodes Visited: *([0-9]+).*/\1/')"

  # fallback if parsing fails
  value="${value:-NA}"
  time_sec="${time_sec:-NA}"
  greedy_time="${greedy_time:-NA}"
  solvep1_time="${solvep1_time:-NA}"
  omega_time="${omega_time:-NA}"
  gap="${gap:-NA}"
  nodes="${nodes:-NA}"

  echo "$base,$time_sec,$greedy_time,$solvep1_time,$omega_time,$value,$gap,$nodes" >> "$OUTPUT_CSV"

  rm -f "$tmp_out"
done

echo "Done. Results written to $OUTPUT_CSV"

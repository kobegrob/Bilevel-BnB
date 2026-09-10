#!/usr/bin/env bash

set -e

cd "$HOME/Documents/BnB/Bilevel-BnB"

./run_bilevel_batch.sh n_100 Bilevel_100.csv
./run_bilevel_heuristic.sh n_100 heuristic_100.csv
./run_bilevel_NonLinearLog.sh n_100 NonLinearLog_100.csv
./run_bilevel_StrongerLBNLLog.sh n_100 StongerLBNLLog_100.csv
./run_bilevel_StrongerLBNLQuad.sh n_100 StongerLBNLQuad_100.csv
./run_bilevel_NonLinearQuad.sh n_100 NonLinearQuad_100.csv

echo "All commands finished successfully."
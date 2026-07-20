#!/usr/bin/env bash

set -e

cd "$HOME/Documents/BnB/Bilevel-BnB"

./run_bilevel_NonLinearQuad.sh n_10 NonLinearQuad_10.csv
./run_bilevel_NonLinearQuad.sh n_25 NonLinearQuad_25.csv
./run_bilevel_NonLinearLog.sh n_10 NonLinearLog_10.csv
./run_bilevel_NonLinearLog.sh n_25 NonLinearLog_25.csv
./run_bilevel_StrongerLBNLLog.sh n_25 StongerLBNLLog_25.csv
./run_bilevel_StrongerLBNLQuad.sh n_50_2 StongerLBNLQuad_50.csv
./run_bilevel_StrongerLBNLQuad.sh n_100_2 StongerLBNLQuad_100.csv
./run_bilevel_NonLinearQuad.sh n_50_2 NonLinearQuad_50.csv
./run_bilevel_NonLinearQuad.sh n_100_2 NonLinearQuad_100.csv
./run_bilevel_StrongerLBNLLog.sh n_50_2 StongerLBNLLog_50.csv
./run_bilevel_StrongerLBNLLog.sh n_100_2 StongerLBNLLog_100.csv
./run_bilevel_NonLinearLog.sh n_50_2 NonLinearLog_50.csv
./run_bilevel_NonLinearLog.sh n_100_2 NonLinearLog_100.csv

echo "All commands finished successfully."
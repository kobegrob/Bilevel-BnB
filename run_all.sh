#!/usr/bin/env bash

set -e

cd "$HOME/Documents/BnB/Bilevel-BnB"

./run_bilevel_StrongerLBNLLog.sh n_10 StongerLBNLLog_10.csv
./run_bilevel_StrongerLBNLLog.sh n_25 StongerLBNLLog_25.csv
./run_bilevel_StrongerLBNLLog.sh n_50 StongerLBNLLog_50.csv
./run_bilevel_StrongerLBNLLog.sh n_100 StongerLBNLLog_100.csv
./run_bilevel_StrongerLBNLQuad.sh n_10 StongerLBNLQuad_10.csv
./run_bilevel_StrongerLBNLQuad.sh n_25 StongerLBNLQuad_25.csv
./run_bilevel_StrongerLBNLQuad.sh n_50 StongerLBNLQuad_50.csv
./run_bilevel_StrongerLBNLQuad.sh n_100 StongerLBNLQuad_100.csv
echo "All commands finished successfully."
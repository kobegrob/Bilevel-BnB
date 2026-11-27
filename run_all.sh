#!/usr/bin/env bash

set -e

cd "$HOME/Documents/BnB/Bilevel-BnB"
./run_bilevel_batch.sh n_10 n_10.csv
./run_bilevel_batch.sh n_25 n_25.csv
./run_bilevel_batch.sh n_50 n_50.csv
./run_bilevel_batch.sh n_100 n_100.csv
./run_bilevel_batch.sh BnB_Test BnB_Test.csv
./run_bilevel_NonLinearLog.sh n_10 NonLinearLog_10.csv
./run_bilevel_NonLinearLog.sh n_25 NonLinearLog_25.csv
./run_bilevel_NonLinearLog.sh n_50 NonLinearLog_50.csv
./run_bilevel_NonLinearLog.sh n_100 NonLinearLog_100.csv
./run_bilevel_NonLinearQuad.sh n_10 NonLinearQuad_10.csv
./run_bilevel_NonLinearQuad.sh n_25 NonLinearQuad_25.csv
./run_bilevel_NonLinearQuad.sh n_50 NonLinearQuad_50.csv
./run_bilevel_NonLinearQuad.sh n_100 NonLinearQuad_100.csv
./run_bilevel_StrongerLBNLLog.sh n_10 StongerLBNLLog_10.csv
./run_bilevel_StrongerLBNLLog.sh n_25 StongerLBNLLog_25.csv
./run_bilevel_StrongerLBNLLog.sh n_50 StongerLBNLLog_50.csv
./run_bilevel_StrongerLBNLLog.sh n_100 StongerLBNLLog_100.csv
./run_bilevel_StrongerLBNLQuad.sh n_10 StongerLBNLQuad_10.csv
./run_bilevel_StrongerLBNLQuad.sh n_25 StongerLBNLQuad_25.csv
./run_bilevel_StrongerLBNLQuad.sh n_50 StongerLBNLQuad_50.csv
./run_bilevel_StrongerLBNLQuad.sh n_100 StongerLBNLQuad_100.csv
echo "All commands finished successfully."
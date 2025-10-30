#!/usr/bin/bash

set -euo pipefail

THR="${1:-}" # Check THR is provided
TEST="${2:-}" # Check THR is provided

if [[ -z "$THR" ]];then
    echo "Missing length threshold"
    exit 1
fi

if [[ -z "$TEST" ]];then
    echo "Missing switch: 0 real data, 1 test"
    exit 1
fi

if [[ "$TEST" == 1 ]]; then
    input_folder="../data/test"
    input_results_folder="../results/Test"
else
    input_folder="../data"
    input_results_folder="../results"
fi

#for f in G_i1 G_i2 G_i3 G_i4 G_i6 G_i7
for f in G_i1 G_i3 G_i4
do
    echo $f
    python filter_by_qname_sam.py -i ${input_folder}/${f}_${THR}.sam -o ${input_folder}/${f}_${THR}_compatible.sam --qnames ${input_results_folder}/index_${f}_${THR}_compatible.csv --keep --verbose
    python filter_by_qname_sam.py -i ${input_folder}/${f}_${THR}.sam -o ${input_folder}/${f}_${THR}_noncompatible.sam --qnames ${input_results_folder}/index_${f}_${THR}_noncompatible.csv --keep --verbose
done

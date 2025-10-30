#!/usr/bin/env bash

set -euo pipefail

THR="${1:-}"
TEST="${2:-}"

if [[ -z "$THR" ]];then
    echo "Missing length threshold"
    exit 1
fi
if [[ -z "$TEST" ]];then
    echo "Missing switch"
    exit 1
fi

echo "Convert to vector"
./convert_to_vector.sh ${THR} ${TEST}
echo "Get summary"
./get_summary_table.sh ${THR} ${TEST}
echo "Plot mutations"
./plot_mutations.py ${THR} ${TEST}
echo "Classification reads"
./classification_reads.sh ${THR} ${TEST}
./stats_classification.py ${THR} ${TEST}
echo "Split data"
./split_index.sh ${THR} ${TEST}
./filter_by_qname_sam.sh ${THR} ${TEST}

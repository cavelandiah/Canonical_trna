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

./convert_to_vector.sh ${THR} ${TEST}
./get_summary_table.sh ${THR} ${TEST}
./plot_mutations.py ${THR} ${TEST}
./classification_reads.sh ${THR} ${TEST}
./stats_classification.py ${THR} ${TEST}
./split_index.sh ${THR} ${TEST}
./filter_by_qname_sam.sh ${THR} ${TEST}

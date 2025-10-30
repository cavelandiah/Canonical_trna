#!/usr/bin/env bash

set -euo pipefail

THR="${1:-}" # Check THR is provided

if [[ -z "$THR" ]];then
    echo "Missing length threshold"
    exit 1
fi

./convert_to_vector.sh ${THR}
./get_summary_table.sh ${THR}
./plot_mutations.py ${THR}
./classification_reads.sh ${THR}
./stats_classification.py ${THR} "1"
./split_index.sh ${THR}
./filter_by_qname_sam.sh ${THR} "1"

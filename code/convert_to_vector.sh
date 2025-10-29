#!/usr/bin/env bash

set -euo pipefail

THR="${1:-}" # Check THR is provided

if [[ -z "$THR" ]];then
    echo "Missing length threshold"
    exit 1
fi

#for f in G_i1 G_i2 G_i3 G_i4 G_i6 G_i7
for f in G_i1 G_i3 G_i4
do
    echo $f
    ./convert_to_vector.py $f ${THR} 1
done

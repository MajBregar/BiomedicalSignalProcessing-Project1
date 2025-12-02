#!/bin/bash

in_dir="$1"      # input dataset directory
out_dir="$2"     # output directory for .mat/.hea/.dat

# ensure absolute paths
in_dir="$(realpath "$in_dir")"
out_dir="$(realpath "$out_dir")"

mkdir -p "$out_dir"

# change into input directory so wfdb2mat recognizes record names
cd "$in_dir" || exit 1

for f in *.hea; do
    rec="${f%.*}"    # record name without extension
    echo "$rec"

    # wfdb2mat writes: recm.mat, recm.info, rec.hea, rec.dat
    wfdb2mat -r "$rec"

    # move all generated files for this record into output dir
    mv "${rec}"* "$out_dir"/
done

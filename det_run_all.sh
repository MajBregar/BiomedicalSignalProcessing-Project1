#!/bin/bash

trap "echo 'Interrupted'; exit 1" INT

#clear all asc and qrs files
rm -f converted_ltst/*.asc converted_ltst/*.qrs

#rerun for every file
record_dir_path="converted_ltst"
for f in ${record_dir_path}/*.dat
do
    f=$(basename "$f")
    record_id="${f%.*}"
    record_path="${record_dir_path}/${record_id}"

    echo "Processing Record: $record_path"
    samp_rate=$(sampfreq "$record_path")

    cd "$(dirname "$0")/src" || exit 1
    octave --no-gui --quiet --eval "Run_Detector(\"../converted_ltst/$record_id\", $samp_rate, false)"
    cd ..
done


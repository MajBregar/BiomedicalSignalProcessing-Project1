#!/bin/bash

trap "echo 'Interrupted'; exit 1" INT

#clear all asc and qrs files
rm -f converted_mitbih/*.asc converted_mitbih/*.qrs

#rerun for every file
record_dir_path="converted_mitbih"
for f in ${record_dir_path}/*.dat
do
    f=$(basename "$f")
    record_id="${f%.*}"
    record_path="${record_dir_path}/${record_id}"

    echo "Processing Record: $record_path"
    samp_rate=$(sampfreq "$record_path")

    cd "$(dirname "$0")/src" || exit 1
    octave --no-gui --quiet --eval "Run_Detector(\"../converted_mitbih/$record_id\", $samp_rate, false)"
    cd ..
done

#Now you can copy average Se and +P from results.txt
#Columns in results.txt that are of your interest are following:
# (note that the detector does not distinguish between different
# types of heart-beats, it just needs to detect them - so it should
# detect all N, V, F, etc. heart-beats as N, hence, FP=Nn+Vn+Fn,
# for heartbeat types see:
# - (Nn+Vn+Fn) = true positive (N, V, or F detected as N)
# - (No+Vo+Fo) = false negative (failed to detect N, V, or F)
# - On = false positive (heartbeat was detected where there is none)
# Now you can calculate Se, +P using formulas from lectures.
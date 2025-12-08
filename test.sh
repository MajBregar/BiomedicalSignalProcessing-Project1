#!/bin/bash

trap "echo 'Interrupted'; exit 1" INT

record_dir_path="converted_ltst"
counter=0

for f in "${record_dir_path}"/*.dat
do
    f=$(basename "$f")
    record_id="${f%.*}"

    echo "$record_id"

    ((counter++))
done
echo "Number of Records: $counter"
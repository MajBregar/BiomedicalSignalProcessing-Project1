#!/bin/bash

record_dir_path="converted_mitbih"
qrs_output_dir="mitbih_original_qrs"

rm -f eval1.txt eval2.txt
mkdir -p "$qrs_output_dir"

cd "$record_dir_path"

for f in *.asc
do
    record_id="${f%.*}"

    wrann -r "$record_id" -a qrs < "$record_id.asc"

    mv "${record_id}.qrs" "../$qrs_output_dir/"
done



for f in *.asc
do
    record_id="${f%.*}"
    cp "../$qrs_output_dir/${record_id}.qrs" .
    bxb -r "$record_id" -a atr qrs -l "../eval1.txt" "../eval2.txt"
    rm -f "${record_id}.qrs"
done

cd ..
sumstats eval1.txt eval2.txt > eval_benchmark.txt
rm -f eval1.txt eval2.txt

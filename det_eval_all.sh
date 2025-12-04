#!/bin/bash

rm -f eval1.txt eval2.txt

record_dir_path="converted_mitbih"
cd converted_mitbih

for f in *.asc
do
    f=$(basename "$f")
    record_id="${f%.*}"
    
    wrann -r $record_id -a qrs < "$record_id.asc"
    bxb -r $record_id -a atr qrs -l ../eval1.txt ../eval2.txt
done

cd ..
sumstats eval1.txt eval2.txt > eval_benchmark.txt
rm -f eval1.txt eval2.txt

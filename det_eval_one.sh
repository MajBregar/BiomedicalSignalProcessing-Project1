#!/bin/bash

record_id="${1:-100}"
database_dir="${2:-converted_mitbih}"

#clear statistics
rm -f eval1.txt eval2.txt

#evaluate record
cd "./$database_dir"
wrann -r $record_id -a qrs < "$record_id.asc"
bxb -r $record_id -a atr qrs -l ../eval1.txt ../eval2.txt


#final stats
cd ..
sumstats eval1.txt eval2.txt > eval_benchmark.txt
rm -f eval1.txt eval2.txt

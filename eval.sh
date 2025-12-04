#!/bin/bash

record_id="100"

#clear statistics
rm -f eval1.txt eval2.txt

#evaluate record
cd converted_mitbih
wrann -r $record_id -a qrs < "$record_id.asc"
bxb -r $record_id -a atr qrs -l ../eval1.txt ../eval2.txt


#final stats
cd ..
sumstats eval1.txt eval2.txt > eval_benchmark.txt
rm -f eval1.txt eval2.txt

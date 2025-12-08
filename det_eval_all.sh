#!/bin/bash

rm -f eval1.txt eval2.txt

record_dir_path="converted_ltst"
cd converted_ltst

for f in *.asc
do
    f=$(basename "$f")
    record_id="${f%.*}"
    
    echo $record_id
    wrann -r $record_id -a qrs < "$record_id.asc"
    bxb -r $record_id -a atr qrs -l ../eval1.txt ../eval2.txt
done

cd ..
sumstats eval1.txt eval2.txt > eval_benchmark.txt
rm -f eval1.txt eval2.txt

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
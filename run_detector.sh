#!/bin/bash

rm eval1.txt
rm eval2.txt
FILES=/path/to/directory/with/records/*.dat
for f in $FILES
do
 f=$(basename $f)
 f=${f%.*}

 echo $f
 wfdb2mat -r $f #convert to Matlab format
done

#Run algorithm in Matlab. Output should be annotations in text files
#with WFDB annotator structure. See Matlab frame on the web classroom.
for f in $FILES
do
 f=$(basename $f)
 f=${f%.*}

 echo $f
 #convert text annotator to WFDB format
 wrann -r $f -a qrs < $f".asc"
 #evaluate using reference annotations atr and your .asc files
 bxb -r $f -a atr qrs -l eval1.txt eval2.txt
done
sumstats eval1.txt eval2.txt > results.txt #final statistics
#Now you can copy average Se and +P from results.txt
#Columns in results.txt that are of your interest are following:
# (note that the detector does not distinguish between different
# types of heart-beats, it just needs to detect them - so it should
# detect all N, V, F, etc. heart-beats as N, hence, FP=Nn+Vn+Fn,
# for heartbeat types see:
# https://archive.physionet.org/physiobank/annotations.shtml)
# - (Nn+Vn+Fn) = true positive (N, V, or F detected as N)
# - (No+Vo+Fo) = false negative (failed to detect N, V, or F)
# - On = false positive (heartbeat was detected where there is none)
# Now you can calculate Se, +P using formulas from lectures.
#!/bin/sh

#PATH = ""
# topics="20 30 50 80 100 150 200 250 300 400"
# iterations="0 1 2 3 4 5 6 7 8 9"
topics="20 30"
iterations="0 1"
for a in $topics
do
	echo "$a"
	for i in $iterations
	do 
		echo "$i"
		./sctm ../input/unigram/trainFolder$i/abagf.AT.txt ../input/unigram/trainFolder$i/cbagf.AT.txt ../output/sctm$a/$i $a sctm 0
		./sctm ../input/unigram/testFolder$i/abagf.AT.txt ../input/unigram/testFolder$i/cbagf.AT.txt ../output/sctm$a/$i $a sctm 1 > ../output/sctm$a/$i/output_$a_$i.txt
	done
done
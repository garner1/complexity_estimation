#!/usr/bin/env bash

# Compute the smoothed Good-Toulmin estimator.
# Based on arXiv:1511.07428.(http://www.pnas.org/content/113/47/13283.full)
# see also https://arxiv.org/abs/1707.03854

file=$1				# full path to q10_chr-loc-strand-umi-pcr input file
step=$2				# step size on the x-axis of the accumulation plot
sd_samples=$3			# number of samples for the smoothing distro
field=$4			# field number where to find pcr duplicates count
extension=$5			# how many times the size of the initial sample you want to sequence now

#EXTRACT THE HISTOGRAM OF HISTOGRAM 
cut -f"$field" $file | LC_ALL=C sort -n | LC_ALL=C uniq -c | awk '{print $1"\t"$2}' | LC_ALL=C sort -k2,2 > phi_i.tsv
max=$(cat phi_i.tsv | datamash max 2)
seq $max | LC_ALL=C sort > all_i.tsv
join -v 2 -1 2 -2 1 phi_i.tsv all_i.tsv|awk '{print 0"\t"$1}' > phi_missingI.tsv
cat phi_i.tsv phi_missingI.tsv | LC_ALL=C sort -k2,2n | cut -f1 > prevalences.tsv # INTEGRATE WITH THE FULL LIST OF COUNTS FOR I

n=$(cat $file | datamash sum "$field")	# the number of samples for the initial experiment: THE SUM OF THE DUPLICATES
#t=$(echo $n | awk '{print int(log($1)/log(10))}') # t defined as the floor of the natural log of n: SET THE MAXIMUM FOLD-CHANGE FOR THE ACCUMULATION CURVE
t=$extension

python runSGT.py $n $t $sd_samples $step # THE OUTPUT IS accumulatation_curve.tsv

cat accumulation_curve.tsv | tr -d '(' | tr -d ')' | tr -d ' ' |tr ',' '\t' | awk -v offset=$n '{print $1"\t"offset+offset*$1"\t"$2}' > aux && mv aux accumulation_curve.tsv

title=$(echo $file|rev|cut -d'/' -f1|rev)
titlepdf=$(echo $file|rev|cut -d'/' -f-3|rev|tr '/' '_')
cp accumulation_curve.tsv "$titlepdf"
cat "$titlepdf" | gnuplot -p -e "set term pdf;set output '$titlepdf.pdf';set title '$title';set xlabel 'number of reads';set ylabel 'complexity' ;plot '/dev/stdin' u 2:3"
xpdf $titlepdf.pdf
rm phi_i.tsv all_i.tsv phi_missingI.tsv prevalences.tsv # cleaning

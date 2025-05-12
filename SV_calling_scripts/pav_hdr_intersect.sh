#!/bin/bash


#Bedtools intersect to remove pav variants within HDRs
#bedtools start is 0 based, so subtract 1 from POS


#Convert filtered PAV variants to bed format (combine chromosome and strain so Bedtools treats each strain individually)
awk 'BEGIN {OFS="\t"} NR == 1 {print $1, "start", "end", "line"} NR > 1 {print $1"_"$4, $2-1, $2, NR - 1}' vfilt_longer.tsv > vfilt_longer.bed

#Convert HDRs to bed format (also combine chromosome and strain)
awk 'BEGIN {OFS="\t"} NR == 1 {print $0} NR > 1 {print $1"_"$4, $2, $3}' hdrs.tsv > hdrs.bed

#Run bedtools intersect 
bedtools intersect -a vfilt_longer.bed -b hdrs.bed -v > vfilt_longer_no_hdrs.bed
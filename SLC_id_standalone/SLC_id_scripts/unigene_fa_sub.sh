#!/usr/bin/env bash

### not working

#cd /home/shanedenecke/Dropbox/wp5_midgut_uptake/ABC_phylo/ABC_id/Nv/annotation/proteome
## first argument fasta file
## 2nd argument file with list of genes


## sort fasta according to length
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  $1 |\
    awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
    sort -k1,1rn | cut -f 2- | tr "\t" "\n" > sorted.fa

#extract longest fasta from each record
touch unigene.fa


IFS=$'\n'; 
for next in $(cat $2)
do 
  grep -F -A 1 -m 1 ${next} sorted.fa >> unigene.fa
done

cat unigene.fa

rm unigene.fa
rm sorted.fa

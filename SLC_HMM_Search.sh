#!/usr/bin/env bash

## argument 1= database directory
## argument 2= new proteome to search
## argument 3= New_database folder to create

#A=/data2/shane/Documents/SLC_id/HomSap_Database
#B=/data2/shane/Documents/SLC_id/proteomes/HUMAN_DenPon_unigene.faa
#C=/data2/shane/Documents/SLC_id/Human_search/HUMAN_DenPon_unigene.faa


## reset bash and create new output directory
source ~/.bashrc
mkdir $3
cd $3
##mkdir $C 
#cd $C

echo 'Performing HMM search'
## search HMM profiles against target proteome using hmm search
rm -rf ./hmm_outputs
mkdir ./hmm_outputs
for i in $1/hmm_profiles/*; do
#for i in $A/hmm_profiles/*; do
  base=$(echo $(basename $i))
  hmmsearch --notextw -E 20 $i $2 > ./hmm_outputs/$base.hmmoutput
  #hmmsearch --notextw -E 20 $i $B > ./hmm_outputs/$base.hmmoutput
done


## parse hmm outputs into usable format
rm -rf ./hmm_clean
mkdir ./hmm_clean
for i in ./hmm_outputs/*; do
  base=$(echo $(basename $i))
  cat $i | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > ./hmm_clean/$base.table
done


## extract fasta sequences from All HMM hits
rm -rf ./SLC_fa
mkdir ./SLC_fa
for i in ./hmm_clean/*.table; do
  base=$(echo $(basename $i))
  cut -f 9 $i | sed 's/\s+//g'| /data2/shane/Applications/custom/unigene_fa_sub.sh $2 - > ./SLC_fa/$base'.fa'
  #cut -f 9 $i | sed 's/\s+//g' | /data2/shane/Applications/custom/unigene_fa_sub.sh $B - > ./SLC_fa/$base'.fa'
done
find ./SLC_fa/* -size 0 -delete

##perform blast with HMM hits as queries and the source genomes SLC_mark.fa proteome as a target 
echo 'Blast away'
rm -rf ./recip_blast
mkdir ./recip_blast
for i in ./SLC_fa/*.fa; do
  base=$(echo $(basename $i))
  #blastp -query $i -db $A/reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 > ./recip_blast/$base'_blast.tsv'
  blastp -query $i -db $1/reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 6 -max_hsps 1 -num_threads 24 > ./recip_blast/$base'_blast.tsv'
done
find ./recip_blast/* -size 0 -delete


## Run R script to output table 
mkdir ./prelim_summary
#Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_Family_Sort.R $A'/SLC_source_dict.csv' >  ./prelim_summary/Family_sort_preliminary.csv
Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_Family_Sort.R $1'/SLC_source_dict.csv' >  ./prelim_summary/Family_sort_preliminary.csv



## Filter based on lengths of human SLC gene and on TM domains
mkdir length_analysis

### make fasta from reciprocal blast results
cut -d ',' -f 1 ./prelim_summary/Family_sort_preliminary.csv | /data2/shane/Applications/custom/unigene_fa_sub.sh $2 - > ./length_analysis/preliminary_SLC.fa
#cut -d ',' -f 1 ./prelim_summary/Family_sort_preliminary.csv | /data2/shane/Applications/custom/unigene_fa_sub.sh $B - > ./length_analysis/preliminary_SLC.fa
#B=/data2/shane/Documents/SLC_id/proteomes/AcrEch_unigene.faa
### Create table of each gene with corresponding number of TM domains

cut -d ',' -f 2 ./prelim_summary/Family_sort_preliminary.csv | sed '1d' > ./length_analysis/length_families.txt ### why remove last line
#grep ">" ./length_analysis/preliminary_SLC_TMHMM.fa | sed 's/>//g'

awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ./length_analysis/preliminary_SLC.fa > ./length_analysis/names_lengths.txt
grep ">" ./length_analysis/names_lengths.txt | perl -pe  's/^>(.+$)/$1/;'| cut -d ' ' -f 1  > ./length_analysis/all_proteins.txt
grep -E "^[0-9]" ./length_analysis/names_lengths.txt > ./length_analysis/all_lengths.txt
paste -d',' ./length_analysis/all_proteins.txt ./length_analysis/all_lengths.txt ./length_analysis/length_families.txt > ./length_analysis/gene_lengths.txt
Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_length_filter.R > ./length_analysis/total_slc_table.csv




#### Produce final fasta file
rm -f ./length_analysis/SLC_final.faa
for i in $(cat ./length_analysis/total_slc_table.csv | cut -d ',' -f 1)
do
grep -A 1 $i ./length_analysis/preliminary_SLC.fa >> ./length_analysis/SLC_final.faa
done



## final output
mkdir final_output
cp ./length_analysis/total_slc_table.csv ./final_output/total_slc_table.csv 
cp ./length_analysis/SLC_final.faa ./final_output/SLC_final.faa
Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_dictionary_format.R > ./final_output/SLC_final_output.csv




#/data2/shane/Applications/custom/tmhmm_filter.sh ./length_analysis/preliminary_SLC.fa 0 > ./length_analysis/preliminary_SLC_TMM_table.txt

##Re generate dictionary filtered for TMM values
#rm -f ./length_analysis/SLC_TMM_filter_codes.csv
#for i in $(cat ./length_analysis/preliminary_SLC_TMM_table.txt | cut -f 1)
#do
#grep $i ./prelim_summary/Family_sort_preliminary.csv >> ./length_analysis/SLC_TMM_filter_codes.csv
#done
#sed  -i '1i code,name' ./length_analysis/SLC_TMM_filter_codes.csv

## regenerate fasta file
#cut -d ',' -f 1 ./length_analysis/SLC_TMM_filter_codes.csv | sed '1d' | /data2/shane/Applications/custom/unigene_fa_sub.sh ./length_analysis/preliminary_SLC.fa - > ./length_analysis/TMM_Filter_SLC.fa


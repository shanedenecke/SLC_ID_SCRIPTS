#!/usr/bin/env bash

## argument 1= database directory
## argument 2= new proteome to search
## argument 3= New_database folder to create

#B=~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa
#A=~/Documents/SLC_id/Human_HMM_SLC
#C=~/Documents/SLC_id/Drosophila_Database/Hs_to_Dm_Search


## reset bash and create new output directory
source ~/.bashrc
mkdir $3
cd $3
##mkdir $C ### not making proper directory
#cd $C


## search HMM profiles against target proteome using hmm search
rm -rf ./hmm_outputs
mkdir ./hmm_outputs
rm $1/hmm_profile/*hmmoutput
for i in $1/hmm_profiles/*; do
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
  cut -f 9 $i | sed 's/\s+//g'| ~/Applications/custom/unigene_fa_sub.sh $2 - > ./SLC_fa/$base'.fa'
  #cut -f 9 $i | sed 's/\s+//g' | ~/Applications/custom/unigene_fa_sub.sh $B - > ./SLC_fa/$base'.fa'
done

##perform blast with HMM hits as queries and the source genomes SLC_mark.fa proteome as a target 
rm -rf ./recip_blast
mkdir ./recip_blast
for i in ./SLC_fa/*.fa; do
  base=$(echo $(basename $i))
  #blastp -query $i -db $A/reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 > ./recip_blast/$base'_blast.tsv'
  blastp -query $i -db $1/reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 > ./recip_blast/$base'_blast.tsv'
done
find ./recip_blast/* -size 0 -delete


## Run R script to output table 
mkdir ./prelim_summary
#Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_Family_Sort.R $A'/SLC_source_dict.csv' >  ./prelim_summary/Family_sort_preliminary.csv
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_Family_Sort.R $1'/SLC_source_dict.csv' >  ./prelim_summary/Family_sort_preliminary.csv



## Filter based on lengths of human SLC geen family
mkdir length_analysis
cut -d ',' -f 1 ./prelim_summary/Family_sort_preliminary.csv | ~/Applications/custom/unigene_fa_sub.sh $2 - > ./length_analysis/preliminary_SLC.fa
#cut -d ',' -f 1 ./prelim_summary/Family_sort_preliminary.csv | ~/Applications/custom/unigene_fa_sub.sh $B - > ./length_analysis/preliminary_SLC.fa
cut -d ',' -f 2 ./prelim_summary/Family_sort_preliminary.csv | sed '1d' > ./length_analysis/length_families.txt
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ./length_analysis/preliminary_SLC.fa > ./length_analysis/names_lengths.txt
grep ">" ./length_analysis/names_lengths.txt | perl -pe  's/^>(.+$)/$1/;'| cut -d ' ' -f 1  > ./length_analysis/all_proteins.txt
grep -E "^[0-9]" ./length_analysis/names_lengths.txt > ./length_analysis/all_lengths.txt
paste -d',' ./length_analysis/all_proteins.txt ./length_analysis/all_lengths.txt ./length_analysis/length_families.txt > ./length_analysis/gene_lengths.txt
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_length_filter.R > ./length_analysis/total_slc_table.csv

## final output
mkdir final_output
mv ./length_analysis/total_slc_table.csv ./final_output/total_slc_table.csv 
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_table_format.R > ./final_output/SLC_final_output.csv

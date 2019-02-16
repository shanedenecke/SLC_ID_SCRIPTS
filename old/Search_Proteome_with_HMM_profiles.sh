#!/usr/bin/env bash

source ~/.bashrc


## argument 1= database directory
## argument 2= new proteome to search
## argument 3= New_database folder to create


#A=$f
#B=$i
#C=~/Documents/SLC_id/ITERATIVE_Bg


cd $1
#cd $A
rm -R ./*FIND
rm -R $3/*


## search HMM profiles against Helicoverpa proteome
mkdir ./hmm_outputs_FIND
rm ./hmm_outputs_FIND/*
rm ./hmm_profiles/*.hmmoutput
for i in ./hmm_profiles/*
do
  hmmsearch --notextw -E 20 $i $2 > $i.hmmoutput
  #hmmsearch --notextw -E 20 $i $B > $i.hmmoutput
done
mv ./hmm_profiles/*.hmmoutput ./hmm_outputs_FIND


## parse hmm outputs
mkdir ./hmm_clean_FIND
for i in ./hmm_outputs_FIND/*
do
cat $i | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > $i.table #| perl -pe 's/(^.+) type=protein.+$/$1/;' > $i.table
done
mv ./hmm_outputs_FIND/*.table ./hmm_clean_FIND




## extract fasta sequences
mkdir SLC_fa_FIND
rm  ./SLC_fa_FIND/*
rm ./hmm_clean_FIND/*Hs.fa
for i in ./hmm_clean_FIND/*.table
do
  cut -f 9 $i | sed 's/\s+//g' > temp.txt
  ~/Applications/custom/unigene_fa_sub.sh $2 ./temp.txt > $i'.fa'
  #~/Applications/custom/unigene_fa_sub.sh $B ./temp.txt > $i'.fa'
  rm temp.txt
done
mv ./hmm_clean_FIND/*.fa ./SLC_fa_FIND/


##perform blast
mkdir recip_blast_FIND
rm ./recip_blast_FIND/*blast.tsv
rm ./SLC_fa_FIND/*blast.tsv
for i in ./SLC_fa_FIND/*.fa
do
  blastp -query $i -db ./reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 > $i'_blast.tsv'
done
mv ./SLC_fa_FIND/*blast.tsv ./recip_blast_FIND/
find ./recip_blast_FIND/* -size 0 -delete


## Run R script to output table 
mkdir prelim_summary_FIND
#Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_id_draft_summary.R $2 >  ./prelim_summary/preliminary_summary.csv
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_id_draft_summary.R $(pwd)'/SLC_dict.csv' >  ./prelim_summary_FIND/preliminary_summary.csv



## Get length table for original fasta
mkdir length_analysis_FIND
cd length_analysis_FIND
cut -d ',' -f 1 ../prelim_summary_FIND/preliminary_summary.csv | ~/Applications/custom/unigene_fa_sub.sh $2 - > ./preliminary_SLC.fa
#cut -d ',' -f 1 ../prelim_summary_FIND/preliminary_summary.csv | ~/Applications/custom/unigene_fa_sub.sh $B - > ./preliminary_SLC.fa
cut -d ',' -f 2 ../prelim_summary_FIND/preliminary_summary.csv | sed '1d' > ./length_families.txt
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ./preliminary_SLC.fa > ./names_lengths.txt
grep ">" ./names_lengths.txt | perl -pe  's/^>(.+$)/$1/;'| cut -d ' ' -f 1  > ./all_proteins.txt
grep -E "^[0-9]" ./names_lengths.txt > ./all_lengths.txt
paste -d',' ./all_proteins.txt ./all_lengths.txt ./length_families.txt > gene_lengths.txt
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_length_compare_general.R > final_slc_list.csv
cd ..

## final transport
mkdir dictionary_create_FIND
mv ./length_analysis_FIND/final_slc_list.csv ./dictionary_create_FIND/FINAL_LIST.csv
cut -f 1 ./dictionary_create_FIND/FINAL_LIST.csv| cut -d ' ' -f 1 > ./dictionary_create_FIND/final_code.csv
cut -d ',' -f 2 ./dictionary_create_FIND/FINAL_LIST.csv > ./dictionary_create_FIND/final_slcs.csv
paste -d',' ./dictionary_create_FIND/final_code.csv ./dictionary_create_FIND/final_slcs.csv > ./dictionary_create_FIND/reduced.csv
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_dictionary_format.R > ./dictionary_create_FIND/SLC_dict.csv
cd ..
mkdir $3 

cd $1
for i in *FIND
do
mv $i $3
done

cd $3
cp ./dictionary_create_FIND/SLC_dict.csv ./SLC_dict.csv



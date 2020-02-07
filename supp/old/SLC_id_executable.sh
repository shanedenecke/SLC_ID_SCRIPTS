#!/usr/bin/env bash

source ~/.bashrc

cd /home/shanedenecke/Dropbox/wp7_prodrug/SLC_id


#B=~/Documents/omics_data/Drosophila_melanogaster/unigene/Dm_unigene.faa
#C=~/Dropbox/wp7_prodrug/SLC_id/general_reference/Dm_update/DM_SLC_dict_post_HsHMM.csv
#D=~/Documents/omics_data/Helicoverpa_armigera/unigene/Ha_unigene.faa

### Arguments as follow.
## 1) Absolute path of your starting species unigene fasta file
## 2) Absolute path for your starting species dictionary
## 3) Absolute path for your target species unigene fasta file

newdir=$(basename $3)
#newdir=$(basename $D)
mkdir $newdir
cd $newdir


## rename fasta files with tag for SLC family
mkdir ./reference_proteome
~/Applications/custom/fasta_rename.py $1 $2 > ./reference_proteome/proteome_SLC_mark.fa

##Make blast dbs for target reference proteome
makeblastdb -in ./reference_proteome/proteome_SLC_mark.fa -parse_seqids -dbtype prot

## extract list of SLC names to be used in future steps
mkdir ./list
cat $2 | tr -s ' ' ',' | csvcut -c name  > ./list/SLC_names_final.txt

## extract fasta of gene lists from Dmel
~/Applications/custom/unigene_fa_sub.sh ./reference_proteome/proteome_SLC_mark.fa  ./list/SLC_names_final.txt > ./list/SLC_genes.fa


## divide each family into independent fasta
mkdir ./family_fasta 
rm ./family_fasta/*
IFS=$'\n'; 
for next in $(cat /home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/general_reference/SLC_families.txt)
do 
  grep -A 1 ${next}"_" ./list/SLC_genes.fa | sed '/--/d' > ./family_fasta/${next}.fa
done

## align all sequences
mkdir ./muscle_alignments
for i in ./family_fasta/*
do
~/Applications/muscle3.8.31_i86linux64 -in $i -out $i.aln
/home/shanedenecke/Applications/trimal-trimAl/source/trimal -in $i.aln -out $i.trimmed
done
mv ./family_fasta/*.trimmed ./muscle_alignments
rm ./family_fasta/*.aln

## build HMM profiles for each family
mkdir ./hmm_profiles
for i in ./muscle_alignments/*
do
hmmbuild $i.hmm $i
done
mv ./muscle_alignments/*.hmm ./hmm_profiles


## search HMM profiles against Helicoverpa proteome
mkdir ./hmm_outputs
rm ./hmm_outputs/*
rm ./hmm_profiles/*.hmmoutput
for i in ./hmm_profiles/*
do
  hmmsearch --notextw -E 20 $i $3 > $i.hmmoutput
done
mv ./hmm_profiles/*.hmmoutput ./hmm_outputs


## parse hmm outputs
mkdir ./hmm_clean
for i in ./hmm_outputs/*
do
cat $i | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > $i.table #| perl -pe 's/(^.+) type=protein.+$/$1/;' > $i.table
done
mv ./hmm_outputs/*.table ./hmm_clean




## extract fasta sequences
mkdir SLC_fa
rm  ./SLC_fa/*
rm ./hmm_clean/*Hs.fa
for i in ./hmm_clean/*.table
do
  cut -f 9 $i | sed 's/\s+//g' > temp.txt
  ~/Applications/custom/unigene_fa_sub.sh $3 ./temp.txt > $i'.fa'
  rm temp.txt
done
mv ./hmm_clean/*.fa ./SLC_fa/


##perform blast
mkdir recip_blast
rm ./recip_blast/*blast.tsv
rm ./SLC_fa/*blast.tsv
for i in ./SLC_fa/*.fa
do
  blastp -query $i -db ./reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 > $i'_blast.tsv'
done
mv ./SLC_fa/*blast.tsv ./recip_blast/
find ./recip_blast/* -size 0 -delete


## Run R script to output table 
mkdir prelim_summary
Rscript ~/Dropbox/wp7_prodrug/SLC_id/SLC_id_scripts/SLC_id_draft_summary.R $2 >  ./prelim_summary/preliminary_summary.csv


## Get length table for original fasta
mkdir length_analysis
cd length_analysis

cut -d ',' -f 1 ../prelim_summary/preliminary_summary.csv | ~/Applications/custom/unigene_fa_sub.sh $3 - > ./preliminary_SLC.fa
cut -d ',' -f 2 ../prelim_summary/preliminary_summary.csv | sed '1d' > ./length_families.txt
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ./preliminary_SLC.fa > ./names_lengths.txt
grep ">" ./names_lengths.txt | perl -pe  's/^>(.+$)/$1/;' > ./all_proteins.txt
grep -E "^[0-9]" ./names_lengths.txt > ./all_lengths.txt

paste -d',' ./all_proteins.txt ./all_lengths.txt ./length_families.txt > gene_lengths.txt
Rscript ~/Dropbox/wp7_prodrug/SLC_id/SLC_id_scripts/SLC_length_compare_general.R > final_slc_list.csv
cd ..

## final transport
mkdir FINAL_ANALYSIS
mv ./length_analysis/final_slc_list.csv ./FINAL_ANALYSIS/FINAL_LIST.csv


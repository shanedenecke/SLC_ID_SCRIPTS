#!/usr/bin/env bash

### ###################### Create HMM profiles for each SLC family in species ###################
## currently stuck on server. Error is that /lib/x86_64-linux-gnu/libm.so.6: version `GLIBC_2.27' required by trimAL

## argument 1= Predicted unigene protein dataset
## argument 2= SLC dictionary with codes and names
## argument 3= Name of database


### This script will create a HMM database for each SLC family in a given species. The inputs are a uniprotein set, 
### an SLC dictionary (gene codes and SLC names) and a new directory that will contain the database


## use these for debugging

#A=~/Documents/SLC_id/proteomes/model/Dm_unigene.faa
#B=$d
#C=~/Documents/SLC_id/Dm_Final_Database
################################################################################


## reset bash
source ~/.bashrc

## Create and set working directory
mkdir $3
cd $3

#mkdir $C
#cd $C


## rename fasta files with tag for SLC family
mkdir reference_proteome
~/Applications/custom/fasta_rename.py $1 $2 > ./reference_proteome/proteome_SLC_mark.fa
#~/Applications/custom/fasta_rename.py $A $B > ./reference_proteome/proteome_SLC_mark.fa

#~/Applications/custom/fasta_rename.py ~/Documents/omics_data/Homo_sapiens/unigene/Hs_unigene.fa ~/Dropbox/wp7_prodrug/SLC_id/general_reference/Hs_SLC_dict.csv > ./reference_proteome/proteome_SLC_mark.fa

##Make blast dbs for target reference proteome
makeblastdb -in ./reference_proteome/proteome_SLC_mark.fa -parse_seqids -dbtype prot

## extract list of SLC names to be used in future steps
mkdir ./list
cat $2 | tr -s ' ' ',' | csvcut -c name  > ./list/SLC_names_final.txt
#cat $B | tr -s ' ' ',' | csvcut -c name  > ./list/SLC_names_final.txt

## extract fasta of gene lists from Dmel
~/Applications/custom/unigene_fa_sub.sh ./reference_proteome/proteome_SLC_mark.fa  ./list/SLC_names_final.txt > ./list/SLC_genes.fa


## divide each family into independent fasta
mkdir ./family_fasta 
rm ./family_fasta/*
IFS=$'\n'; 
for next in $(cat ~/Documents/SLC_id/general_reference/SLC_info/SLC_families.txt)
do 
  grep -A 1 ${next}"_" ./list/SLC_genes.fa | sed '/--/d' > ./family_fasta/${next}.fa
done

## align all sequences
mkdir ./muscle_alignments
for i in ./family_fasta/*
do
~/Applications/muscle3.8.31_i86linux64 -in $i -out $i.aln
~/Applications/trimal-trimAl/source/trimal -in $i.aln -out $i.trimmed
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


cd $3
#mkdir dictionary_create_FIND
cp $2 $(pwd $3)'/SLC_dict.csv' ## needs to change name
#cp $B $(pwd $C)'/dictionary_create_FIND/SLC_dict.csv' ## needs to change name
#cp $B $(pwd $C'/SLC_dict.csv')
#cp ~/Documents/SLC_id/general_reference/Hs_SLC_dict.csv ~/Documents/SLC_id/Human_HMM_SLC'/SLC_dict.csv'

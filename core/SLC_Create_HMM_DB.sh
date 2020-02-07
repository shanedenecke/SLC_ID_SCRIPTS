#!/usr/bin/env bash

### ###################### Create HMM profiles for each SLC family in species ###################

## argument 1= Predicted unigene protein dataset
## argument 2= SLC dictionary with codes and names
## argument 3= Name of database


### This script will create a HMM database for each SLC family in a given species. The inputs are a uniprotein set, 
### an SLC dictionary (gene codes and SLC names) and a new directory that will contain the database


## use these for debugging

#A=$H/GENERAL_REFERENCE/model_proteomes/HomSap_unigene.faa
#B=$H/GENERAL_REFERENCE/model_SLC_info/HomSap_SLC_dict.csv
#C=$H/HomSap_Database
################################################################################
#

## reset bash
#source ~/.bashrc

## Create and set working directory
mkdir $3
cd $3

#mkdir $C
#cd $C


## rename fasta files with tag for SLC family
mkdir reference_proteome
$H/SLC_ID_SCRIPTS/general_scripts/fasta_rename.py $1 $2 > ./reference_proteome/proteome_SLC_mark.fa
#$H/SLC_ID_SCRIPTS/general_scripts/fasta_rename.py $A $B > ./reference_proteome/proteome_SLC_mark.fa
#/data2/shane/Applications/custom/fasta_rename.py $A $B > ./reference_proteome/proteome_SLC_mark.fa

#/data2/shane/Applications/custom/fasta_rename.py /data2/shane/Documents/omics_data/Homo_sapiens/unigene/Hs_unigene.fa /data2/shane/Dropbox/wp7_prodrug/SLC_id/general_reference/Hs_SLC_dict.csv > ./reference_proteome/proteome_SLC_mark.fa

##Make blast dbs for target reference proteome
makeblastdb -in ./reference_proteome/proteome_SLC_mark.fa -parse_seqids -dbtype prot

## extract list of SLC names to be used in future steps
mkdir ./list
#cat $2 | tr -s ' ' ',' | csvcut -c name  > ./list/SLC_names_final.txt
#cat $B | tr -s ' ' ',' | csvcut -c name  > ./list/SLC_names_final.txt

num=$(head -1 $2 | tr ',' '\n' | cat -n | grep "name")
#num=$(head -1 $B | tr ',' '\n' | cat -n | grep "name")
num2=$(echo $num | cut -f 1 | sed -E 's/\s+//g') 
cut -d',' -f $num2 $2 | tail -n+2 > ./list/SLC_names_final.txt
#cut -d',' -f $num2 $B | tail -n+2 > ./list/SLC_names_final.txt

## extract fasta of gene lists from Dmel
/data2/shane/Applications/custom/unigene_fa_sub.sh ./reference_proteome/proteome_SLC_mark.fa  ./list/SLC_names_final.txt > ./list/SLC_genes.fa


## divide each family into independent fasta
mkdir ./family_fasta 
rm -rf ./family_fasta/*
IFS=$'\n'; 
for next in $(cat $H/GENERAL_REFERENCE/family_species_lists_phylo/SLC_families.txt)
do 
  grep -A 1 ${next} ./list/SLC_genes.fa | sed '/--/d' > ./family_fasta/${next}.fa
done

## align all sequences
mkdir ./muscle_alignments
for i in ./family_fasta/*
do
seq=$(cat $i | wc -l)

if [ $seq -gt 2 ]
  then
    mafft --thread 24 $i > $i.aln
    #/home/pioannidis/Programs/trimAl/source/trimal -in $i.aln -out $i.trimmed
    $H'/SLC_ID_SCRIPTS/general_scripts/trimAl/source/trimal' -in $i.aln -out $i.trimmed
  else
    cat $i > $i.trimmed
  fi
done

mv ./family_fasta/*.trimmed ./muscle_alignments
rm -rf ./family_fasta/*.aln

## build HMM profiles for each family
mkdir ./hmm_profiles
for i in ./muscle_alignments/*
do
#hmmbuild --cpu $threads $i.hmm $i
hmmbuild $i.hmm $i
done
mv ./muscle_alignments/*.hmm ./hmm_profiles


cd $3
#mkdir dictionary_create_FIND
cp $2 $(pwd $3)'/SLC_source_dict.csv' ## needs to change name
#cp $B $(pwd $C)'/dictionary_create_FIND/SLC_dict.csv' ## needs to change name
#cp $B $(pwd $C'/SLC_dict.csv')
cd $H
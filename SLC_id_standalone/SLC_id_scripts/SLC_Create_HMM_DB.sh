#!/usr/bin/env bash

### ###################### Create HMM profiles for each SLC family in species ###################

### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteome: Predicted unigene protein dataset of organism
  -dict: SLC dictionary with codes and names
  -out: Name of database to be created. This will output a folder with this name
  -threads: threads "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteome) PROTEOME="$2"; shift 2;;
    -dict) DICT="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -out) OUT="$2"; shift 2;;
  esac
done

#PROTEOME=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_reference/DroMel_unigene.faa
#DICT=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_reference/HomSap_SLC_dict.csv
#THREADS=14
#OUT=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/DroMel_Database
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname "$SCRIPT_DIR")"

###
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname "$SCRIPT_DIR")"


## Create and set working directory
mkdir $OUT
cd $OUT

## rename fasta files with tag for SLC family
mkdir reference_proteome
$SCRIPT_DIR/fasta_rename.py $PROTEOME $DICT > ./reference_proteome/proteome_SLC_mark.fa

##Make blast dbs for target reference proteome
makeblastdb -in ./reference_proteome/proteome_SLC_mark.fa -parse_seqids -dbtype prot

## extract list of SLC names to be used in future steps
mkdir ./list
cut -d',' -f 2 $DICT | tail -n+2 > ./list/SLC_names_final.txt

## extract fasta of gene lists from Dmel
$SCRIPT_DIR/unigene_fa_sub.sh ./reference_proteome/proteome_SLC_mark.fa  ./list/SLC_names_final.txt > ./list/SLC_genes.fa


## divide each family into independent fasta
mkdir ./family_fasta 
rm -rf ./family_fasta/*
IFS=$'\n'; 
for next in $(cat $SOURCE_DIR/SLC_id_reference/SLC_Families.txt)
do 
  grep -A 1 ${next} ./list/SLC_genes.fa | sed '/--/d' > ./family_fasta/${next}.fa
done

## align all sequences
mkdir ./alignments
for i in ./family_fasta/*
do
seq=$(cat $i | wc -l)

if [ $seq -gt 2 ]
  then
    mafft-linsi --quiet --thread $THREADS $i > $i.aln
    ~/Applications/trimal/source/trimal -automated1 -in $i.aln -out $i.trimmed
  else
    cat $i > $i.trimmed
  fi
done

mv ./family_fasta/*.trimmed ./alignments
rm -rf ./family_fasta/*.aln

## build HMM profiles for each family
mkdir ./hmm_profiles
for i in ./alignments/*
do
hmmbuild --cpu $THREADS $i.hmm $i
done
mv ./alignments/*.hmm ./hmm_profiles

### copy starting dictionary over
cp $DICT $OUT/SLC_source_dict.csv

cd $OUT


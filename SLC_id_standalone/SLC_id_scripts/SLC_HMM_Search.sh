#!/usr/bin/env bash

### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -database: starting database. Pathw to folder of model species you will use to search target species
  -target: path to target species proteome
  -out: Name of output folder. Will create a new SLC database in this folder for the target species. 
  -threads: threads"
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -database) DATABASE="$2"; shift 2;;
    -target) TARGET="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -out) OUT="$2"; shift 2;;
  esac
done

### Set scripts directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname "$SCRIPT_DIR")"

#For Debugging
#TARGET=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_reference/DroMel_unigene.faa
#DATABASE=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/HomSap_Database
#THREADS=14
#OUT=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/Hs_to_DroMel_Search
#SCRIPT_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_scripts/
#SOURCE_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/

#echo 'Target is '$TARGET
#echo 'Database is '$DATABASE
#echo 'Out is '$OUT
#echo 'Script dir '$SCRIPT_DIR
#echo 'Source dir '$SOUCRE_DIR



## Create new output directory
mkdir $OUT
cd $OUT

echo 'Performing HMM search'
## search HMM profiles against target proteome using hmm search
rm -rf ./hmm_outputs
mkdir ./hmm_outputs
for i in $DATABASE/hmm_profiles/*; do
  base=$(echo $(basename $i))
  hmmsearch --notextw -E 20 $i $TARGET > ./hmm_outputs/$base.hmmoutput
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
  cut -f 9 $i | sed 's/\s+//g'| $SCRIPT_DIR/unigene_fa_sub.sh $TARGET - > ./SLC_fa/$base'.fa'
done
find ./SLC_fa/* -size 0 -delete

##perform blast with HMM hits as queries and the source genomes SLC_mark.fa proteome as a target 
echo 'Blast away'
rm -rf ./recip_blast
mkdir ./recip_blast
for i in ./SLC_fa/*.fa; do
  base=$(echo $(basename $i))
  blastp -query $i -db $DATABASE/reference_proteome/proteome_SLC_mark.fa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 6 -max_hsps 1 -num_threads $THREADS > ./recip_blast/$base'_blast.tsv'
done
find ./recip_blast/* -size 0 -delete


## Run R script to output table 
mkdir ./prelim_summary
Rscript $SCRIPT_DIR/SLC_Family_Sort.R $DATABASE'/SLC_source_dict.csv' >  ./prelim_summary/Family_sort_preliminary.csv

## Filter based on lengths of human SLC gene
mkdir length_analysis

### make fasta from reciprocal blast results
cut -d ',' -f 1 ./prelim_summary/Family_sort_preliminary.csv | $SCRIPT_DIR/unigene_fa_sub.sh $TARGET - > ./length_analysis/preliminary_SLC.fa
cut -d ',' -f 2 ./prelim_summary/Family_sort_preliminary.csv | sed '1d' > ./length_analysis/length_families.txt 
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ./length_analysis/preliminary_SLC.fa > ./length_analysis/names_lengths.txt
grep ">" ./length_analysis/names_lengths.txt | perl -pe  's/^>(.+$)/$1/;'| cut -d ' ' -f 1  > ./length_analysis/all_proteins.txt
grep -E "^[0-9]" ./length_analysis/names_lengths.txt > ./length_analysis/all_lengths.txt
paste -d',' ./length_analysis/all_proteins.txt ./length_analysis/all_lengths.txt ./length_analysis/length_families.txt > ./length_analysis/gene_lengths.txt
Rscript $SCRIPT_DIR/SLC_length_filter.R $SOURCE_DIR > ./length_analysis/total_slc_table.csv

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
Rscript $SCRIPT_DIR/SLC_dictionary_format.R > ./final_output/SLC_final_output.csv

cd -

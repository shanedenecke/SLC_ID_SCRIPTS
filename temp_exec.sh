#!/usr/bin/env bash

cd /data2/shane/Documents/SLC_id
mkdir genome_score

######################## 0) Download and Clean Sequences 
#./SLC_id_scripts/SLC_Proteome_prepare.sh

######################## 1) Build Human Database
./SLC_id_scripts/SLC_Create_HMM_DB.sh /data2/shane/Documents/SLC_id/general_reference/model_proteomes/HomSap_unigene.faa /data2/shane/Documents/SLC_id/general_reference/SLC_info/HomSap_SLC_dict.csv /data2/shane/Documents/SLC_id/HomSap_Database



######################## 2) Bulild good quality Drosophila database
mkdir Dm_Database_Generate

## search from humans 
./SLC_id_scripts/SLC_HMM_Search.sh /data2/shane/Documents/SLC_id/HomSap_Database /data2/shane/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa /data2/shane/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Search

## Xref with flybase
Rscript ./SLC_id_scripts/SLC_Flybase_human_SLCxref.R > ./Dm_Database_Generate/SLC_source_dict_flybaseXref.csv

## Make final Drosophila database
./SLC_id_scripts/SLC_Create_HMM_DB.sh /data2/shane/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa /data2/shane/Documents/SLC_id/Dm_Database_Generate/SLC_source_dict_flybaseXref.csv /data2/shane/Documents/SLC_id/DroMel_Database


########################  3) Search species with Human database
mkdir Human_search
for i in /data2/shane/Documents/SLC_id/proteomes/*.faa; do
a=$(basename $i) 
echo 'HUMAN SEARCH '$a
./SLC_id_scripts/SLC_HMM_Search.sh /data2/shane/Documents/SLC_id/HomSap_Database $i /data2/shane/Documents/SLC_id/Human_search/'HUMAN_'$a
done

########################  4) Search other species with Drosohpila database
mkdir Drosophila_search
for i in /data2/shane/Documents/SLC_id/proteomes/*.fa*; do
b=$(echo $(basename $i) | cut -d '_' -f 1) 
echo 'DROSOPHILA SEARCH '$b
./SLC_id_scripts/SLC_HMM_Search.sh /data2/shane/Documents/SLC_id/DroMel_Database $i /data2/shane/Documents/SLC_id/Drosophila_search/'DROSOPHILA_'$b
done


########################  5) collate searches for each species by Xref with Human and Drosophila searches

Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_crossref_human_dros_searches.R

####################### 6) summarize counts of SLC tables

#Rscript ./SLC_id_scripts/SLC_family_count_combine.R
Rscript ./SLC_id_scripts/SLC_id_PPP.R
#!/usr/bin/env bash

cd ~/Documents/SLC_id/


######################## 0) Download and Clean Sequences 
python3 ./SLC_id_scripts/SLC_download_clean_genomes.py
cp ./general_reference/model_proteomes/HarArm_unigene.faa ./proteomes/
find ./proteomes -type f -empty -delete

######################## 1) Build Human Database
./SLC_id_scripts/SLC_Create_HMM_DB.sh ~/Documents/SLC_id/general_reference/model_proteomes/HomSap_unigene.faa ~/Documents/SLC_id/general_reference/SLC_info/HomSap_SLC_dict_new.csv ~/Documents/SLC_id/HomSap_Database



######################## 2) Bulild good quality Drosophila database
mkdir Dm_Database_Generate

## search from humans 
./SLC_id_scripts/SLC_HMM_Search.sh ~/Documents/SLC_id/HomSap_Database ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Search

## Make Drosophila database from Human search
./SLC_id_scripts/SLC_Create_HMM_DB.sh ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Search/final_output/SLC_final_output.csv ~/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Database

# Iteratative search Drosophila
./SLC_id_scripts/SLC_HMM_Search.sh ~/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Database ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Dm_Database_Generate/DroMel_iterative_search

## Xref with flybase SLC calls
Rscript ./SLC_id_scripts/SLC_Flybase_human_SLCxref.R > ./Dm_Database_Generate/SLC_source_dict_flybaseXref.csv


## Make final Drosophila database
./SLC_id_scripts/SLC_Create_HMM_DB.sh ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Dm_Database_Generate/SLC_source_dict_flybaseXref.csv ~/Documents/SLC_id/DroMel_Database


########################  3) Search species with Human database
mkdir Human_search
for i in ~/Documents/SLC_id/proteomes/*.faa; do
a=$(basename $i) 
./SLC_id_scripts/SLC_HMM_Search.sh ~/Documents/SLC_id/HomSap_Database $i ~/Documents/SLC_id/Human_search/'HUMAN_'$a
done

########################  4) Search other species with Drosohpila database
mkdir Drosophila_search
for i in ~/Documents/SLC_id/proteomes/*.fa*; do
b=$(echo $(basename $i) | cut -d '_' -f 1) 
./SLC_id_scripts/SLC_HMM_Search.sh ~/Documents/SLC_id/DroMel_Database $i ~/Documents/SLC_id/Drosophila_search/'DROSOPHILA_'$b
done


########################  5) collate searches for each species by Xref with Human and Drosophila searches

Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_crossref_human_dros_searches.R


######################## 6) Iteratively search each species

## create databases for each species
mkdir iterative_database
for i in  ~/Documents/SLC_id/proteomes/*.faa; do
  c=$(echo $(basename $i) | cut -d '_' -f 1) 
  d=$(ls ~/Documents/SLC_id/Human_Drosophila_crossref/$c*)
  ./SLC_id_scripts/SLC_Create_HMM_DB.sh $i $d ~/Documents/SLC_id/iterative_database/'iterative_database_'$c
done


## use iterative databases to recursivley search genomes
mkdir iterative_search
mkdir final_SLC_dicts
for i in  ~/Documents/SLC_id/proteomes/*.faa; do
e=$(echo $(basename $i) | cut -d '_' -f 1) 
f=$(ls -d ~/Documents/SLC_id/iterative_database/iterative_database*$e)
echo 'now performing '$e' analysis'
./SLC_id_scripts/SLC_HMM_Search.sh $f $i ~/Documents/SLC_id/iterative_search/'iterative_search_'$e
cp ~/Documents/SLC_id/iterative_search/'iterative_search_'$e/final_output/SLC_final_output.csv ~/Documents/SLC_id/final_SLC_dicts/$e'Final_SLC_table.csv'
done


###################### 7) summarize counts of SLC tables
#python3 ./SLC_id_scripts/SLC_summary_count.py

###################### 8) Extract sequences from each relevant species. Perform alignment and phylogeny
./SLC_id_scripts/SLC_id_Align_and_tree.sh
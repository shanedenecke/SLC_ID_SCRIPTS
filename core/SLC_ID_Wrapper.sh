#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/SLC_id'
PHYLO=$H/GENERAL_REFERENCE/input_arguments/Phylo_list.txt
SLC_FAM=$H/GENERAL_REFERENCE/input_arguments/SLC_Families.txt
SPEC=$H/GENERAL_REFERENCE/input_arguments/target_species.tsv
THREADS=12

cd $H
mkdir genome_score
mkdir intermediate

######################## 0) Download and Clean Sequences 
source ./SLC_ID_SCRIPTS/core/SLC_Proteome_prepare.sh

######################## 1) Build Human Database
source ./SLC_ID_SCRIPTS/core/SLC_Create_HMM_DB.sh $H/GENERAL_REFERENCE/model_proteomes/HomSap_unigene.faa $H/GENERAL_REFERENCE/model_SLC_info/HomSap_SLC_dict.csv $H/HomSap_Database



######################## 2) Bulild good quality Drosophila database
mkdir Dm_Database_Generate

## search from humans 
source ./SLC_ID_SCRIPTS/HMM_Search/SLC_HMM_Search.sh $H/HomSap_Database $H/GENERAL_REFERENCE/model_proteomes/DroMel_unigene.faa $H/Dm_Database_Generate/Hs_to_DroMel_Search

## Xref with flybase
Rscript ./SLC_ID_SCRIPTS/core/SLC_Flybase_human_SLCxref.R $H > ./Dm_Database_Generate/SLC_source_dict_flybaseXref.csv

## Make final Drosophila database
source ./SLC_ID_SCRIPTS/core/SLC_Create_HMM_DB.sh $H/GENERAL_REFERENCE/model_proteomes/DroMel_unigene.faa $H/Dm_Database_Generate/SLC_source_dict_flybaseXref.csv $H/DroMel_Database
mv Dm_Database_Generate intermediate

########################  3) Search species with Human database
mkdir Human_search
for i in $H/proteomes/*.faa; do
  a=$(basename $i) 
  echo 'HUMAN SEARCH '$a
  source ./SLC_ID_SCRIPTS/HMM_Search/SLC_HMM_Search.sh $H/HomSap_Database $i $H/Human_search/'HUMAN_'$a
done

########################  4) Search other species with Drosohpila database
mkdir Drosophila_search
for i in $H/proteomes/*.faa; do
  b=$(echo $(basename $i) | cut -d '_' -f 1) 
  echo 'DROSOPHILA SEARCH '$b
  source ./SLC_ID_SCRIPTS/HMM_Search/SLC_HMM_Search.sh $H/DroMel_Database $i $H/Drosophila_search/'DROSOPHILA_'$b
done


########################  5) collate searches for each species by Xref with Human and Drosophila searches

Rscript ./SLC_ID_SCRIPTS/core/SLC_crossref_human_dros_searches.R $H



####################### 6) Post Process SLC tables
Rscript ./SLC_ID_SCRIPTS/Post_ID_summary/SLC_id_PPP.R $H
Rscript ./SLC_ID_SCRIPTS/Post_ID_summary/SLC_Figures.R $H ### Maybe needs edits. Can't tell with low number of families

### clean up folder
mv DroMel_Database ./intermediate
mv HomSap_Database ./intermediate
mv Drosophila_search ./intermediate
mv Human_search ./intermediate
mv TMHMM_filter ./intermediate/
mv preliminary_SLC_dicts intermediate
rm -r TMHMM_*

###################### 7) Extract sequences from each relevant species. Perform alignment and phylogeny
source ./SLC_ID_SCRIPTS/Align_Tree/SLC_id_Align_and_tree.sh



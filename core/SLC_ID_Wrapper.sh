#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/SLC_id'
cd $H
mkdir genome_score

######################## 0) Download and Clean Sequences 
source ./SLC_ID_SCRIPTS/SLC_Proteome_prepare.sh

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
rm -rf Dm_Database_Generate

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

### clean up folder
mkdir intermediate
mv DroMel_Database ./intermediate
mv HomSap_Database ./intermediate
mv Drosophila_search ./intermediate
mv Human_search ./intermediate

####################### 6) summarize counts of SLC tables

#Rscript ./SLC_ID_SCRIPTS/SLC_family_count_combine.R
### Ultrametric tree generate
Rscript ./SLC_ID_SCRIPTS/Post_ID_summary/SLC_id_PPP.R $H
Rscript ./SLC_ID_SCRIPTS/Post_ID_summary/SLC_Figures.R $H

###################### 7) Extract sequences from each relevant species. Perform alignment and phylogeny
source ./SLC_ID_SCRIPTS/SLC_id_Align_and_tree.sh


##################### 8) Prepare Ultrametric Tree
#./SLC_id/scripts/CAFE/SLC_ultrametric_tree_prepare.sh
./SLC_ID_SCRIPTS/CAFE/Ultrametric_tree_generate.sh


## RUN CAFE
./SLC_ID_SCRIPTS/CAFE/CAFE_run_full.sh
Rscript ./SLC_ID_SCRIPTS/CAFE/CAFE_figures.R

##################### 9) Prepare shiny material
#Rscript ./SLC_id/scripts/shiny/SLC_shiny_prepare.R

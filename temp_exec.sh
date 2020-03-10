#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/SLC_id'
PHYLO=$H/GENERAL_REFERENCE/input_arguments/Phylo_list.txt
SLC_FAM=$H/GENERAL_REFERENCE/input_arguments/SLC_Families.txt
SPEC=$H/GENERAL_REFERENCE/input_arguments/target_species.tsv
THREADS=32


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

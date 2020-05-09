#!/bin/bash
H=/mnt/disk/shane/Transporter_ID/SLC_id
PHYLO=$H/GENERAL_REFERENCE/CAFE/Phylo_list.txt
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
THREADS=14


cd $H

########################  3) Search species with Human database
rm ./proteomes/HomSap_unigene.faa
rm ./proteomes/DroMel_unigene.faa
mkdir Human_search
for i in $H/proteomes/*.faa; do
  a=$(basename $i) 
  echo 'HUMAN SEARCH '$a
  source ./SLC_ID_SCRIPTS/SLC_id/SLC_HMM_Search.sh $H/HomSap_Database $i $H/Human_search/'HUMAN_'$a
done

########################  4) Search other species with Drosohpila database
mkdir Drosophila_search
for i in $H/proteomes/*.faa; do
  b=$(echo $(basename $i) | cut -d '_' -f 1) 
  echo 'DROSOPHILA SEARCH '$b
  source ./SLC_ID_SCRIPTS/SLC_id/SLC_HMM_Search.sh $H/DroMel_Database $i $H/Drosophila_search/'DROSOPHILA_'$b
done

#!/usr/bin/env bash

cd /data2/shane/Documents/SLC_id


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

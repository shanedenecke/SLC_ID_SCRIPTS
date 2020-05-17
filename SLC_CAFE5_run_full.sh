#!/usr/bin/env bash

## create CAFE script files
rm -r ./CAFE/scripts
rm -r ./CAFE/logfiles
rm -r ./CAFE/outputs


#mkdir ./CAFE/scripts
#mkdir ./CAFE/logfiles
mkdir ./CAFE/outputs

for i in ./CAFE/CAFE_tables/*.tsv
do  
  b=$(echo $(basename $i) | sed 's/_SLC_CAFE_table.tsv//g')
  cafexp -i $i -o ./CAFE/outputs/$b -t ./CAFE/clean_raxml_trees/$b'_tree_ultrametric.nwk' 
done  


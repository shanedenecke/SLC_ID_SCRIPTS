#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline
PHYLO=$H/GENERAL_REFERENCE/CAFE/Phylo_list.txt
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
THREADS=10


cd $H
for i in ./phylogeny/phylip/*.phy
do
  b=$(echo $(basename $i) | cut -d '_' -f 1,2) 
  raxfile=$i
  raxdir=$H/phylogeny/raxml_trees
  
  ### hash next two lines if you don't want to actually make the trees. 
  ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $b'.nwk' -w $raxdir 
  cp $H/phylogeny/raxml_trees/RAxML_bipartitions."$b".nwk ./phylogeny/clean_raxml_trees/
done

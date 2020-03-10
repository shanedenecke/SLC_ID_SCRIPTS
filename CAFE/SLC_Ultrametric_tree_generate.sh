#!/usr/bin/env bash

mkdir CAFE
mkdir ./CAFE/Ultrametric_tree
mkdir ./CAFE/species_lists
mkdir ./CAFE/clean_raxml_trees


for i in ./GENERAL_REFERENCE/CAFE/*preliminary_species.txt; do
  b=$(echo $(basename $i) | sed 's/.+CAFE.//g' | sed 's/_preliminary_species.txt//g')
  cp $i ./CAFE/species_lists/$b'_preliminary_species.txt'
  grep -f ./Final_raw_outputs/Poor_quality_species.txt ./CAFE/species_lists/$b'_preliminary_species.txt' > ./CAFE/species_lists/$b'_final_species.txt'
done


## Create ultrametric Trees
for i in ./CAFE/species_lists/*final*
do
  b=$(echo $(basename $i) | sed 's/.+CAFE.//g' | sed 's/_final_species.txt//g')
  
  if [ $b = "Arthropod" ]; then
    /data2/shane/Applications/custom/OrthoDB/one_to_one_ID_exec.py -node "Metazoa" -taxid $i -output seq
  else
    /data2/shane/Applications/custom/OrthoDB/one_to_one_ID_exec.py -node "Arthropod" -taxid $i -output seq
  fi
  
  ### rename identified foler of sequences 
  fulltemp=$(readlink -f './CAFE/Ultrametric_tree/'$b'_og_sequences')
  mv og_sequences $fulltemp
  
  
  ### perform alignments for all one to ones 
  for x in  $fulltemp/*
  do
    mafft --quiet --thread $THREADS $x > $x'.aln'
    /data2/shane/Applications/trimAl/source/trimal -in $x'.aln' -out $x'.aln.trimm'
    /data2/shane/Applications/custom/fasta_2_phylip.sh $x'.aln.trimm' | sed '1d' > $x'.aln.trimm.phy'
  done
  
  #### merge all phylip files
  Rscript ./ABC_ID_SCRIPTS/ABC_Phylip_merge.R $fulltemp
  
  #### make trees 
  if [ $b = "Arthropod" ]; then
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp -o 6239_0
  elif [ $b = 'Lepidopteran' ]; then
	echo 'lep'
	#/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp -o 7029_0
  elif [ $b = "Hemipteran" ]; then
	sed -i 's/J/A/g' $fulltemp/Full_species.phy
	sed -i 's/\./A/g' $fulltemp/Full_species.phy
	echo "Hemi"
	#/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp -o 7227_0
  elif [ $b = 'Arachnid' ]; then
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp
 fi
 cp $fulltemp'/RAxML_bipartitions.'$b'.tre' ./CAFE/clean_raxml_trees/
done

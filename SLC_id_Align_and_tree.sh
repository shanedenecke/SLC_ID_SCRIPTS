################### ALIGN AND TREE #####################


## copy Dromel and Homsap tables to final dictionary file
#Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_Rename_Dm_SLCs.R > ./final_SLC_dicts/DroMelFinal_SLC_table.csv

#cp ./DroMel_Database/SLC_source_dict.csv ./final_SLC_dicts/DroMelFinal_SLC_table.csv
#cp ./HomSap_Database/SLC_source_dict.csv ./final_SLC_dicts/HomSapFinal_SLC_table.csv

mkdir ./phylogeny
mkdir ./phylogeny/renamed_dicts
mkdir ./phylogeny/SLC_fa
rm ./phylogeny/renamed_dicts/*
rm ./phylogeny/SLC_fa/*
specieslist=SpoFru,DroMel,NezVir,HelArm,TriCas
for i in ${specieslist//,/ }
do
  cp ./Final_outputs/Final_SLC_dicts/$i'_final_SLC_table.csv' ./phylogeny/renamed_dicts/
  num=$(head -1 ./Final_outputs/Final_SLC_dicts/$i'_final_SLC_table.csv' | tr ',' '\n' | cat -n | sed -E 's/\s+/,/g' | sed -E 's/^,//g' | grep "code")
  num2=$(echo $num | cut -d',' -f 1) 
  cut -d',' -f $num2 ./Final_outputs/Final_SLC_dicts/$i'_final_SLC_table.csv' | sed 's/"//g' > cod.txt
  if [ $i == 'DroMel' ] || [ $i == 'HomSap' ]
  then
    ~/Applications/Custom_Applications/unigene_fa_sub.sh ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_reference/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa'
  else
    ~/Applications/Custom_Applications/unigene_fa_sub.sh ./proteomes/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa' #||  ~/Applications/Custom_Applications/unigene_fa_sub.sh ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_reference/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa'
  fi
   
   ~/Applications/Custom_Applications/fasta_rename.py ./phylogeny/SLC_fa/$i'_SLC.faa' ./phylogeny/renamed_dicts/$i'_final_SLC_table.csv' >> ./phylogeny/SLC_fa/combined_renamed.faa
rm cod.txt
done

  
  
  
mkdir ./phylogeny/SLC_byfam  
mkdir ./phylogeny/alignments
mkdir ./phylogeny/trimms
mkdir ./phylogeny/phylip
cat ./GENERAL_REFERENCE/SLC_Families.txt | while read i
do
  grep -E -A 1 $i ./phylogeny/SLC_fa/combined_renamed.faa | sed '/--/d' > './phylogeny/SLC_byfam/'$i'phylo_subset.faa'
  mafft-linsi --thread $THREADS './phylogeny/SLC_byfam/'$i'phylo_subset.faa' > './phylogeny/alignments/'$i'phylo_subset.faa.aln'
  ~/Applications/trimal/source/trimal -automated1 -phylip_paml -in './phylogeny/alignments/'$i'phylo_subset.faa.aln' -out './phylogeny/trimms/'$i'phylo_subset.faa.aln.trimm'
  Rscript ~/Applications/Custom_Applications/Phylip_duplicate.R './phylogeny/trimms/'$i'phylo_subset.faa.aln.trimm' > './phylogeny/phylip/'$i'phylo_subset.faa.aln.trimm.phy'
  sed -i 's/\;/_/g' './phylogeny/phylip/'$i'phylo_subset.faa.aln.trimm.phy'
done

mkdir phylogeny/raxml_trees
mkdir phylogeny/clean_raxml_trees
rm phylogeny/raxml_trees/*
rm phylogeny/clean_raxml_trees/*
for i in ./phylogeny/phylip/*.phy
do
  b=$(echo $(basename $i) | cut -d '_' -f 1,2) 
  raxfile=$i
  raxdir=$H/phylogeny/raxml_trees
  
  ### hash next two lines if you don't want to actually make the trees. 
  ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $b'.nwk' -w $raxdir 
  cp $H/phylogeny/raxml_trees/RAxML_bipartitions."$b".nwk ./phylogeny/clean_raxml_trees/
done


mv phylogeny ./intermediate/

#Rscript 


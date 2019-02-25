################### ALIGN AND TREE #####################


## copy Dromel and Homsap tables to final dictionary file
Rscript ~/Documents/SLC_id/SLC_id_scripts/SLC_Rename_Dm_SLCs.R > ./final_SLC_dicts/DroMelFinal_SLC_table.csv
cp ./HomSap_Database/SLC_source_dict.csv ./final_SLC_dicts/HomSapFinal_SLC_table.csv


#### Make new directory for model insect dictionaries
mkdir renamed_SLC_dicts
mkdir marked_SLC_for_aligment
rm ./marked_SLC_for_aligment/SLC_all.fa
touch ./marked_SLC_for_aligment/SLC_all.fa

grep -A 1 "SLC_" ./DroMel_Database/reference_proteome/proteome_SLC_mark.fa | perl -pe 's/>/>DroMel_/' >> ./marked_SLC_for_aligment/SLC_all.fa
grep -A 1 "SLC_" ./HomSap_Database/reference_proteome/proteome_SLC_mark.fa | perl -pe 's/>/>HomSap_/'  >> ./marked_SLC_for_aligment/SLC_all.fa

## append species name to each SLC name and move to new directory
cat ./general_reference/SLC_info/important_species_codes.txt | while read i
do
  cut -d ',' -f 1 ./final_SLC_dicts/$i* > temp.txt ### Put species prefix on each SLC 
  cat ./final_SLC_dicts/$i* | cut -d ',' -f 2 | sed -e "s/^SLC/$i\_SLC/" | paste --delimiter ',' temp.txt - > ./renamed_SLC_dicts/$i'_Annotated_SLC_dict.csv' ### Put species prefix on each SLC 
  ~/Applications/custom/fasta_rename.py ./proteomes/$i* ./renamed_SLC_dicts/$i'_Annotated_SLC_dict.csv' | grep -A 1 "SLC_" >> ./marked_SLC_for_aligment/SLC_all.fa
  rm temp.txt
done
find ./renamed_SLC_dicts/* -size 0 -delete

## Divide each SLC family by SLC_###_ marker align and tree
mkdir SLC_phylogeny
cat ./general_reference/SLC_info/SLC_families.txt | while read i
do
  a=./SLC_phylogeny/$i'.fa'
  grep -A 1 $i ./marked_SLC_for_aligment/SLC_all.fa | sed '/--/d' > $a ## isolate sequences in fasta format
  ~/Applications/muscle3.8.31_i86linux64 -in $a -out $a'.aln'
  /home/pioannidis/Programs/trimAl/source/trimal -in $a'.aln' -out $a'.aln.trimm'
  #~/Applications/trimal-trimAl/source/trimal -in $a'.aln' -out $a'.aln.trimm' ## LOCAL
  ~/Applications/custom/fasta_2_phylip.sh $a'.aln.trimm' > $a'.aln.trimm.phy'
  
  raxfile=$(readlink -f  $a'.aln.trimm.phy')
  raxdir=$(readlink -f ./SLC_phylogeny)
  rm ./SLC_phylogeny/RAxML*
  ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 2 -T 10 -m PROTGAMMAAUTO -s $raxfile -n 'ABC_phylo' -w $raxdir 
  #~/Applications/standard-RAxML-master/raxmlHPC-AVX -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s $raxfile -n $i'.tre' -w $raxdir ## LOCAL

done






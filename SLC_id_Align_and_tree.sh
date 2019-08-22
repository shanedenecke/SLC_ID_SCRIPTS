################### ALIGN AND TREE #####################
cd /data2/shane/Documents/SLC_id

## copy Dromel and Homsap tables to final dictionary file
#Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_Rename_Dm_SLCs.R > ./final_SLC_dicts/DroMelFinal_SLC_table.csv

#cp ./DroMel_Database/SLC_source_dict.csv ./final_SLC_dicts/DroMelFinal_SLC_table.csv
#cp ./HomSap_Database/SLC_source_dict.csv ./final_SLC_dicts/HomSapFinal_SLC_table.csv

mkdir phylogeny
mkdir ./phylogeny/renamed_dicts
mkdir ./phylogeny/SLC_fa
rm ./phylogeny/renamed_dicts/*
rm ./phylogeny/SLC_fa/*
cat ./general_reference/SLC_info/Phylo_list.txt | while read i
do
  cp ./final_SLC_dicts/$i'Final_SLC_table.csv' ./phylogeny/renamed_dicts/
  csvcut -c name ./phylogeny/renamed_dicts/$i'Final_SLC_table.csv' | sed 's/"//g' > nam.txt
  csvcut -c code ./phylogeny/renamed_dicts/$i'Final_SLC_table.csv' | sed 's/"//g' > cod.txt
  paste  -d '_' nam.txt cod.txt | sed -e "s/^/"$i"_/g" | sed -e '1s/.*/name/g' > newcol.txt
  paste -d ',' newcol.txt cod.txt > ./phylogeny/renamed_dicts/$i'Final_SLC_table.csv'
  
  if [ $i == 'DroMel' ] || [ $i == 'HomSap' ]
  then
    /data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa'
  else
    /data2/shane/Applications/custom/unigene_fa_sub.sh ./proteome_clean/clean_fasta/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa' || /data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa'
  fi
  /data2/shane/Applications/custom/fasta_rename.py ./phylogeny/SLC_fa/$i'_SLC.faa' ./phylogeny/renamed_dicts/$i'Final_SLC_table.csv' >> ./phylogeny/SLC_fa/combined_renamed.faa
done
  
  
  
mkdir ./phylogeny/SLC_byfam  
mkdir ./phylogeny/alignments
mkdir ./phylogeny/trimms
mkdir ./phylogeny/phylip
cat ./general_reference/SLC_info/SLC_families.txt | while read i
do
#echo $i
#done
  grep -A 1 $i ./phylogeny/SLC_fa/combined_renamed.faa | sed '/--/d' > './phylogeny/SLC_byfam/'$i'phylo_subset.faa'
  mafft --threadtb 24 './phylogeny/SLC_byfam/'$i'phylo_subset.faa' > './phylogeny/alignments/'$i'phylo_subset.faa.aln'
  /home/pioannidis/Programs/trimAl/source/trimal -in './phylogeny/alignments/'$i'phylo_subset.faa.aln' -out './phylogeny/trimms/'$i'phylo_subset.faa.aln.trimm'
  /data2/shane/Applications/custom/fasta_2_phylip.sh './phylogeny/trimms/'$i'phylo_subset.faa.aln.trimm' > './phylogeny/phylip/'$i'phylo_subset.faa.aln.trimm.phy'
done

  
mkdir SLC_phylogeny
for i in /data2/shane/Documents/SLC_id/phylogeny/phylip/*.phy
do
  b=$(echo $(basename $i) | cut -d '_' -f 1,2) 
  raxfile=$i
  raxdir=/data2/shane/Documents/SLC_id/SLC_phylogeny/
  #rm ./SLC_phylogeny/RAxML*
  /data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T 36 -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
  #/data2/shane/Applications/standard-RAxML-master/raxmlHPC-AVX -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s $raxfile -n $i'.tre' -w $raxdir ## LOCAL
done
  







#### Make new directory for model insect dictionaries
mkdir renamed_SLC_dicts
mkdir marked_SLC_for_aligment
mkdir SLC_align

rm ./marked_SLC_for_aligment/SLC_all.fa
touch ./marked_SLC_for_aligment/SLC_all.fa

Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/DroMel_phylo_rename.R 
/data2/shane/Applications/custom/fasta_rename.py /data2/shane/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa /data2/shane/Documents/SLC_id/SLC_align/updated_Dros_names.csv > /data2/shane/Documents/SLC_id/marked_SLC_for_aligment/Dros_renamed.fa

grep -A 1 "SLC_" /data2/shane/Documents/SLC_id/marked_SLC_for_aligment/Dros_renamed.fa | perl -pe 's/>/>DroMel_/' >> ./marked_SLC_for_aligment/SLC_all.fa
grep -A 1 "SLC_" ./HomSap_Database/reference_proteome/proteome_SLC_mark.fa | perl -pe 's/>/>HomSap_/'  >> ./marked_SLC_for_aligment/SLC_all.fa



## append species name to each SLC name and move to new directory
cat ./general_reference/SLC_info/important_species_codes.txt | while read i
do
  cut -d ',' -f 1 ./final_SLC_dicts/$i* > temp.txt ### Put species prefix on each SLC 
  cat ./final_SLC_dicts/$i* | cut -d ',' -f 2 | sed -e "s/^SLC/$i\_SLC/" | paste --delimiter ',' temp.txt - > ./renamed_SLC_dicts/$i'_Annotated_SLC_dict.csv' ### Put species prefix on each SLC 
  /data2/shane/Applications/custom/fasta_rename.py ./proteomes/$i* ./renamed_SLC_dicts/$i'_Annotated_SLC_dict.csv' | grep -A 1 "SLC_" >> ./marked_SLC_for_aligment/SLC_all.fa
  rm temp.txt
done
find ./renamed_SLC_dicts/* -size 0 -delete

## Divide each SLC family by SLC_###_ marker align and tree

cat ./general_reference/SLC_info/SLC_families.txt | while read i
do
  a=/data2/shane/Documents/SLC_id/SLC_align/$i'.fa'
  grep -A 1 $i ./marked_SLC_for_aligment/SLC_all.fa | sed '/--/d' > $a ## isolate sequences in fasta format
  #/data2/shane/Applications/muscle3.8.31_i86linux64 -in $a -out $a'.aln'
  mafft --threadtb 24 $a > $a'.aln'
  /home/pioannidis/Programs/trimAl/source/trimal -in $a'.aln' -out $a'.aln.trimm'
  #/data2/shane/Applications/trimal-trimAl/source/trimal -in $a'.aln' -out $a'.aln.trimm' ## LOCAL
  /data2/shane/Applications/custom/fasta_2_phylip.sh $a'.aln.trimm' > $a'.aln.trimm.phy'
done

mkdir ./phylogeny/SLC_phylogeny
for i in /data2/shane/Documents/SLC_id/SLC_align/*.phy
do
  b=$(echo $(basename $i) | cut -d '_' -f 1,2) 
  raxfile=$i
  raxdir=/data2/shane/Documents/SLC_id/phylogeny/SLC_phylogeny/
  #rm ./SLC_phylogeny/RAxML*
  /data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T 36 -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
  #/data2/shane/Applications/standard-RAxML-master/raxmlHPC-AVX -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s $raxfile -n $i'.tre' -w $raxdir ## LOCAL
done



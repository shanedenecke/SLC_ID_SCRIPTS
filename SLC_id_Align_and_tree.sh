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
  cp ./real_final_SLC_dicts/$i'_final_SLC_table.csv' ./phylogeny/renamed_dicts/
  #csvcut -c name ./real_final_SLC_tables/$i'_final_SLC_table.csv' | sed 's/"//g' > nam.txt
  #csvcut -c code ./real_final_SLC_tables/$i'_final_SLC_table.csv' | sed 's/"//g' > cod.txt
  num=$(head -1 ./real_final_SLC_dicts/$i'_final_SLC_table.csv' | tr ',' '\n' | cat -n | sed -E 's/\s+/,/g' | sed -E 's/^,//g' | grep "code")
  num2=$(echo $num | cut -d',' -f 1) 
  cut -d',' -f $num2 ./real_final_SLC_dicts/$i'_final_SLC_table.csv' | sed 's/"//g' > cod.txt
  #awk -F "," '{print $code}' ./real_final_SLC_tables/$i'_final_SLC_table.csv' | head
  #paste  -d '_' nam.txt cod.txt | sed -e "s/^/"$i"_/g" | sed -e '1s/.*/name/g' > newcol.txt
  #paste -d ',' newcol.txt cod.txt > ./phylogeny/renamed_dicts/$i'Final_SLC_table.csv'
  if [ $i == 'DroMel' ] || [ $i == 'HomSap' ]
  then
    /data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa'
  else
    /data2/shane/Applications/custom/unigene_fa_sub.sh ./proteomes/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa' || /data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/$i'_unigene.faa' cod.txt > ./phylogeny/SLC_fa/$i'_SLC.faa'
  fi
  /data2/shane/Applications/custom/fasta_rename.py ./phylogeny/SLC_fa/$i'_SLC.faa' ./phylogeny/renamed_dicts/$i'_final_SLC_table.csv' >> ./phylogeny/SLC_fa/combined_renamed.faa
rm cod.txt
done

  
  
  
mkdir ./phylogeny/SLC_byfam  
mkdir ./phylogeny/alignments
mkdir ./phylogeny/trimms
mkdir ./phylogeny/phylip
#cat ./general_reference/SLC_info/SLC_families.txt | while read i
cat ./general_reference/SLC_info/SLC_families.txt | while read i
do
  grep -E -A 1 $i ./phylogeny/SLC_fa/combined_renamed.faa | sed '/--/d' > './phylogeny/SLC_byfam/'$i'phylo_subset.faa'
  mafft --thread 24 './phylogeny/SLC_byfam/'$i'phylo_subset.faa' > './phylogeny/alignments/'$i'phylo_subset.faa.aln'
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
  /data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T 240 -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
  #/data2/shane/Applications/standard-RAxML-master/raxmlHPC-AVX -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s $raxfile -n $i'.tre' -w $raxdir ## LOCAL
done
  







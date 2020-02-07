################### ALIGN AND TREE #####################
cd /data2/shane/Documents/SLC_id

## copy Dromel and Homsap tables to final dictionary file
#Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/SLC_Rename_Dm_SLCs.R > ./final_SLC_dicts/DroMelFinal_SLC_table.csv

#cp ./DroMel_Database/SLC_source_dict.csv ./final_SLC_dicts/DroMelFinal_SLC_table.csv
#cp ./HomSap_Database/SLC_source_dict.csv ./final_SLC_dicts/HomSapFinal_SLC_table.csv

mkdir phylogeny_bayer
mkdir ./phylogeny_bayer/renamed_dicts
mkdir ./phylogeny_bayer/SLC_fa
rm ./phylogeny_bayer/renamed_dicts/*
rm ./phylogeny_bayer/SLC_fa/*
cat ./general_reference/SLC_info/Phylo_list_bayer.txt | while read i
do
  cp ./real_final_SLC_dicts/$i'_final_SLC_table.csv' ./phylogeny_bayer/renamed_dicts/
  #csvcut -c name ./real_final_SLC_tables/$i'_final_SLC_table.csv' | sed 's/"//g' > nam.txt
  #csvcut -c code ./real_final_SLC_tables/$i'_final_SLC_table.csv' | sed 's/"//g' > cod.txt
  num=$(head -1 ./real_final_SLC_dicts/$i'_final_SLC_table.csv' | tr ',' '\n' | cat -n | sed -E 's/\s+/,/g' | sed -E 's/^,//g' | grep "code")
  num2=$(echo $num | cut -d',' -f 1) 
  cut -d',' -f $num2 ./real_final_SLC_dicts/$i'_final_SLC_table.csv' | sed 's/"//g' > cod.txt
  #awk -F "," '{print $code}' ./real_final_SLC_tables/$i'_final_SLC_table.csv' | head
  #paste  -d '_' nam.txt cod.txt | sed -e "s/^/"$i"_/g" | sed -e '1s/.*/name/g' > newcol.txt
  #paste -d ',' newcol.txt cod.txt > ./phylogeny_bayer/renamed_dicts/$i'Final_SLC_table.csv'
  if [ $i == 'DroMel' ] || [ $i == 'HomSap' ]
  then
    /data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/$i'_unigene.faa' cod.txt > ./phylogeny_bayer/SLC_fa/$i'_SLC.faa'
  else
    /data2/shane/Applications/custom/unigene_fa_sub.sh ./proteomes/$i'_unigene.faa' cod.txt > ./phylogeny_bayer/SLC_fa/$i'_SLC.faa' || /data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/$i'_unigene.faa' cod.txt > ./phylogeny_bayer/SLC_fa/$i'_SLC.faa'
  fi
  /data2/shane/Applications/custom/fasta_rename.py ./phylogeny_bayer/SLC_fa/$i'_SLC.faa' ./phylogeny_bayer/renamed_dicts/$i'_final_SLC_table.csv' >> ./phylogeny_bayer/SLC_fa/combined_renamed.faa
rm cod.txt
done

  
  
  
mkdir ./phylogeny_bayer/SLC_byfam  
mkdir ./phylogeny_bayer/alignments
mkdir ./phylogeny_bayer/trimms
mkdir ./phylogeny_bayer/phylip
#cat ./general_reference/SLC_info/SLC_families.txt | while read i
cat ./general_reference/SLC_info/SLC_families.txt | while read i
do
  grep -E -A 1 $i ./phylogeny_bayer/SLC_fa/combined_renamed.faa | sed '/--/d' > './phylogeny_bayer/SLC_byfam/'$i'phylo_subset.faa'
  mafft --thread 24 './phylogeny_bayer/SLC_byfam/'$i'phylo_subset.faa' > './phylogeny_bayer/alignments/'$i'phylo_subset.faa.aln'
  /home/pioannidis/Programs/trimAl/source/trimal -in './phylogeny_bayer/alignments/'$i'phylo_subset.faa.aln' -out './phylogeny_bayer/trimms/'$i'phylo_subset.faa.aln.trimm'
  /data2/shane/Applications/custom/fasta_2_phylip.sh './phylogeny_bayer/trimms/'$i'phylo_subset.faa.aln.trimm' > './phylogeny_bayer/phylip/'$i'phylo_subset.faa.aln.trimm.phy'
done

Rscript ./SLC_id_scripts/Phylip_duplicate.R
  
mkdir SLC_phylogeny_bayer
for i in /data2/shane/Documents/SLC_id/phylogeny_bayer/phylip/*.phy
do
  b=$(echo $(basename $i) | cut -d '_' -f 1,2) 
  raxfile=$i
  raxdir=/data2/shane/Documents/SLC_id/SLC_phylogeny_bayer/
  #rm ./SLC_phylogeny_bayer/RAxML*
  /data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 200 -T 24 -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
  #/data2/shane/Applications/standard-RAxML-master/raxmlHPC-AVX -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s $raxfile -n $i'.tre' -w $raxdir ## LOCAL
done
  







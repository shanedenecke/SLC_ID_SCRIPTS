



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
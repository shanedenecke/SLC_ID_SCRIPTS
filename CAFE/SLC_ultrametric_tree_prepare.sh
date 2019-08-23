cd /data2/shane/Documents/SLC_id

mkdir ultrametric_tree
cd ultrametric_tree

cat ../general_reference/non_model_proteomes/raw_fasta/* > ./arthropod_raw_orthodb.faa


cat ../general_reference/non_model_proteomes/odb10v0_OG2genes.tab | cut -f 2 | cut -d'_' -f 2 > ./odb10_genes.txt
cat ../general_reference/non_model_proteomes/odb10v0_OG2genes.tab | cut -f 2 | cut -d'_' -f 1 > ./odb10_taxid.txt
cat ../general_reference/non_model_proteomes/odb10v0_OG2genes.tab | cut -f 1 > ./odb10_groups.txt

Rscript ../SLC_id_scripts/CAFE/SLC_ultrametric_arthropod_sub.R

python3 ../SLC_id_scripts/CAFE/SLC_ultrameric_1t1_sub.py


## ok one random shiny prep command in here
cat ./proteomes/* ./general_reference/model_proteomes/*.faa > ./shiny_prep/all_proteomes.faa


/data2/shane/Applications/custom/unigene_fa_sub.sh ./arthropod_raw_orthodb.faa ./unicodes_fa_subset.txt > ./Raw_ultrametric.faa
/data2/shane/Applications/custom/fasta_rename.py ./Raw_ultrametric.faa ./renaming_dictionary.csv > ./Renamed_ultrameric.faa



cd /data2/shane/Documents/SLC_id/ultrametric_tree
mkdir OG_Fastas
mkdir OG_aln
mkdir OG_trimm
mkdir OG_phy
cat oto_ish_arthropod_orthologues.txt | while read i
do
grep -A 1 $i ./Renamed_ultrameric.faa | sed '/--/d' > './OG_Fastas/'$i'.faa'
mafft --threadtb 24 './OG_Fastas/'$i'.faa' > './OG_Fastas/'$i'.faa.aln'
/home/pioannidis/Programs/trimAl/source/trimal -in './OG_Fastas/'$i'.faa.aln' -out './OG_Fastas/'$i'.faa.aln.trimm'
/data2/shane/Applications/custom/fasta_2_phylip.sh './OG_Fastas/'$i'.faa.aln.trimm' > './OG_Fastas/'$i'.faa.aln.trimm.phy'
done
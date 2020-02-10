cd $H

mkdir ultrametric_tree
cd ultrametric_tree

cat ../GENERAL_REFERENCE/non_model_proteomes/orthoDB_fasta/* > ./arthropod_raw_orthodb.faa
cat ../GENERAL_REFERENCE/non_model_proteomes/non_orthoDB_fasta/* >> ./arthropod_raw_orthodb.faa

cat ../GENERAL_REFERENCE/CAFE/odb10v0_OG2genes.tab | cut -f 2 | cut -d'_' -f 2 > ./odb10_genes.txt
cat ../GENERAL_REFERENCE/CAFE/odb10v0_OG2genes.tab | cut -f 2 | cut -d'_' -f 1 > ./odb10_taxid.txt
cat ../GENERAL_REFERENCE/CAFE/odb10v0_OG2genes.tab | cut -f 1 > ./odb10_groups.txt

Rscript ../SLC_ID_SCRIPTS/CAFE/SLC_ultrametric_arthropod_sub.R

python3 ../SLC_ID_SCRIPTS/CAFE/SLC_ultrameric_1t1_sub.py


$H/SLC_ID_SCRIPTS/general_scripts/unigene_fa_sub.sh ./arthropod_raw_orthodb.faa ./unicodes_fa_subset.txt > ./Raw_ultrametric.faa
$H/SLC_ID_SCRIPTS/general_scripts/fasta_rename.py ./Raw_ultrametric.faa ./renaming_dictionary.csv > ./Renamed_ultrameric.faa



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
$H/SLC_ID_SCRIPTS/general_scripts/fasta_2_phylip.sh './OG_Fastas/'$i'.faa.aln.trimm' > './OG_Fastas/'$i'.faa.aln.trimm.phy'
done

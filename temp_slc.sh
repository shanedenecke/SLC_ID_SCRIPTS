

cd /data2/shane/Documents/SLC_id/ultrametric_tree
mkdir OG_Fastas
cat oto_ish_arthropod_orthologues.txt | while read i
do
grep -A 1 $i ./Renamed_ultrameric.faa | sed '/--/d' > './OG_Fastas/'$i'.faa'
mafft --threadtb 24 './OG_Fastas/'$i'.faa' > './OG_Fastas/'$i'.faa.aln'
/home/pioannidis/Programs/trimAl/source/trimal -in './OG_Fastas/'$i'.faa.aln' -out './OG_Fastas/'$i'.faa.aln.trimm'
/data2/shane/Applications/custom/fasta_2_phylip.sh './OG_Fastas/'$i'.faa.aln.trimm' > './OG_Fastas/'$i'.faa.aln.trimm.phy'
done



cd /data2/shane/Documents/SLC_id/general_reference/orthoDB_process
rm ./output/*


## copy files 
cp /data/panos/db/odb10/arthropoda/Rawdata/* ./raw_fasta

### CREATE TaxID code names (Manual filter to remove non-insect species)
#ls ./raw_fasta | perl -pe 's/([0-9]+_0).fs/$1/g' > ./reference/insect_taxid_codes.tx
#rm ./reference/insect_taxid_codes_names.txt
#while read i
#do
#grep $i /data/panos/db/odb10/odb10v0_species.tab | cut -f 1,2 | sed 's/_0//g' | grep -v 'virus' |\
#perl -pe 's/(^[0-9]+)\t([A-Z]..)[a-z]+ ([a-z]{3}).+$/$1\t$2$3/' >> ./reference/insect_taxid_codes_names.txt
#done < ./reference/insect_taxid_codes.txt
#################### DO NOT RUN every time
#perl -i -pe 's/(.+[a-z]) ([a-z]+$)/$1_$2/g' ./reference/Taxid_key.tsv 


while IFS=$'\t' read -r col1 col2
do
#echo $col1 
#echo $col2
#done < ./reference/Taxid_key.tsv
cp /data/panos/db/odb10/arthropoda/Rawdata/$col1* ./output/$col2'.fasta'
grep -e "^$col1" ./reference/odb10v0_genes.tab | cut -f 1,3 > ./output/$col2'_OrthoDB_key.txt'
echo -e "code,name" | cat - ./output/$col2'_OrthoDB_key.txt' | perl -pe 's/\t/,/g' > ./output/$col2'_OrthoDB_key_DICT.txt'
/data2/shane/Applications/fasta_rename.py ./output/$col2'.fasta' ./output/$col2'_OrthoDB_key_DICT.txt' > ./output/$col2'_unigene.faa'
done < ./reference/Taxid_key.tsv

cd /data2/shane/Documents/SLC_id
cp /data2/shane/Documents/SLC_id/general_reference/orthoDB_process/output/*_unigene.faa /data2/shane/Documents/SLC_id/proteomes/

cd /data2/shane/Documents/SLC_id/
mkdir proteome_clean
## copy files 

mkdir ./proteome_clean/raw_fasta
cp /data/panos/db/odb10/arthropoda/Rawdata/* ./proteome_clean/raw_fasta

mkdir ./proteome_clean/clean_fasta
while IFS=$'\t' read -r col1 col2
do
cp /data/panos/db/odb10/arthropoda/Rawdata/$col1* ./proteome_clean/clean_fasta/$col2'.fasta'
grep -e "^$col1" ./general_reference/non_model_proteomes/odb10v0_genes.tab | cut -f 1,3 > ./proteome_clean/clean_fasta/$col2'_OrthoDB_key.txt'
echo -e "code,name" | cat - ./proteome_clean/clean_fasta/$col2'_OrthoDB_key.txt' | perl -pe 's/\t/,/g' > ./proteome_clean/clean_fasta/$col2'_OrthoDB_key_DICT.txt'
/data2/shane/Applications/fasta_rename.py ./proteome_clean/clean_fasta/$col2'.fasta' ./proteome_clean/clean_fasta/$col2'_OrthoDB_key_DICT.txt' > ./proteome_clean/clean_fasta/$col2'_unigene.faa'
done < ./general_reference/non_model_proteomes/Taxid_OrthoDB_key.tsv 


#### Remove empty duplicate sequences
find ./proteome_clean/clean_fasta/ -size  0 -print0 | xargs -0 rm --
for i in ./proteome_clean/clean_fasta/*_unigene.faa
do
sed 's/\.//g' $i | sed 's/\*//g' | sed 's/ /_/g' > temp.fasta
awk '/^>/{id=$0;getline;arr[id]=$0}END{for(id in arr)printf("%s\n%s\n",id,arr[id])}' temp.fasta > $i
done

mkdir ./proteome_clean/proteome_check/
cp ./proteome_clean/clean_fasta/* ./proteome_clean/proteome_check/
for i in ./proteome_clean/proteome_check/*
do
makeblastdb -in $i -parse_seqids -dbtype prot
done

echo 'Total numberof proteomes is ' $(ls ./proteome_clean/proteome_check/* | grep -E "*unigene.faa$" | wc -l)
echo 'Numberof GOOD proteomes is ' $(ls ./proteome_clean/proteome_check/* | grep -E "*psq$" | wc -l)

for i in ./proteome_clean/proteome_check/*_unigene.faa
do
sed -i 's/\\//g' $i
done

mkdir proteomes
#cp /data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/clean_fasta/*_unigene.faa /data2/shane/Documents/SLC_id/proteomes/
cp ./proteome_clean/proteome_check/*_unigene.faa ./proteomes/
cp /data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/uniprot_other/* /data2/shane/Documents/SLC_id/proteomes/




################## BIN

############## Make Taxid_OrhoDB_key.tsv file
### CREATE TaxID code names (Manual filter to remove non-insect species)
#ls ./raw_fasta | perl -pe 's/([0-9]+_0).fs/$1/g' > ./reference/insect_taxid_codes.tx
#rm ./reference/insect_taxid_codes_names.txt
#while read i
#do
#grep $i /data/panos/db/odb10/odb10v0_species.tab | cut -f 1,2 | sed 's/_0//g' | grep -v 'virus' |\
#perl -pe 's/(^[0-9]+)\t([A-Z]..)[a-z]+ ([a-z]{3}).+$/$1\t$2$3/' >> ./reference/insect_taxid_codes_names.txt
#done < ./reference/insect_taxid_codes.txt

#perl -pe 's/(.+[A-z])_([a-z]+$)/$1_$2/g' /data2/shane/Documents/SLC_id/general_reference/orthoDB_process/reference/Taxid_key.tsv
#perl -pe 's/(^[0-9]+)\t([A-Z]..)[a-z]+_([a-z]{3}).+$/$1\t$2$3/' /data2/shane/Documents/SLC_id/general_reference/orthoDB_process/reference/Taxid_key.tsv


#################### DO NOT RUN every time

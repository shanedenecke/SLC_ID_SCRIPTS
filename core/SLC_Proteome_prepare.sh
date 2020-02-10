cd $H
mkdir proteome_clean


#Subset all orthoDB gene files for only insect ones 
#while read i; do grep $i ./GENERAL_REFERENCE/non_model_proteomes/keys/odb10v0_genes.tab >> ./GENERAL_REFERENCE/non_model_proteomes/Insect_OrthoDB.genes.tab; done < ./GENERAL_REFERENCE/non_model_proteomes/keys/OrthoDB_taxids.tsv
# R script ./SLC_id_scripts/supp/orthodb_table_prepare.R

mkdir ./proteome_clean/clean_fasta
while IFS=$'\t' read -r col1 col2
do
cp ./GENERAL_REFERENCE/non_model_proteomes/orthoDB_fasta/$col1* ./proteome_clean/clean_fasta/$col2'.fasta'
grep -e "^$col1" ./GENERAL_REFERENCE/non_model_proteomes/keys/odb10v0_genes.tab | cut -f 1,3 > ./proteome_clean/clean_fasta/$col2'_OrthoDB_key.txt'
echo -e "code,name" | cat - ./proteome_clean/clean_fasta/$col2'_OrthoDB_key.txt' | perl -pe 's/\t/,/g' > ./proteome_clean/clean_fasta/$col2'_OrthoDB_key_DICT.txt'
/data2/shane/Applications/custom/fasta_rename_fast_only_exact.py ./proteome_clean/clean_fasta/$col2'.fasta' ./proteome_clean/clean_fasta/$col2'_OrthoDB_key_DICT.txt' > ./proteome_clean/clean_fasta/$col2'_unigene.faa'
done < $SPEC



#### Remove empty or duplicate sequences
find ./proteome_clean/clean_fasta/ -size  0 -print0 | xargs -0 rm --
for i in ./proteome_clean/clean_fasta/*_unigene.faa
do
sed 's/\.//g' $i | sed 's/\*//g' | sed 's/ /_/g' | sed 's/\\//g' | tr '[:lower:]' '[:upper:]' > temp.fasta
awk '/^>/{id=$0;getline;arr[id]=$0}END{for(id in arr)printf("%s\n%s\n",id,arr[id])}' temp.fasta > $i
done
rm temp.fasta

#### makeblastDB for each proteome
for i in ./proteome_clean/clean_fasta/*_unigene.faa
do
makeblastdb -in $i -parse_seqids -dbtype prot
done

### report number of good and total proteomes
echo 'Total numberof proteomes is ' $(ls ./proteome_clean/clean_fasta/* | grep -E "*unigene.faa$" | wc -l)
echo 'Numberof GOOD proteomes is ' $(ls ./proteome_clean/clean_fasta/* | grep -E "*psq$" | wc -l)


mkdir proteomes
mv ./proteome_clean/clean_fasta/*_unigene.faa ./proteomes/
cp ./GENERAL_REFERENCE/non_model_proteomes/non_orthoDB_fasta/* ./proteomes/

rm -rf ./proteome_clean/

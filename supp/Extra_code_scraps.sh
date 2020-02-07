

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


#for i in ./proteome_clean/clean_fasta/*_unigene.faa
#do
#sed -i 's/\\//g' $i
#done
#################### DO NOT RUN every time

cd /data2/shane/Documents/SLC_id
cat ./proteomes/* ./general_reference/model_proteomes/*.faa > ./shiny_prep/all_proteomes.faa
/data2/shane/Applications/custom/unigene_fa_sub.sh ./shiny_prep/all_proteomes.faa ./shiny_prep/slc_codes.txt > ./shiny_prep/SLC_all_raw.faa
/data2/shane/Applications/custom/fasta_rename.py ./shiny_prep/SLC_all_raw.faa ./shiny_prep/Rename_SLC_dict.csv > ./shiny_prep/Renamed_SLC.faa
cd /data2/shane/Documents/SLC_id/


/data2/shane/Applications/custom/unigene_fa_sub.sh ./general_reference/model_proteomes/DroMel_unigene.faa ./Pipeline_compare/Unique_test_codes.txt > ./Pipeline_compare/New_drosophila_SLC.faa

/data2/shane/Applications/custom/tmhmm_filter.sh ./Pipeline_compare/New_drosophila_SLC.faa 4 > ./Pipeline_compare/New_drosophila_SLC_TMM_table.txt


/home/pioannidis/Programs/tmhmm-2.0c/bin/tmhmm  ./Pipeline_compare/New_drosophila_SLC.faa | grep "Number of predicted"

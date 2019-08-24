 cd /data2/shane/Documents/SLC_id

#/data2/shane/Applications/custom/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/Rename_SLC_dict.csv > ./TMHMM_filter/Renamed_unfiltered_SLC.faa
#/data2/shane/Applications/custom/tmhmm_filter.sh ./TMHMM_filter/Renamed_unfiltered_SLC.faa 0 > ./TMHMM_filter/SLC_TMHMM_scores.txt

/home/pioannidis/Programs/tmhmm-2.0c/bin/tmhmm  ./TMHMM_filter/Renamed_unfiltered_SLC.faa | grep "Number of predicted" | perl -pe 's/..(.+) Number of predicted TMHs:\s+(\S+)/$1\t$2/g' > ./TMHMM_filter/SLC_TMHMM_scores.txt
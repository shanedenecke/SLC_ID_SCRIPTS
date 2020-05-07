H=/mnt/disk/shane/Transporter_ID/Comparative_SLC_Arthropod
PHYLO=$H/GENERAL_REFERENCE/CAFE/Phylo_list.txt
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
THREADS=14



cd $H

mkdir BUSCO
mkdir BUSCO/clean_summary

cd BUSCO
for i in ../proteomes/*; do
  b=$(echo $(basename $i) | sed 's/_unigene.faa//g')
  mkdir ./BUSCO/$i
  python3 ~/Applications/busco/bin/busco -c $THREADS --config ~/Applications/busco/config/myconfig.ini -m proteins -i $i -o $b -f -l arthropoda_odb10
  cp $b/*.txt > ./clean_summary/$b'_clean_summary.txt'
done
python3 ~/Applications/Custom_Applications/BUSCO_parse.py -dir ./clean_summary/ > $H/BUSCO/BUSCO_final_summary.tsv


cd $H

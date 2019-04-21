
######################## 6) Iteratively search each species

## create databases for each species
mkdir iterative_database
for i in  /data2/shane/Documents/SLC_id/proteomes/*.faa; do
  c=$(echo $(basename $i) | cut -d '_' -f 1) 
  echo 'Iterative Database '$c
  d=$(ls /data2/shane/Documents/SLC_id/Human_Drosophila_crossref/$c*)
  ./SLC_id_scripts/SLC_Create_HMM_DB.sh $i $d /data2/shane/Documents/SLC_id/iterative_database/'iterative_database_'$c
done


## use iterative databases to recursivley search genomes
mkdir iterative_search
mkdir final_SLC_dicts
for i in  /data2/shane/Documents/SLC_id/proteomes/*.faa; do
e=$(echo $(basename $i) | cut -d '_' -f 1) 
f=$(ls -d /data2/shane/Documents/SLC_id/iterative_database/iterative_database*$e)
echo 'Iterative search '$e
./SLC_id_scripts/SLC_HMM_Search.sh $f $i /data2/shane/Documents/SLC_id/iterative_search/'iterative_search_'$e
cp /data2/shane/Documents/SLC_id/iterative_search/'iterative_search_'$e/final_output/SLC_final_output.csv /data2/shane/Documents/SLC_id/final_SLC_dicts/$e'Final_SLC_table.csv'
done


###################### 7) summarize counts of SLC tables
python3 ./SLC_id_scripts/SLC_summary_count.py

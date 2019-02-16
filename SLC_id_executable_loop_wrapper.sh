#### TESTING GITHUB333333333333333

cd ~/Documents/SLC_id/
### SLC draft summary has untested changes

######################## 0) Download and Clean Sequences (Musca Manduca and Glossina still not working)
python3 ./SLC_id_scripts/download_clean_genomes.py
cp ./general_reference/model_proteomes/HarArm_unigene.faa ./proteomes/
find ./proteomes -type f -empty -delete

######################## 1) Build Human Database
./SLC_id_scripts/Create_SLC_HMM_database_exec.sh ~/Documents/SLC_id/general_reference/model_proteomes/HomSap_unigene.faa ~/Documents/SLC_id/general_reference/SLC_info/HomSap_SLC_dict.csv ~/Documents/SLC_id/Human_HMM_SLC



######################## 2) Bulild good quality Drosophila database
mkdir Drosophila_Database

## search from humans
./SLC_id_scripts/Search_Proteome_with_HMM_profiles.sh ~/Documents/SLC_id/Human_HMM_SLC ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Drosophila_Database/Hs_to_Dm_Search

## Xref with flybase SLC calls
Rscript ./SLC_id_scripts/Flybase_human_SLCxref.R > ./Drosophila_Database/SLC_dict_flybaseXref.csv

# create Drosophila iteration 1 database
./SLC_id_scripts/Create_SLC_HMM_database_exec.sh ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Drosophila_Database/SLC_dict_flybaseXref.csv ~/Documents/SLC_id/Drosophila_Database/Post_Flybase_Xref_database

# Iteratative search Drosophila proteome problem*************************
./SLC_id_scripts/Search_Proteome_with_HMM_profiles.sh ~/Documents/SLC_id/Drosophila_Database/Post_Flybase_Xref_database ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Drosophila_Database/Dm_Final_Iterative

## Make final Drosophila database
./SLC_id_scripts/Create_SLC_HMM_database_exec.sh ~/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa ~/Documents/SLC_id/Drosophila_Database/Dm_Final_Iterative/SLC_dict.csv ~/Documents/SLC_id/Dm_Final_Database

## Rename Dm database with R script
Rscript ./SLC_id_scripts/Rename_Dm_SLCs.R > ./Dm_Final_Database/SLC_dict2.csv
rm ~/Documents/SLC_id/Dm_Final_Database/SLC_dict.csv
mv ~/Documents/SLC_id/Dm_Final_Database/SLC_dict2.csv ~/Documents/SLC_id/Dm_Final_Database/SLC_dict.csv


########################  3) Search species with Human database
mkdir Human_search
for i in ~/Documents/SLC_id/proteomes/*.faa; do
a=$(basename $i) 
./SLC_id_scripts/Search_Proteome_with_HMM_profiles.sh ~/Documents/SLC_id/Human_HMM_SLC $i ~/Documents/SLC_id/Human_search/'HUMAN_'$a
done

########################  4) Search other species with Drosohpila database
mkdir Drosophila_search
for i in ~/Documents/SLC_id/proteomes/*.fa*; do
b=$(echo $(basename $i) | cut -d '_' -f 1) 
./SLC_id_scripts/Search_Proteome_with_HMM_profiles.sh ~/Documents/SLC_id/Dm_Final_Database $i ~/Documents/SLC_id/Drosophila_search/'DROSOPHILA_'$b
done


########################  5) collate searches for each species by Xref with Human and Drosophila searches

Rscript ~/Documents/SLC_id/SLC_id_scripts/crossref_human_dros_searches.R


######################## 6) Iteratively search each species

## create databases for each species
mkdir iterative_database
for i in  ~/Documents/SLC_id/proteomes/*.fa
do
c=$(echo $(basename $i) | cut -d '_' -f 1) 
d=$(ls /home/shanedenecke/Documents/SLC_id/Human_Drosophila_crossref/$c*)
./SLC_id_scripts/Create_SLC_HMM_database_exec.sh $i $d ~/Documents/SLC_id/iterative_database/'iterative_database_'$c
done


## use iterative databases to search genome
mkdir iterative_search
mkdir final_SLC_dicts
for i in  ~/Documents/SLC_id/proteomes/*.fa
do
e=$(echo $(basename $i) | cut -d '_' -f 1) 
f=$(ls -d /home/shanedenecke/Documents/SLC_id/iterative_database/iterative_database*$e/)

./SLC_id_scripts/Search_Proteome_with_HMM_profiles.sh $f $i ~/Documents/SLC_id/iterative_search/'iterative_search_'$e

cp ~/Documents/SLC_id/iterative_search/'iterative_search_'$e/SLC_dict.csv ~/Documents/SLC_id/final_SLC_dicts/$e'Final_SLC_dict.csv'
done


###################### 7) summarize counts of SLC tables
python3 ./SLC_id_scripts/SLC_summary_count.py

###################### 8) Extract sequences from each relevant species. Perform alignment and phylogeny

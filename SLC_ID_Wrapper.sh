#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline
PHYLO=$H/GENERAL_REFERENCE/CAFE/Phylo_list.txt
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
THREADS=14

### 2) Set home directory for pipeline
cd $H

### 3) Get proteomes into folder
mkdir -p proteomes
source ./SLC_ID_SCRIPTS/SLC_Proteome_prepare.sh #Relies on parsing OrthoDB files and copies manually curated proteome files

#### Run SLC_id standalone
nohup ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/SLC_id.sh -proteomes $H/proteomes -busco_thresh 75 -threads $THREADS -outdir $H -metadata ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_reference/Arthropod_species_metadata.tsv &

Rscript ./SLC_ID_SCRIPTS/SLC_id/SLC_Figures.R

### clean up folder
mv DroMel_Database ./intermediate
mv HomSap_Database ./intermediate
mv Drosophila_search ./intermediate
mv Human_search ./intermediate
mv TMHMM_filter ./intermediate/
mv preliminary_SLC_dicts intermediate


###################### 7) Extract sequences from each relevant species. Perform alignment and phylogeny
mkdir SLC_phylogeny
mkdir SLC_phylogeny/raxml_trees
#source ./SLC_ID_SCRIPTS/Align_Tree/SLC_id_Align_and_tree.sh
cp ./GENERAL_REFERENCE/phylo_premade/* ./SLC_phylogeny/raxml_trees
for i in ./SLC_phylogeny/raxml_trees/*.tre; do Rscript ~/Applications/Custom_Applications/ggtree_clean_phylogeny.R -r $i -s $PHYLO -o ./SLC_phylogeny/raxml_trees/; done

############# Build Ultrametric Trees
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees
#source ./SLC_ID_SCRIPTS/CAFE/SLC_Ultrametric_tree_generate ### run to build ultrametric trees. Will take a while
cp ./GENERAL_REFERENCE/CAFE/premade_trees/* ./CAFE/clean_raxml_trees/
Rscript ./SLC_ID_SCRIPTS/CAFE/SLC_CAFE_prep.R
source ./SLC_ID_SCRIPTS/CAFE/SLC_CAFE_run_full.sh
Rscript ./SLC_ID_SCRIPTS/CAFE/SLC_CAFE_figures.R



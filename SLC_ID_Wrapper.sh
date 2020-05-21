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
nohup ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/SLC_id.sh -proteomes $H/proteomes -busco_thresh 75 -threads $THREADS -outdir $H -metadata /mnt/disk/shane/Transporter_ID/SLC_id_pipeline/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv &

### Create figures
Rscript ./SLC_ID_SCRIPTS/SLC_Figures.R

#### PERFORM CAFE

## Make directories
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees

## Create species phylogeneis
#source ./ABC_ID_SCRIPTS/SLC_species_phylogeny.sh
cp ./GENERAL_REFERENCE/CAFE/ultrametric_tree_backup/*.nwk ./CAFE/clean_raxml_trees/

#### Run CAFE
./SLC_ID_SCRIPTS/SLC_CAFE_prep.R
./SLC_ID_SCRIPTS/SLC_CAFE5_run_full.sh
./SLC_ID_SCRIPTS/SLC_CAFE5_figures.R


###################### 7) Extract sequences from each relevant species. Perform alignment and phylogeny
mkdir SLC_phylogeny
mkdir SLC_phylogeny/raxml_trees
#source ./SLC_ID_SCRIPTS/Align_Tree/SLC_id_Align_and_tree.sh
cp ./GENERAL_REFERENCE/phylo_premade/* ./SLC_phylogeny/raxml_trees
for i in ./SLC_phylogeny/raxml_trees/*.tre; do Rscript ~/Applications/Custom_Applications/ggtree_clean_phylogeny.R -r $i -s $PHYLO -o ./SLC_phylogeny/raxml_trees/; done


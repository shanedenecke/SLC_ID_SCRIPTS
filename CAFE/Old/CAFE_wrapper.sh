#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/SLC_id'
cd $H


##################### 8) Prepare Ultrametric Tree
source ./SLC_ID_SCRIPTS/CAFE/Ultrametric_tree_generate.sh #### Still not entirely clean


## RUN CAFE
source ./SLC_ID_SCRIPTS/CAFE/CAFE_run_full.sh
Rscript ./SLC_ID_SCRIPTS/CAFE/CAFE_figures.R $H

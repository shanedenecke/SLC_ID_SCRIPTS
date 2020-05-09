#!/usr/bin/env bash


### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteomes: Path to folder containing one or more Arthropod proteomes that you wish to search
  -busco_thresh: Threshold that you wish to set for BUSCO completeness scores. e.g. 75 means that only proteomes with >75 BUSCO score will be considered for the analysis
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically"
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteomes) PROTEOMES="$2"; shift 2;;
    -busco_thresh) BUSCO_THRESH="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
  esac
done


### Establish directory for scripts and reference files common to this shell script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname "$SCRIPT_DIR")"

### For debugging
PROTEOMES=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/proteomes
BUSCO_THRESH=75
THREADS=14
OUTDIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline
SCRIPT_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_scripts/
SOURCE_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/


#### Start analysis
cd $OUTDIR
mkdir SLC_search_intermediates
mkdir $OUTDIR/filtered_proteomes/
### Perform BUSCO assessment of proteomes
if [ "$BUSCO_THRESH" -gt 0 ]; then
	mkdir BUSCO
	mkdir BUSCO/clean_summary
	mkdir BUSCO/junk
	mkdir BUSCO/target_proteomes

	for i in $PROTEOMES/*; do
	  speciesname=$(echo $(basename $i) | cut -d '_' -f 1)
	  mkdir $OUTDIR/BUSCO/$speciesname
	  busco -c $THREADS --config ~/Applications/busco/config/myconfig.ini -m proteins -i $i -o $speciesname --out_path $OUTDIR/BUSCO/ -l arthropoda_odb10
	  cp $OUTDIR/BUSCO/"$speciesname"/*.txt $OUTDIR/BUSCO/clean_summary/"$speciesname"_clean_summary.txt
	done

	### Parse BUSCO files and filter out crap
	$SCRIPT_DIR/BUSCO_parse.py -dir $OUTDIR/BUSCO/clean_summary/ -thresh $busco_thresh> $OUTDIR/BUSCO/BUSCO_final_summary.tsv
	ls $OUTDIR/BUSCO/ | grep -E '^[[:alpha:]]{6}$' | while read i; do mv $OUTDIR/BUSCO/$i* $OUTDIR/BUSCO/junk/;done ### move all BUSCO junk to junk folder
	cat $OUTDIR/BUSCO/BUSCO_final_summary.tsv | cut -f 6 | tail -n +2 | while read i; do cp $PROTEOMES/"$i"* $OUTDIR/filtered_proteomes/; done ### copy all good proteomes to final folder
else 
	cp $PROTEOMES $OUTDIR/filtered_proteomes/
fi


#### 5) Build Human database
$SCRIPT_DIR/SLC_Create_HMM_DB.sh -proteome $SOURCE_DIR/SLC_id_reference/HomSap_unigene.faa -dict $SOURCE_DIR/SLC_id_reference/HomSap_SLC_dict.csv -out $OUTDIR/HomSap_Database -threads $THREADS


### 6) Bulild good quality Drosophila database
$SCRIPT_DIR/SLC_HMM_Search.sh -database $OUTDIR/HomSap_Database -target $SOURCE_DIR/SLC_id_reference/DroMel_unigene.faa -out $OUTDIR/Hs_to_DroMel_Search -threads $THREADS ## Target Drosophila proteome with human database
$SCRIPT_DIR/SLC_Flybase_human_SLCxref.R $OUTDIR > $OUTDIR/Hs_to_DroMel_Search/SLC_source_dict_flybaseXref.csv ## Xref Human search with known SLCs from Flybase
$SCRIPT_DIR/SLC_Create_HMM_DB.sh -proteome $SOURCE_DIR/SLC_id_reference/DroMel_unigene.faa -dict $OUTDIR/Hs_to_DroMel_Search/SLC_source_dict_flybaseXref.csv -out $OUTDIR/DroMel_Database -threads $THREADS ## create Drosophila database
mv Hs_to_DroMel_Search SLC_search_intermediates

### 7) Search species with Human and Drosophila databases
rm $OUTDIR/filtered_proteomes/HomSap_unigene.faa $OUTDIR/filtered_proteomes/DroMel_unigene.faa ### just in case these ended up in your target proteomes folder somehow
mkdir $OUTDIR/Human_search
mkdir $OUTDIR/Drosophila_search
for i in $OUTDIR/filtered_proteomes/*.faa; do
  target_species=$(echo $(basename $i) | cut -d '_' -f 1) 
  echo 'Target Species is '$target_species
  $SCRIPT_DIR/SLC_HMM_Search.sh -database $OUTDIR/HomSap_Database -target $i -out $OUTDIR/Human_search/'HUMAN_'$target_species -threads $THREADS
  $SCRIPT_DIR/SLC_HMM_Search.sh -database $OUTDIR/DroMel_Database -target $i -out $OUTDIR/Drosophila_search/'DROSOPHILA_'$target_species -threads $THREADS
done
$SCRIPT_DIR/SLC_crossref_human_dros_searches.R ### Crossreference Human and Drosophila searches and output to "preliminary SLC dicts" folder



####################### 6) Post Process SLC tables
$SCRIPT_DIR/SLC_id_Pre_TMHMM.R
tmhmm  ./TMHMM_filter/Renamed_unfiltered_SLC.faa > ./TMHMM_filter/SLC_TMHMM_full.txt
cat ./TMHMM_filter/SLC_TMHMM_full.txt | grep 'Number of predicted' | perl -pe 's/^..([A-z].+) Number of predicted TMHs:\s+(\S)/$1\t$2/g' > ./TMHMM_filter/SLC_TMHMM_scores.txt


##### Clean everything up


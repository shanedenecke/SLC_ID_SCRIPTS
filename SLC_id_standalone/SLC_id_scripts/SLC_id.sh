#!/usr/bin/env bash


### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteomes: Path to folder containing one or more Arthropod proteomes that you wish to search. Protoemes should be labled witht a 6 letter abbreviation followed by '_unigene.faa' e.g. DroMel_unigene.faa for Drosophila melanogaster 
  -busco_thresh: Threshold that you wish to set for BUSCO completeness scores. e.g. 75 means that only proteomes with >75 BUSCO score will be considered for the analysis
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
  -metadata: A tab separated table containign the 'Species_name' and the 'abbreviation' of the targeted species
  example
  ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/SLC_id.sh -proteomes $H/proteomes -busco_thresh 75 -threads $THREADS -outdir $H -metadata /mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
  "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteomes) PROTEOMES="$2"; shift 2;;
    -busco_thresh) BUSCO_THRESH="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -metadata) META="$2"; shift 2;;
  esac
done


### Establish directory for scripts and reference files common to this shell script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname "$SCRIPT_DIR")"

### For debugging
#PROTEOMES=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/proteomes
#BUSCO_THRESH=75
#THREADS=14
#OUTDIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline
#META=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
#SCRIPT_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/
#SOURCE_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/SLC_id_standalone/

echo 'The proteomes folder is '$PROTEOMES
echo 'The BUSCO THRESHOLD is '$BUSCO_THRESH
echo 'The Output directory is '$OUTDIR
echo 'The source directory is '$SOURCE_DIR

#### Start analysis
cd $OUTDIR
mkdir -p SLC_search_intermediates
mkdir -p $OUTDIR/filtered_proteomes/
### Perform BUSCO assessment of proteomes
if [ "$BUSCO_THRESH" -gt 0 ]; then
	mkdir -p BUSCO
	mkdir -p BUSCO/clean_summary
	mkdir -p BUSCO/junk
	mkdir -p BUSCO/target_proteomes
	
	cd BUSCO
	for i in "$PROTEOMES"/*; do
	  echo $i
	  speciesname=$(echo $(basename $i) | cut -d '_' -f 1)
	  mkdir -p $speciesname
	  busco -c $THREADS -f --config $SOURCE_DIR/SLC_id_reference/myconfig.ini -m proteins -i $i -o $speciesname -l arthropoda_odb10
	  cp $OUTDIR/BUSCO/"$speciesname"/*.txt $OUTDIR/BUSCO/clean_summary/"$speciesname"_clean_summary.txt
	done
	cd ..

	### Parse BUSCO files and filter out crap
	$SCRIPT_DIR/BUSCO_parse.py -dir $OUTDIR/BUSCO/clean_summary/ -thresh $BUSCO_THRESH > $OUTDIR/BUSCO/BUSCO_final_summary.tsv
	$SCRIPT_DIR/BUSCO_parse.py -dir $OUTDIR/BUSCO/clean_summary/ > $OUTDIR/BUSCO/BUSCO_final_summary_unfiltered.tsv
	ls $OUTDIR/BUSCO/ | grep -E '^[[:alpha:]]{6}$' | while read i; do mv $OUTDIR/BUSCO/$i* $OUTDIR/BUSCO/junk/;done ### move all BUSCO junk to junk folder
	cat $OUTDIR/BUSCO/BUSCO_final_summary.tsv | cut -f 6 | tail -n +2 | while read i; do cp $PROTEOMES/"$i"* $OUTDIR/filtered_proteomes/; done ### copy all good proteomes to final folder
else 
	cp $PROTEOMES $OUTDIR/filtered_proteomes/
fi


#### 5) Build Human database
$SCRIPT_DIR/SLC_Create_HMM_DB.sh -proteome $SOURCE_DIR/SLC_id_reference/HomSap_unigene.faa -dict $SOURCE_DIR/SLC_id_reference/HOMSAP_SLC_DICT.csv -out $OUTDIR/HomSap_Database -threads $THREADS


### 6) Bulild good quality Drosophila database
$SCRIPT_DIR/SLC_HMM_Search.sh -database $OUTDIR/HomSap_Database -target $SOURCE_DIR/SLC_id_reference/DroMel_unigene.faa -out $OUTDIR/Hs_to_DroMel_Search -threads $THREADS ## Target Drosophila proteome with human database
$SCRIPT_DIR/SLC_Flybase_human_SLCxref.R --out $OUTDIR > $OUTDIR/Hs_to_DroMel_Search/SLC_source_dict_flybaseXref.csv ## Xref Human search with known SLCs from Flybase
$SCRIPT_DIR/SLC_Create_HMM_DB.sh -proteome $SOURCE_DIR/SLC_id_reference/DroMel_unigene.faa -dict $OUTDIR/Hs_to_DroMel_Search/SLC_source_dict_flybaseXref.csv -out $OUTDIR/DroMel_Database -threads $THREADS ## create Drosophila database

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
$SCRIPT_DIR/SLC_id_Pre_TMHMM.R --meta $META
tmhmm  ./TMHMM_filter/Renamed_unfiltered_SLC.faa > ./TMHMM_filter/SLC_TMHMM_full.txt
cat ./TMHMM_filter/SLC_TMHMM_full.txt | grep 'Number of predicted' | perl -pe 's/^..([A-z].+) Number of predicted TMHs:\s+(\S)/$1\t$2/g' > ./TMHMM_filter/SLC_TMHMM_scores.txt
$SCRIPT_DIR/SLC_id_post_TMHMM.R --meta $META

##### Clean everything up
mv $OUTDIR/Hs_to_DroMel_Search $OUTDIR/SLC_search_intermediates
mv $OUTDIR/DroMel_Database $OUTDIR/SLC_search_intermediates
mv $OUTDIR/HomSap_Database $OUTDIR/SLC_search_intermediates
mv $OUTDIR/TMHMM* $OUTDIR/SLC_search_intermediates
mv $OUTDIR/BUSCO $OUTDIR/SLC_search_intermediates
mv $OUTDIR/filtered_proteomes $OUTDIR/SLC_search_intermediates
mv $OUTDIR/preliminary_SLC_dicts $OUTDIR/SLC_search_intermediates
mv $OUTDIR/*_search $OUTDIR/SLC_search_intermediates
cp $OUTDIR/SLC_search_intermediates/BUSCO/BUSCO_final_summary.tsv $OUTDIR/Final_outputs/
cp $OUTDIR/SLC_search_intermediates/BUSCO/BUSCO_final_summary_unfiltered.tsv $OUTDIR/Final_outputs/

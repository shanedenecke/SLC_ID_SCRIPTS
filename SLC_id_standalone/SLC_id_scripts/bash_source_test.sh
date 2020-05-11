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
  -metadata: A tab separated table containign the 'Species_name' and the 'abbreviation' of the targeted species"
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
#SCRIPT_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/
#SOURCE_DIR=/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/SLC_id_standalone/

echo 'The proteomes folder is '$PROTEOMES
echo 'The BUSCO THRESHOLD is '$BUSCO_THRESH
echo 'The Output directory is '$OUTDIR
echo 'The source directory is '$SOURCE_DIR

for i in "$PROTEOMES"/*; do
	  head $i
done


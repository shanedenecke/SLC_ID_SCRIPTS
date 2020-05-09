#!/usr/bin/env bash

### ###################### Create HMM profiles for each SLC family in species ###################

### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteome: Predicted unigene protein dataset of organism
  -dict: SLC dictionary with codes and names
  -out: Name of database to be created. This will output a folder with this name
  -threads: threads "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteome) PROTEOME="$2"; shift 2;;
    -dict) DICT="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -out) OUT="$2"; shift 2;;
  esac
done

echo $PROTEOME
echo $OUT
echo $DICT

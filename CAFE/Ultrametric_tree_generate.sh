mkdir $H/CAFE
cd $H/CAFE
mkdir Ultrametric_tree_CAFE
mkdir $H/CAFE/trees/
### SETUP basic data
DataList='Hemipteran!!!Arachnid!!!Arthropod'
Field_Separator=$IFS
# set comma as internal field separator for the string list
IFS=!!!

## Create ultrametric Trees
cd Ultrametric_tree_CAFE
mkdir odb_files
cp $H/GENERAL_REFERENCE/CAFE/odb10v0_OG2genes.33208.tab ./odb_files/
cp $H/GENERAL_REFERENCE/CAFE/odb10v0_OG2genes.6656.tab ./odb_files/
cp $H/GENERAL_REFERENCE/CAFE/odb_all_metazoa.faa ./odb_files/

##### NEED TO FIX
mkdir species_lists
cp $H/GENERAL_REFERENCE/CAFE/*species_list.txt ./species_lists

for i in $DataList
do
  mkdir $i
  cd $i
  if [ $i = "Arthropod" ] || [ $i = "Arachnid" ]; then 
  #echo 'hello'
       $H/SLC_ID_SCRIPTS/CAFE/find_orthologs_from_mapping_data.pl ../odb_files/odb10v0_OG2genes.33208.tab  ../species_lists/$i'_species_list.txt' > orthology.txt
    else
    #echo 'goodbye'
       $H/SLC_ID_SCRIPTS/CAFE/find_orthologs_from_mapping_data.pl  ../odb_files/odb10v0_OG2genes.6656.tab  ../species_lists/$i'_species_list.txt' > orthology.txt
  fi
  
  
  ## move all .names files to new directory
  mkdir sc_fasta
  mv *names sc_fasta
  cd sc_fasta
  
  #Generate Fasta file with all species in list
  fgrep -A1 -f ../../species_lists/$i'_species_list.txt' ../../odb_files/odb_all_metazoa.faa | grep -v "\-\-" > all_species.f.faa
  
  #Extract fasta sequences
  for x in *names; do $H/SLC_ID_SCRIPTS/CAFE/extract_specific_fasta_seqs_v3.pl $x all_species.f.faa > `basename $x .names`.fs; done 
  mkdir fasta
  mv *.fs ./fasta
  #Move each OG to its own directory
  
  #for x in *fs; do mkdir $x.raxml; mv $x $x.raxml; mv `basename $x .fs`.names $x.raxml; done
  #for x in *raxml; do cd $x; ../create_tree_from_fasta_file.sh.pl `ls *.fs` 7227_0,7070_0,7460_0 16 > create_tree.stdout_stderr 2>&1; cd ../; done &
  
  ##make final phylip file
  #mkdir final_tree
  #rm ./final_tree/all_trimms_horiz.trimm.phy
  #touch ./final_tree/all_trimms_horiz.trimm.phy
  
  ### make individual phylips
  for x in ./fasta/*.fs
  do
  #echo $x
  #done
    b=$(echo $(basename $x) | cut -d'.' -f1)
    
    mkdir $b'_aln'
    #cd $b
    cat $x | perl -pe 's/(^>[0-9]+_0).+$/$1/g' > $b'_aln'/$b'_renamed.fasta'
    mafft --thread 12 $b'_aln'/$b'_renamed.fasta' > $b'_aln'/$b'_renamed.fasta.aln'
    $H/SLC_ID_SCRIPTS/general_scripts/trimAl/source/trimal -in $b'_aln'/$b'_renamed.fasta.aln' -out $b'_aln'/$b'_renamed.fasta.aln.trimm'
    $H/SLC_ID_SCRIPTS/general_scripts/fasta_2_phylip.sh $b'_aln'/$b'_renamed.fasta.aln.trimm' | sed '1d' > $b'_aln'/$b'_renamed.fasta.aln.trimm.phy'
    
    #cat ./final_tree/all_trimms_horiz.trimm.phy > temp_hold.txt
    #paste $b'_aln'/$b'_renamed.fasta.aln.trimm.phy'  temp_hold.txt >> ./final_tree/all_trimms_horiz.trimm.phy
  done
  
  ## RUN R SCRIPT Phylip_merge.R
 Rscript $H/SLC_ID_SCRIPTS/CAFE/Phylip_merge.R
  
  if [ $i = "Arthropod" ] || [ $i = "Arachnid" ]; then
       	#echo 'arth'
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T 24 -m PROTGAMMAAUTO -s Full_species.phy -n $i -o 6239_0
  elif [ $i = 'Lepidopteran' ]; then
	#echo 'lep'       
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T 24 -m PROTGAMMAAUTO -s Full_species.phy -n $i -o 7029_0
  elif [ $i = "Hemipteran" ]; then
 	#echo 'hemi'       
	sed -i 's/J/A/g' Full_species.phy
	sed -i 's/\./A/g' Full_species.phy
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T 24 -m PROTGAMMAAUTO -s Full_species.phy -n $i -o 7227_0
 fi
    cp 'RAxML_bipartitions.'$i $H/CAFE/trees/'raxml_tree_'$i'.tre'
  cd ../../
done

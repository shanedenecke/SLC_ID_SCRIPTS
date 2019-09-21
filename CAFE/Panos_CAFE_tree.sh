#28-Aug-2019
/data2/shane/Documents/Ultrametric_tree_CAFE



#For a CAFE analysis, one has to first build a ML phylogeny of the species he's interested in. Then, an R package will be used to convert this ML phylogeny into a chronogram (aka ultrametric tree). Finally, he will give CAFE (a) this chronogram, and (b) a file containing the number of genes for each gene family, in each species.
#A) Build the chronogram
#PWD: /data2/panos/cafe_SLC/Lepidoptera/
#I obtained the list of species Shane is interested in. He is using binomial names, but I need taxon IDs, which I got from one of ODB's flat files.

#grep "Manduca sexta" /data/panos/db/odb10/odb10v0_species.tab | grep -v virus | cut -f1 >> species.list
#I continued by adding all species one by one (too lazy to script it)! Just be aware that its better to keep the "_0" bit in the end of the taxon ID, even though I'm not sure why it is included in ODB10 (it wasn't there in v9). It makes executing the next steps easier.
#Only keep the OGs corresponding to Arthropoda (taxon ID = 6656); the *genes.tab file contains ALL ODB10 OGs and is therefore, huge! Filtering will greatly speed up the process.
#grep -P "at6656\t" odb10v0_OG2genes.tab | perl -pe 's/at\d+//;' > odb10v0_OG2genes.6656.tab

#Run the orthology analysis using the same script you were using for ODB9 data. As long as the taxon IDs contain the aforementioned "_0" bit this script is happy!

mkdir Arthropod2
cd Arthropod2
../scripts/find_orthologs_from_mapping_data.pl ../odb_files/odb10v0_OG2genes.33208.tab /data2/shane/Documents/Ultrametric_tree_CAFE/species_lists/Arthropod_species_list.txt > orthology.txt

#Move the single-copy OGs to a separate directory

mkdir sc_fasta
mv *names sc_fasta
cd sc_fasta

#Generate the fasta file containing all genes (proteins) from all species included in the species.list file. Beware that the following command might include additional species (e.g. taxon IDs 7227_0 AND 17227_0 - if the latter actually exists), but it's okay since the unnecessary one(s) will never be used. The main reason for running this command is to (greatly) reduce the size of the input all_species.faa file and speed-up the process. If you're not interested in reducing runtime (which would be weird) then it's fine if you just use the unfiltered all_species.faa file.

fgrep -A1 -f ../../species_lists/Arthropod_species_list.txt /data/panos/db/odb10/metazoa/all_species.faa | grep -v "\-\-" > all_species.f.faa

Extract fasta sequences

for x in *names; do ../extract_specific_fasta_seqs_v3.pl $x all_species.f.faa > `basename $x .names`.fs; done 

#Move each OG to its own directory

for x in *fs; do mkdir $x.raxml; mv $x $x.raxml; mv `basename $x .fs`.names $x.raxml; done
#for x in *raxml; do cd $x; ../create_tree_from_fasta_file.sh.pl `ls *.fs` 7227_0,7070_0,7460_0 16 > create_tree.stdout_stderr 2>&1; cd ../; done &
mkdir final_tree
rm ./final_tree/all_trimms_horiz.trimm.phy
touch ./final_tree/all_trimms_horiz.trimm.phy
for i in ./fasta/*.fs
do
#echo $i
#done
  b=$(echo $(basename $i) | cut -d'.' -f1)
  
  mkdir $b'_aln'
  #cd $b
  cat $i | perl -pe 's/(^>[0-9]+_0).+$/$1/g' > $b'_aln'/$b'_renamed.fasta'
  mafft --thread 12 $b'_aln'/$b'_renamed.fasta' > $b'_aln'/$b'_renamed.fasta.aln'
  /home/pioannidis/Programs/trimAl/source/trimal -in $b'_aln'/$b'_renamed.fasta.aln' -out $b'_aln'/$b'_renamed.fasta.aln.trimm'
  /data2/shane/Applications/custom/fasta_2_phylip.sh $b'_aln'/$b'_renamed.fasta.aln.trimm' | sed '1d' > $b'_aln'/$b'_renamed.fasta.aln.trimm.phy'
  
  #cat ./final_tree/all_trimms_horiz.trimm.phy > temp_hold.txt
  #paste $b'_aln'/$b'_renamed.fasta.aln.trimm.phy'  temp_hold.txt >> ./final_tree/all_trimms_horiz.trimm.phy
done
  
## RUN R SCRIPT Phylip_merge.R
cd ..
/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T 24 -m PROTGAMMAAUTO -s 'Full_lepi.phy' -n 'Lepidopteran_species.tre' -o 7227_0,7070_0,7460_0
mv RAxML_bipartitions.Lepidopteran_species.tre /data2/shane/Documents/SLC_id/CAFE/trees/raxml_treeLepi.tre






#Concatenate alignments and convert to the Phylip format
#/data2/shane/Applications/custom/fasta_2_phylip.sh all_trimms.trimm > all_trimms.trimm.phy
#cd final_tree
#all_trimms.trimm
#../filter_trees.pl -p raxml -b 0.01 -o BS_ge_0.01
#cd BS_ge_0.01
#convert.sh `ls *.trimmed` > `ls *.trimmed`.phy
#Run the phylogenetic analysis

raxmlHPC-PTHREADS-AVX -T 24 -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s `ls *.phy` -w `pwd` -n `ls *.phy`.tre -o 7227_0,7070_0,7460_0 > stdout_stderr 2>&1 &

Use R to convert this tree to a chronogram

library(ape)
read.tree("RAxML_bipartitions.concat_genes.fs.aln.trimmed.phy.tre") -> my_tree
chronopl(my_tree, lambda=0.1) -> my_tree2
write.tree(my_tree2, file="tree_chronogram.tre")


B) Run CAFE

Take Shane's file with counts for each family and replace the binomial species names with their taxon IDs. Then filter out the species that are not in the chronogram.

PWD: /data2/panos/cafe_SLC/Arthropoda/sc_fasta/BS_ge_0.01/cafe/
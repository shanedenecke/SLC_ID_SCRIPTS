#cd $H 
#cd CAFE

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/CAFE/Hemiptera_species.txt -ortho_algo Orthofinder -outgroups "DroMel BomMor ApiMel AedAeg " -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/Lepidoptera_species.txt -ortho_algo Orthofinder -outgroups "ApiMel" -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/Diptera_species.txt -ortho_algo OrthoDB -outgroups "BomMor" -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/CAFE/Arthropod_species.txt -ortho_algo Orthofinder -outgroups "CaeEle" -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ./GENERAL _REFERENCE/CAFE/Insect_species.txt -ortho_algo Orthofinder -outgroups "DapPul" -threads 14

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ./GENERAL_REFERENCE/CAFE/Arachnid_species.txt -ortho_algo Orthofinder -outgroups "ApiMel DroMel FraOcc MyzPer DapPul" -threads 14

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ./GENERAL_REFERENCE/CAFE/ArachInsect_species.txt -ortho_algo Orthofinder -outgroups "None" -threads 10

#cp ./*/rax_output/RAxML_bipartitions.*.nwk ./clean_raxml_trees/

cd /data2/shane/Documents/SLC_id/CAFE

### SETUP basic data
DataList='Lepidopteran___Hemipteran___Arachnid___Arthropod'
Field_Separator=$IFS
# set comma as internal field separator for the string list
IFS=___

## Create ultrametric Trees
cd Ultrametric_tree_CAFE

for i in $DataList
do
  mkdir $i
  cd $i
  
  if [ $i = "Arthropod" ] || [ $i = "Arachnid" ]; then 
  #echo 'hello'
       /data2/shane/Documents/SLC_id/SLC_id_scripts/CAFE/find_orthologs_from_mapping_data.pl ../odb_files/odb10v0_OG2genes.33208.tab ../species_lists/$i'_species_list.txt' > orthology.txt
    else
    #echo 'goodbye'
       /data2/shane/Documents/SLC_id/SLC_id_scripts/CAFE/find_orthologs_from_mapping_data.pl ../odb_files/odb10v0_OG2genes.6656.tab ../species_lists/$i'_species_list.txt' > orthology.txt
  fi
  
  
  ## move all .names files to new directory
  mkdir sc_fasta
  mv *names sc_fasta
  cd sc_fasta
  
  #Generate Fasta file with all species in list
  fgrep -A1 -f ../../species_lists/$i'_species_list.txt' ../../odb_files/all_species.faa | grep -v "\-\-" > all_species.f.faa
  
  #Extract fasta sequences
  for x in *names; do /data2/shane/Documents/SLC_id/SLC_id_scripts/CAFE/extract_specific_fasta_seqs_v3.pl $x all_species.f.faa > `basename $x .names`.fs; done 
  mkdir fasta
  mv *.fs ./fasta
  #Move each OG to its own directory
  
  #for x in *fs; do mkdir $x.raxml; mv $x $x.raxml; mv `basename $x .fs`.names $x.raxml; done
  #for x in *raxml; do cd $x; ../create_tree_from_fasta_file.sh.pl `ls *.fs` 7227_0,7070_0,7460_0 16 > create_tree.stdout_stderr 2>&1; cd ../; done &
  
  ##make final phylip file
  mkdir final_tree
  rm ./final_tree/all_trimms_horiz.trimm.phy
  touch ./final_tree/all_trimms_horiz.trimm.phy
  
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
    /home/pioannidis/Programs/trimAl/source/trimal -in $b'_aln'/$b'_renamed.fasta.aln' -out $b'_aln'/$b'_renamed.fasta.aln.trimm'
    /data2/shane/Applications/custom/fasta_2_phylip.sh $b'_aln'/$b'_renamed.fasta.aln.trimm' | sed '1d' > $b'_aln'/$b'_renamed.fasta.aln.trimm.phy'
    
    #cat ./final_tree/all_trimms_horiz.trimm.phy > temp_hold.txt
    #paste $b'_aln'/$b'_renamed.fasta.aln.trimm.phy'  temp_hold.txt >> ./final_tree/all_trimms_horiz.trimm.phy
  done
  
  ## RUN R SCRIPT Phylip_merge.R
 Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/CAFE/Phylip_merge.R
  
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
    cp 'RAxML_bipartitions.'$i /data2/shane/Documents/SLC_id/CAFE/trees/'raxml_tree_'$i'.tre'
  cd ../../
done

cd ../ ## Now puts you at CAFE home directory

Rscript /data2/shane/Documents/SLC_id/SLC_id_scripts/CAFE/CAFE_format_coutns.R

## create CAFE script files

for i in $DataList;
do
  touch ./scripts/$i'_Cafe_Script.cafe'
  echo '#!cafe' >> ./scripts/$i'_Cafe_Script.cafe'
  echo 'load -i /data2/shane/Documents/SLC_id/CAFE/CAFE_tables/'$i'_SLC_CAFE_table.tsv -t 10 -l /data2/shane/Documents/SLC_id/CAFE/logfiles/'$i'_logfile.txt -p 0.01' >> ./scripts/$i'_Cafe_Script.cafe'
  a=$(cat ./trees/$i'_tree_ultrametric.tre')
  echo 'tree '$a >> ./scripts/$i'_Cafe_Script.cafe'
  l=$(cat ./trees/$i'_tree_lambda.txt')
  echo 'lambda -s -t '$l >> ./scripts/$i'_Cafe_Script.cafe'
  echo 'report /data2/shane/Documents/SLC_id/CAFE/outputs/'$i'_cafe_output' >> ./scripts/$i'_Cafe_Script.cafe'
  
  ## run cafe
  cafe ./scripts/$i'_Cafe_Script.cafe'
  python2 ./scripts/Fulton_python_scripts/cafetutorial_report_analysis.py -i outputs/$i'_cafe_output.cafe' -o outputs/$i'_summary.txt'
done

python2 ./scripts/Fulton_python_scripts/cafetutorial_report_analysis.py -i ./outputs/Lepidopteran_cafe_output.cafe -o ./outputs/Lepidopteran_summary.txt
### Arthropod

python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Arth_summary.txt_node.txt \
-t '(((TetUrt:57.81880282,(LimPol:49.71837335,IxoSca:49.71837335):8.100429473):5.959877943,((FolCan:56.07268736,DapPul:56.07268736):3.997986376,((CalSpl:20.20918347,LadFul:20.20918347):33.5137095,((LocMig:41.76599698,ZooNev:41.76599698):9.795154716,(((HalHal:42.30200991,AcyPis:42.30200991):5.261729612,PedHum:47.56373952):2.754776794,((((DanPle:12.92318211,BomMor:12.92318211):28.0946928,(DroMel:28.58784296,AnoGam:28.58784296):12.43003196):3.232169938,(LepDec:21.14607746,TriCas:21.14607746):23.1039674):3.355186797,(ApiMel:23.6822136,NasVit:23.6822136):23.92301805):2.713284662):1.242635384):2.161741264):6.347780776):3.708007023):36.22131924,CaeEle:100)' \
-d '(((TetUrt<0>,(LimPol<2>,IxoSca<4>)<3>)<1>,((FolCan<6>,DapPul<8>)<7>,((CalSpl<10>,LadFul<12>)<11>,((LocMig<14>,ZooNev<16>)<15>,(((HalHal<18>,AcyPis<20>)<19>,PedHum<22>)<21>,((((DanPle<24>,BomMor<26>)<25>,(DroMel<28>,AnoGam<30>)<29>)<27>,(LepDec<32>,TriCas<34>)<33>)<31>,(ApiMel<36>,NasVit<38>)<37>)<35>)<23>)<17>)<13>)<9>)<5>,CaeEle<40>)<39>' \
-y Expansions -o ./outputs/Arthropod_Expansions_tree.png


python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Arth_summary.txt_node.txt \
-t '(((TetUrt:57.81880282,(LimPol:49.71837335,IxoSca:49.71837335):8.100429473):5.959877943,((FolCan:56.07268736,DapPul:56.07268736):3.997986376,((CalSpl:20.20918347,LadFul:20.20918347):33.5137095,((LocMig:41.76599698,ZooNev:41.76599698):9.795154716,(((HalHal:42.30200991,AcyPis:42.30200991):5.261729612,PedHum:47.56373952):2.754776794,((((DanPle:12.92318211,BomMor:12.92318211):28.0946928,(DroMel:28.58784296,AnoGam:28.58784296):12.43003196):3.232169938,(LepDec:21.14607746,TriCas:21.14607746):23.1039674):3.355186797,(ApiMel:23.6822136,NasVit:23.6822136):23.92301805):2.713284662):1.242635384):2.161741264):6.347780776):3.708007023):36.22131924,CaeEle:100)' \
-d '(((TetUrt<0>,(LimPol<2>,IxoSca<4>)<3>)<1>,((FolCan<6>,DapPul<8>)<7>,((CalSpl<10>,LadFul<12>)<11>,((LocMig<14>,ZooNev<16>)<15>,(((HalHal<18>,AcyPis<20>)<19>,PedHum<22>)<21>,((((DanPle<24>,BomMor<26>)<25>,(DroMel<28>,AnoGam<30>)<29>)<27>,(LepDec<32>,TriCas<34>)<33>)<31>,(ApiMel<36>,NasVit<38>)<37>)<35>)<23>)<17>)<13>)<9>)<5>,CaeEle<40>)<39>' \
-y Contractions -o ./outputs/Arthropod_Contractions_tree.png

python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Arth_summary.txt_node.txt \
-t '(((TetUrt:57.81880282,(LimPol:49.71837335,IxoSca:49.71837335):8.100429473):5.959877943,((FolCan:56.07268736,DapPul:56.07268736):3.997986376,((CalSpl:20.20918347,LadFul:20.20918347):33.5137095,((LocMig:41.76599698,ZooNev:41.76599698):9.795154716,(((HalHal:42.30200991,AcyPis:42.30200991):5.261729612,PedHum:47.56373952):2.754776794,((((DanPle:12.92318211,BomMor:12.92318211):28.0946928,(DroMel:28.58784296,AnoGam:28.58784296):12.43003196):3.232169938,(LepDec:21.14607746,TriCas:21.14607746):23.1039674):3.355186797,(ApiMel:23.6822136,NasVit:23.6822136):23.92301805):2.713284662):1.242635384):2.161741264):6.347780776):3.708007023):36.22131924,CaeEle:100)' \
-d '(((TetUrt<0>,(LimPol<2>,IxoSca<4>)<3>)<1>,((FolCan<6>,DapPul<8>)<7>,((CalSpl<10>,LadFul<12>)<11>,((LocMig<14>,ZooNev<16>)<15>,(((HalHal<18>,AcyPis<20>)<19>,PedHum<22>)<21>,((((DanPle<24>,BomMor<26>)<25>,(DroMel<28>,AnoGam<30>)<29>)<27>,(LepDec<32>,TriCas<34>)<33>)<31>,(ApiMel<36>,NasVit<38>)<37>)<35>)<23>)<17>)<13>)<9>)<5>,CaeEle<40>)<39>' \
-y Rapid -o ./outputs/Arthropod_Rapid_tree.png


### Hemipteran

python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Hemi_summary.txt_node.txt \
-t '((FraOcc:95.24237355,(((HomVit:70.0290631,NilLug:70.0290631):9.118605695,(GerBue:61.55875364,((RhoPro:41.9759561,CimLec:41.9759561):8.361451822,(OncFas:30.18902677,HalHal:30.18902677):20.14838116):11.22134572):17.58891515):5.597694879,(BemTab:78.12736443,(DiaCit:71.95282651,((RhoPad:2.567083752,AphGly:2.567083752):1.323422136,(AcyPis:2.747793034,((MyzCer:1.544881883,MyzPer:1.544881883):0.899960938,DiuNox:2.444842821):0.3029502129):1.142712854):68.06232062):6.174537922):6.61799924):10.49700988):4.757626447,(((BomMor:74.58918654,DroMel:74.58918654):10.29767767,TriCas:84.88686421):10.11162783,ApiMel:94.99849204):5.001507956);' \
-d '((FraOcc<0>,(((HomVit<2>,NilLug<4>)<3>,(GerBue<6>,((RhoPro<8>,CimLec<10>)<9>,(OncFas<12>,HalHal<14>)<13>)<11>)<7>)<5>,(BemTab<16>,(DiaCit<18>,((RhoPad<20>,AphGly<22>)<21>,(AcyPis<24>,((MyzCer<26>,MyzPer<28>)<27>,DiuNox<30>)<29>)<25>)<23>)<19>)<17>)<15>)<1>,(((BomMor<32>,DroMel<34>)<33>,TriCas<36>)<35>,ApiMel<38>)<37>)<31>' \
-y Expansions -o ./outputs/Hemipteran_Expansions_tree.png


python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Hemi_summary.txt_node.txt \
-t '((FraOcc:95.24237355,(((HomVit:70.0290631,NilLug:70.0290631):9.118605695,(GerBue:61.55875364,((RhoPro:41.9759561,CimLec:41.9759561):8.361451822,(OncFas:30.18902677,HalHal:30.18902677):20.14838116):11.22134572):17.58891515):5.597694879,(BemTab:78.12736443,(DiaCit:71.95282651,((RhoPad:2.567083752,AphGly:2.567083752):1.323422136,(AcyPis:2.747793034,((MyzCer:1.544881883,MyzPer:1.544881883):0.899960938,DiuNox:2.444842821):0.3029502129):1.142712854):68.06232062):6.174537922):6.61799924):10.49700988):4.757626447,(((BomMor:74.58918654,DroMel:74.58918654):10.29767767,TriCas:84.88686421):10.11162783,ApiMel:94.99849204):5.001507956);' \
-d '((FraOcc<0>,(((HomVit<2>,NilLug<4>)<3>,(GerBue<6>,((RhoPro<8>,CimLec<10>)<9>,(OncFas<12>,HalHal<14>)<13>)<11>)<7>)<5>,(BemTab<16>,(DiaCit<18>,((RhoPad<20>,AphGly<22>)<21>,(AcyPis<24>,((MyzCer<26>,MyzPer<28>)<27>,DiuNox<30>)<29>)<25>)<23>)<19>)<17>)<15>)<1>,(((BomMor<32>,DroMel<34>)<33>,TriCas<36>)<35>,ApiMel<38>)<37>)<31>' \
-y Rapid -o ./outputs/Hemipteran_Rapid_tree.png


python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Hemi_summary.txt_node.txt \
-t '((FraOcc:95.24237355,(((HomVit:70.0290631,NilLug:70.0290631):9.118605695,(GerBue:61.55875364,((RhoPro:41.9759561,CimLec:41.9759561):8.361451822,(OncFas:30.18902677,HalHal:30.18902677):20.14838116):11.22134572):17.58891515):5.597694879,(BemTab:78.12736443,(DiaCit:71.95282651,((RhoPad:2.567083752,AphGly:2.567083752):1.323422136,(AcyPis:2.747793034,((MyzCer:1.544881883,MyzPer:1.544881883):0.899960938,DiuNox:2.444842821):0.3029502129):1.142712854):68.06232062):6.174537922):6.61799924):10.49700988):4.757626447,(((BomMor:74.58918654,DroMel:74.58918654):10.29767767,TriCas:84.88686421):10.11162783,ApiMel:94.99849204):5.001507956);' \
-d '((FraOcc<0>,(((HomVit<2>,NilLug<4>)<3>,(GerBue<6>,((RhoPro<8>,CimLec<10>)<9>,(OncFas<12>,HalHal<14>)<13>)<11>)<7>)<5>,(BemTab<16>,(DiaCit<18>,((RhoPad<20>,AphGly<22>)<21>,(AcyPis<24>,((MyzCer<26>,MyzPer<28>)<27>,DiuNox<30>)<29>)<25>)<23>)<19>)<17>)<15>)<1>,(((BomMor<32>,DroMel<34>)<33>,TriCas<36>)<35>,ApiMel<38>)<37>)<31>' \
-y Contractions -o ./outputs/Hemipteran_Contraction_tree.png


### leipidopteran

python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Lepidopteran_summary.txt_node.txt \
-t '((AcyPis:711.7649398,((FraOcc:585.7314918,LocMig:585.7314918):71.40820331,((TriCas:603.9649631,((PluXyl:253.5494823,(((PapGla:58.22630076,(PapPol:36.61084968,PapXut:36.61084968):21.61545108):145.238105,(PieRap:179.6743168,(DanPle:148.1985816,(MelCin:125.6298448,HelMel:125.6298448):22.56873678):31.47573522):23.79008887):20.0325078,(((ManSex:162.2229771,BomMor:162.2229771):32.59187408,(OpeBru:182.7749885,(HelArm:78.78583246,SpoLit:78.78583246):103.989156):12.03986269):16.81989415,AmyTra:211.6347453):11.86216818):30.05256875):308.3139157,DroMel:561.8633979):42.10156514):34.23023193,ApiMel:638.195195):18.94450016):54.62524461):288.2350602,TetUrt:1000)' \
-d '((AcyPis<0>,((FraOcc<2>,LocMig<4>)<3>,((TriCas<6>,((PluXyl<8>,(((PapGla<10>,(PapPol<12>,PapXut<14>)<13>)<11>,(PieRap<16>,(DanPle<18>,(MelCin<20>,HelMel<22>)<21>)<19>)<17>)<15>,(((ManSex<24>,BomMor<26>)<25>,(OpeBru<28>,(HelArm<30>,SpoLit<32>)<31>)<29>)<27>,AmyTra<34>)<33>)<23>)<9>,DroMel<36>)<35>)<7>,ApiMel<38>)<37>)<5>)<1>,TetUrt<40>)<39>' \
-y Rapid -o ./outputs/Lepidoptera_Rapid_tree.png


python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Lepidopteran_summary.txt_node.txt \
-t '((AcyPis:711.7649398,((FraOcc:585.7314918,LocMig:585.7314918):71.40820331,((TriCas:603.9649631,((PluXyl:253.5494823,(((PapGla:58.22630076,(PapPol:36.61084968,PapXut:36.61084968):21.61545108):145.238105,(PieRap:179.6743168,(DanPle:148.1985816,(MelCin:125.6298448,HelMel:125.6298448):22.56873678):31.47573522):23.79008887):20.0325078,(((ManSex:162.2229771,BomMor:162.2229771):32.59187408,(OpeBru:182.7749885,(HelArm:78.78583246,SpoLit:78.78583246):103.989156):12.03986269):16.81989415,AmyTra:211.6347453):11.86216818):30.05256875):308.3139157,DroMel:561.8633979):42.10156514):34.23023193,ApiMel:638.195195):18.94450016):54.62524461):288.2350602,TetUrt:1000)' \
-d '((AcyPis<0>,((FraOcc<2>,LocMig<4>)<3>,((TriCas<6>,((PluXyl<8>,(((PapGla<10>,(PapPol<12>,PapXut<14>)<13>)<11>,(PieRap<16>,(DanPle<18>,(MelCin<20>,HelMel<22>)<21>)<19>)<17>)<15>,(((ManSex<24>,BomMor<26>)<25>,(OpeBru<28>,(HelArm<30>,SpoLit<32>)<31>)<29>)<27>,AmyTra<34>)<33>)<23>)<9>,DroMel<36>)<35>)<7>,ApiMel<38>)<37>)<5>)<1>,TetUrt<40>)<39>' \
-y Expansions -o ./outputs/Lepidoptera_Expansions_tree.png


python2 ./scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./outputs/Lepidopteran_summary.txt_node.txt \
-t '((AcyPis:711.7649398,((FraOcc:585.7314918,LocMig:585.7314918):71.40820331,((TriCas:603.9649631,((PluXyl:253.5494823,(((PapGla:58.22630076,(PapPol:36.61084968,PapXut:36.61084968):21.61545108):145.238105,(PieRap:179.6743168,(DanPle:148.1985816,(MelCin:125.6298448,HelMel:125.6298448):22.56873678):31.47573522):23.79008887):20.0325078,(((ManSex:162.2229771,BomMor:162.2229771):32.59187408,(OpeBru:182.7749885,(HelArm:78.78583246,SpoLit:78.78583246):103.989156):12.03986269):16.81989415,AmyTra:211.6347453):11.86216818):30.05256875):308.3139157,DroMel:561.8633979):42.10156514):34.23023193,ApiMel:638.195195):18.94450016):54.62524461):288.2350602,TetUrt:1000)' \
-d '((AcyPis<0>,((FraOcc<2>,LocMig<4>)<3>,((TriCas<6>,((PluXyl<8>,(((PapGla<10>,(PapPol<12>,PapXut<14>)<13>)<11>,(PieRap<16>,(DanPle<18>,(MelCin<20>,HelMel<22>)<21>)<19>)<17>)<15>,(((ManSex<24>,BomMor<26>)<25>,(OpeBru<28>,(HelArm<30>,SpoLit<32>)<31>)<29>)<27>,AmyTra<34>)<33>)<23>)<9>,DroMel<36>)<35>)<7>,ApiMel<38>)<37>)<5>)<1>,TetUrt<40>)<39>' \
-y Contractions -o ./outputs/Lepidoptera_Contractions_tree.png

cd $H
THREADS=2



## create CAFE script files

#rm ./CAFE/scripts/*Cafe_Script.*
#rm ./CAFE/outputs/*
#rm ./CAFE/logfiles/*
mkdir ./CAFE/scripts
mkdir ./CAFE/logfiles
mkdir ./CAFE/outputs

for i in ./CAFE/CAFE_tables/*.tsv
do
  b=$(echo $(basename $i) | sed 's/_SLC_CAFE_table.tsv//g')
  
  ##### Make CAFE Script
  rm ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  touch ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  echo '#!cafe' >> ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  echo 'load -i '$H'/CAFE/CAFE_tables/'$b'_SLC_CAFE_table.tsv -t '$THREADS' -l '$H'/CAFE/logfiles/'$b'_SLC_logfile.txt -p 0.05' >> ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  a=$(cat ./CAFE/clean_raxml_trees/$b'_tree_ultrametric.tre')
  echo 'tree '$a >> ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  l=$(cat ./CAFE/clean_raxml_trees/$b'_tree_lambda.txt')
  #echo 'lambda -l '$lambda' -t '$l >> ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  echo 'lambda -s -t '$l >> ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  echo 'report '$H'/CAFE/outputs/'$b'_SLC_cafe_output' >> ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  
  ## run cafe
  /data2/shane/Applications/CAFE/release/cafe ./CAFE/scripts/$b'_SLC_Cafe_Script.cafe'
  python2  $H/SLC_ID_SCRIPTS/CAFE/Fulton_python_scripts/cafetutorial_report_analysis.py -i ./CAFE/outputs/$b'_SLC_cafe_output.cafe' -o ./CAFE/outputs/$b'_SLC_summary.txt'


  t=$(grep 'Tree:' ./CAFE/outputs/$b'_SLC_cafe_output.cafe' | sed 's/Tree://g')
  d=$(grep "IDs of nodes:" ./CAFE/outputs/$b'_SLC_cafe_output.cafe' | sed 's/# IDs of nodes://g')
  
  python2 $H/SLC_ID_SCRIPTS/CAFE/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/$b'_SLC_summary.txt_node.txt' \
  -t $t \
  -d $d \
  -y Expansions -o ./CAFE/outputs/$b'_Expansions_tree.png'
done


############ DRAW TREES


#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Arth_summary.txt_node.txt \
#-t '(((TetUrt:57.81880282,(LimPol:49.71837335,IxoSca:49.71837335):8.100429473):5.959877943,((FolCan:56.07268736,DapPul:56.07268736):3.997986376,((CalSpl:20.20918347,LadFul:20.20918347):33.5137095,((LocMig:41.76599698,ZooNev:41.76599698):9.795154716,(((HalHal:42.30200991,AcyPis:42.30200991):5.261729612,PedHum:47.56373952):2.754776794,((((DanPle:12.92318211,BomMor:12.92318211):28.0946928,(DroMel:28.58784296,AnoGam:28.58784296):12.43003196):3.232169938,(LepDec:21.14607746,TriCas:21.14607746):23.1039674):3.355186797,(ApiMel:23.6822136,NasVit:23.6822136):23.92301805):2.713284662):1.242635384):2.161741264):6.347780776):3.708007023):36.22131924,CaeEle:100)' \
#-d '(((TetUrt<0>,(LimPol<2>,IxoSca<4>)<3>)<1>,((FolCan<6>,DapPul<8>)<7>,((CalSpl<10>,LadFul<12>)<11>,((LocMig<14>,ZooNev<16>)<15>,(((HalHal<18>,AcyPis<20>)<19>,PedHum<22>)<21>,((((DanPle<24>,BomMor<26>)<25>,(DroMel<28>,AnoGam<30>)<29>)<27>,(LepDec<32>,TriCas<34>)<33>)<31>,(ApiMel<36>,NasVit<38>)<37>)<35>)<23>)<17>)<13>)<9>)<5>,CaeEle<40>)<39>' \
#-y Expansions -o ./CAFE/outputs/Arthropod_Expansions_tree.png


#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Arth_summary.txt_node.txt \
#-t '(((TetUrt:57.81880282,(LimPol:49.71837335,IxoSca:49.71837335):8.100429473):5.959877943,((FolCan:56.07268736,DapPul:56.07268736):3.997986376,((CalSpl:20.20918347,LadFul:20.20918347):33.5137095,((LocMig:41.76599698,ZooNev:41.76599698):9.795154716,(((HalHal:42.30200991,AcyPis:42.30200991):5.261729612,PedHum:47.56373952):2.754776794,((((DanPle:12.92318211,BomMor:12.92318211):28.0946928,(DroMel:28.58784296,AnoGam:28.58784296):12.43003196):3.232169938,(LepDec:21.14607746,TriCas:21.14607746):23.1039674):3.355186797,(ApiMel:23.6822136,NasVit:23.6822136):23.92301805):2.713284662):1.242635384):2.161741264):6.347780776):3.708007023):36.22131924,CaeEle:100)' \
#-d '(((TetUrt<0>,(LimPol<2>,IxoSca<4>)<3>)<1>,((FolCan<6>,DapPul<8>)<7>,((CalSpl<10>,LadFul<12>)<11>,((LocMig<14>,ZooNev<16>)<15>,(((HalHal<18>,AcyPis<20>)<19>,PedHum<22>)<21>,((((DanPle<24>,BomMor<26>)<25>,(DroMel<28>,AnoGam<30>)<29>)<27>,(LepDec<32>,TriCas<34>)<33>)<31>,(ApiMel<36>,NasVit<38>)<37>)<35>)<23>)<17>)<13>)<9>)<5>,CaeEle<40>)<39>' \
#-y Contractions -o ./CAFE/outputs/Arthropod_Contractions_tree.png

#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Arth_summary.txt_node.txt \
#-t '(((TetUrt:57.81880282,(LimPol:49.71837335,IxoSca:49.71837335):8.100429473):5.959877943,((FolCan:56.07268736,DapPul:56.07268736):3.997986376,((CalSpl:20.20918347,LadFul:20.20918347):33.5137095,((LocMig:41.76599698,ZooNev:41.76599698):9.795154716,(((HalHal:42.30200991,AcyPis:42.30200991):5.261729612,PedHum:47.56373952):2.754776794,((((DanPle:12.92318211,BomMor:12.92318211):28.0946928,(DroMel:28.58784296,AnoGam:28.58784296):12.43003196):3.232169938,(LepDec:21.14607746,TriCas:21.14607746):23.1039674):3.355186797,(ApiMel:23.6822136,NasVit:23.6822136):23.92301805):2.713284662):1.242635384):2.161741264):6.347780776):3.708007023):36.22131924,CaeEle:100)' \
#-d '(((TetUrt<0>,(LimPol<2>,IxoSca<4>)<3>)<1>,((FolCan<6>,DapPul<8>)<7>,((CalSpl<10>,LadFul<12>)<11>,((LocMig<14>,ZooNev<16>)<15>,(((HalHal<18>,AcyPis<20>)<19>,PedHum<22>)<21>,((((DanPle<24>,BomMor<26>)<25>,(DroMel<28>,AnoGam<30>)<29>)<27>,(LepDec<32>,TriCas<34>)<33>)<31>,(ApiMel<36>,NasVit<38>)<37>)<35>)<23>)<17>)<13>)<9>)<5>,CaeEle<40>)<39>' \
#-y Rapid -o ./CAFE/outputs/Arthropod_Rapid_tree.png


### Hemipteran

#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Hemi_summary.txt_node.txt \
#-t '((FraOcc:95.24237355,(((HomVit:70.0290631,NilLug:70.0290631):9.118605695,(GerBue:61.55875364,((RhoPro:41.9759561,CimLec:41.9759561):8.361451822,(OncFas:30.18902677,HalHal:30.18902677):20.14838116):11.22134572):17.58891515):5.597694879,(BemTab:78.12736443,(DiaCit:71.95282651,((RhoPad:2.567083752,AphGly:2.567083752):1.323422136,(AcyPis:2.747793034,((MyzCer:1.544881883,MyzPer:1.544881883):0.899960938,DiuNox:2.444842821):0.3029502129):1.142712854):68.06232062):6.174537922):6.61799924):10.49700988):4.757626447,(((BomMor:74.58918654,DroMel:74.58918654):10.29767767,TriCas:84.88686421):10.11162783,ApiMel:94.99849204):5.001507956);' \
#-d '((FraOcc<0>,(((HomVit<2>,NilLug<4>)<3>,(GerBue<6>,((RhoPro<8>,CimLec<10>)<9>,(OncFas<12>,HalHal<14>)<13>)<11>)<7>)<5>,(BemTab<16>,(DiaCit<18>,((RhoPad<20>,AphGly<22>)<21>,(AcyPis<24>,((MyzCer<26>,MyzPer<28>)<27>,DiuNox<30>)<29>)<25>)<23>)<19>)<17>)<15>)<1>,(((BomMor<32>,DroMel<34>)<33>,TriCas<36>)<35>,ApiMel<38>)<37>)<31>' \
#-y Expansions -o ./CAFE/outputs/Hemipteran_Expansions_tree.png


#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Hemi_summary.txt_node.txt \
#-t '((FraOcc:95.24237355,(((HomVit:70.0290631,NilLug:70.0290631):9.118605695,(GerBue:61.55875364,((RhoPro:41.9759561,CimLec:41.9759561):8.361451822,(OncFas:30.18902677,HalHal:30.18902677):20.14838116):11.22134572):17.58891515):5.597694879,(BemTab:78.12736443,(DiaCit:71.95282651,((RhoPad:2.567083752,AphGly:2.567083752):1.323422136,(AcyPis:2.747793034,((MyzCer:1.544881883,MyzPer:1.544881883):0.899960938,DiuNox:2.444842821):0.3029502129):1.142712854):68.06232062):6.174537922):6.61799924):10.49700988):4.757626447,(((BomMor:74.58918654,DroMel:74.58918654):10.29767767,TriCas:84.88686421):10.11162783,ApiMel:94.99849204):5.001507956);' \
#-d '((FraOcc<0>,(((HomVit<2>,NilLug<4>)<3>,(GerBue<6>,((RhoPro<8>,CimLec<10>)<9>,(OncFas<12>,HalHal<14>)<13>)<11>)<7>)<5>,(BemTab<16>,(DiaCit<18>,((RhoPad<20>,AphGly<22>)<21>,(AcyPis<24>,((MyzCer<26>,MyzPer<28>)<27>,DiuNox<30>)<29>)<25>)<23>)<19>)<17>)<15>)<1>,(((BomMor<32>,DroMel<34>)<33>,TriCas<36>)<35>,ApiMel<38>)<37>)<31>' \
#-y Rapid -o ./CAFE/outputs/Hemipteran_Rapid_tree.png


#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Hemi_summary.txt_node.txt \
#-t '((FraOcc:95.24237355,(((HomVit:70.0290631,NilLug:70.0290631):9.118605695,(GerBue:61.55875364,((RhoPro:41.9759561,CimLec:41.9759561):8.361451822,(OncFas:30.18902677,HalHal:30.18902677):20.14838116):11.22134572):17.58891515):5.597694879,(BemTab:78.12736443,(DiaCit:71.95282651,((RhoPad:2.567083752,AphGly:2.567083752):1.323422136,(AcyPis:2.747793034,((MyzCer:1.544881883,MyzPer:1.544881883):0.899960938,DiuNox:2.444842821):0.3029502129):1.142712854):68.06232062):6.174537922):6.61799924):10.49700988):4.757626447,(((BomMor:74.58918654,DroMel:74.58918654):10.29767767,TriCas:84.88686421):10.11162783,ApiMel:94.99849204):5.001507956);' \
#-d '((FraOcc<0>,(((HomVit<2>,NilLug<4>)<3>,(GerBue<6>,((RhoPro<8>,CimLec<10>)<9>,(OncFas<12>,HalHal<14>)<13>)<11>)<7>)<5>,(BemTab<16>,(DiaCit<18>,((RhoPad<20>,AphGly<22>)<21>,(AcyPis<24>,((MyzCer<26>,MyzPer<28>)<27>,DiuNox<30>)<29>)<25>)<23>)<19>)<17>)<15>)<1>,(((BomMor<32>,DroMel<34>)<33>,TriCas<36>)<35>,ApiMel<38>)<37>)<31>' \
#-y Contractions -o ./CAFE/outputs/Hemipteran_Contraction_tree.png


### leipidopteran

#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Lepidopteran_summary.txt_node.txt \
#-t '((AcyPis:711.7649398,((FraOcc:585.7314918,LocMig:585.7314918):71.40820331,((TriCas:603.9649631,((PluXyl:253.5494823,(((PapGla:58.22630076,(PapPol:36.61084968,PapXut:36.61084968):21.61545108):145.238105,(PieRap:179.6743168,(DanPle:148.1985816,(MelCin:125.6298448,HelMel:125.6298448):22.56873678):31.47573522):23.79008887):20.0325078,(((ManSex:162.2229771,BomMor:162.2229771):32.59187408,(OpeBru:182.7749885,(HelArm:78.78583246,SpoLit:78.78583246):103.989156):12.03986269):16.81989415,AmyTra:211.6347453):11.86216818):30.05256875):308.3139157,DroMel:561.8633979):42.10156514):34.23023193,ApiMel:638.195195):18.94450016):54.62524461):288.2350602,TetUrt:1000)' \
#-d '((AcyPis<0>,((FraOcc<2>,LocMig<4>)<3>,((TriCas<6>,((PluXyl<8>,(((PapGla<10>,(PapPol<12>,PapXut<14>)<13>)<11>,(PieRap<16>,(DanPle<18>,(MelCin<20>,HelMel<22>)<21>)<19>)<17>)<15>,(((ManSex<24>,BomMor<26>)<25>,(OpeBru<28>,(HelArm<30>,SpoLit<32>)<31>)<29>)<27>,AmyTra<34>)<33>)<23>)<9>,DroMel<36>)<35>)<7>,ApiMel<38>)<37>)<5>)<1>,TetUrt<40>)<39>' \
#-y Rapid -o ./CAFE/outputs/Lepidoptera_Rapid_tree.png


#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Lepidopteran_summary.txt_node.txt \
#-t '((AcyPis:711.7649398,((FraOcc:585.7314918,LocMig:585.7314918):71.40820331,((TriCas:603.9649631,((PluXyl:253.5494823,(((PapGla:58.22630076,(PapPol:36.61084968,PapXut:36.61084968):21.61545108):145.238105,(PieRap:179.6743168,(DanPle:148.1985816,(MelCin:125.6298448,HelMel:125.6298448):22.56873678):31.47573522):23.79008887):20.0325078,(((ManSex:162.2229771,BomMor:162.2229771):32.59187408,(OpeBru:182.7749885,(HelArm:78.78583246,SpoLit:78.78583246):103.989156):12.03986269):16.81989415,AmyTra:211.6347453):11.86216818):30.05256875):308.3139157,DroMel:561.8633979):42.10156514):34.23023193,ApiMel:638.195195):18.94450016):54.62524461):288.2350602,TetUrt:1000)' \
#-d '((AcyPis<0>,((FraOcc<2>,LocMig<4>)<3>,((TriCas<6>,((PluXyl<8>,(((PapGla<10>,(PapPol<12>,PapXut<14>)<13>)<11>,(PieRap<16>,(DanPle<18>,(MelCin<20>,HelMel<22>)<21>)<19>)<17>)<15>,(((ManSex<24>,BomMor<26>)<25>,(OpeBru<28>,(HelArm<30>,SpoLit<32>)<31>)<29>)<27>,AmyTra<34>)<33>)<23>)<9>,DroMel<36>)<35>)<7>,ApiMel<38>)<37>)<5>)<1>,TetUrt<40>)<39>' \
#-y Expansions -o ./CAFE/outputs/Lepidoptera_Expansions_tree.png


#python2 ./CAFE/scripts/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/Lepidopteran_summary.txt_node.txt \
#-t '((AcyPis:711.7649398,((FraOcc:585.7314918,LocMig:585.7314918):71.40820331,((TriCas:603.9649631,((PluXyl:253.5494823,(((PapGla:58.22630076,(PapPol:36.61084968,PapXut:36.61084968):21.61545108):145.238105,(PieRap:179.6743168,(DanPle:148.1985816,(MelCin:125.6298448,HelMel:125.6298448):22.56873678):31.47573522):23.79008887):20.0325078,(((ManSex:162.2229771,BomMor:162.2229771):32.59187408,(OpeBru:182.7749885,(HelArm:78.78583246,SpoLit:78.78583246):103.989156):12.03986269):16.81989415,AmyTra:211.6347453):11.86216818):30.05256875):308.3139157,DroMel:561.8633979):42.10156514):34.23023193,ApiMel:638.195195):18.94450016):54.62524461):288.2350602,TetUrt:1000)' \
#-d '((AcyPis<0>,((FraOcc<2>,LocMig<4>)<3>,((TriCas<6>,((PluXyl<8>,(((PapGla<10>,(PapPol<12>,PapXut<14>)<13>)<11>,(PieRap<16>,(DanPle<18>,(MelCin<20>,HelMel<22>)<21>)<19>)<17>)<15>,(((ManSex<24>,BomMor<26>)<25>,(OpeBru<28>,(HelArm<30>,SpoLit<32>)<31>)<29>)<27>,AmyTra<34>)<33>)<23>)<9>,DroMel<36>)<35>)<7>,ApiMel<38>)<37>)<5>)<1>,TetUrt<40>)<39>' \
#-y Contractions -o ./CAFE/outputs/Lepidoptera_Contractions_tree.png

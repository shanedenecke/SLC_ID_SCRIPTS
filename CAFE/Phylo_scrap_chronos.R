library(data.table)

a=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/Arachnid_summary.txt_anc.txt')
b=a[`Family ID`=='SLC60']

b$`<7>`

b$`SarSca<0>`




((((SarSca<0>,TetUrt<2>)<1>,((IxoSca<4>,(GalOcc<6>,(VarJac<8>,TroMer<10>)<9>)<7>)<5>,(CenScu<12>,(ParTep<14>,SteMim<16>)<15>)<13>)<11>)<3>,((HyaAzt<18>,((LepSal<20>,TigCal<22>)<21>,EurAff<24>)<23>)<19>,(((AcyPis<26>,((BlaGer<28>,BomImp<30>)<29>,(DroMel<32>,TriCas<34>)<33>)<31>)<27>,FolCan<36>)<35>,(DapMag<38>,DapPul<40>)<39>)<37>)<25>)<17>,CaeEle<42>)<41>
  
  library(ape)
tr=read.tree('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/trees/raxml_tree_named_Arthropod.tre')

## hemi 346.079548	232.768338
## thrips and holometabola 334 to 390
## Hexapoda 452 to 509
### Hemiptera + Thysanoptera			384.047421	281.073592
### Diptera			206.278993	107.261303
### Lepidoptera	177.985354	116.448026
### Holometabola	372.428461	317.797467\
### Hexapoda	509.162636	451.538207
### Arthropoda	579.999990	539.090661
### Aparaglossata	102	326.691600	353.049168	301.864480 (hemipterans and holometabola)
## arthropod nematode 800-1000 according to blaire hedgeus book



nodes <- c(); maxes=c()
maxes=c()
mins=c()
if(("CaeEle" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("CaeEle","DroMel")));maxes=c(maxes,1000);mins=c(800)}
if(("AcyPis" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","DroMel")));maxes=c(maxes,353);mins=c(mins,302)}
if(("ApiMel" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ApiMel","DroMel")));maxes=c(maxes,372);mins=c(mins,317)} ## has fossil
if(("AedAeg" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","DroMel")));maxes=c(maxes,206);mins=c(mins,107)} ## has fossil
if(("NilLug" %in% tr$tip.label) & ("AcyPis" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","NilLug")));maxes=c(maxes,346);mins=c(mins,232)}
if(("TetUrt" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("TetUrt","DroMel")));maxes=c(maxes,579);mins=c(mins,539)}
if(("PluXyl" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("PluXyl","BomMor")));maxes=c(maxes,579);mins=c(mins,539)} ## has fossil



mycalibration <- makeChronosCalib(tr, node=c(nodes), age.min=mins,age.max=maxes)
mytimetree <- chronos(tr, lambda = 1, model = "discrete", calibration = mycalibration, control = chronos.control(nb.rate.cat=1))
#mytimetree$edge.length=mytimetree$edge.length*1000
#mytimetree$node.label=NULL
write.tree(mytimetree, file='/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/trees/arth_ult_new.tre')

test=read.tree('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/trees/arth_ult_new.tre')

gp=ggtree(test)#, mrsd = "2010-01-01")
gp=gp+geom_tiplab(size=6)#,face='bold')
gp=gp+geom_nodepoint()
gp=gp+geom_nodelab(hjust=.15)
gp=gp+xlim(0, max(maxes)+100)
gp=gp+theme_tree()
#gp=gp+theme(
print(gp)


## create ultrametric tree
l.tree.ch=chronopl(tr), lambda=0.1)
l.tree.ch$edge.length=l.tree.ch$edge.length*1000
l.tree.ch$node.label=NULL
write.tree(l.tree.ch, file=paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))
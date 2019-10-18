library(data.table)
library(ape)
library(ggplot2)
library(ggtree)

tr=read.tree('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/trees/raxml_tree_named_Arthropod.tre')



nodes <- c(); maxes=c()
maxes=c()
mins=c()
if(("CaeEle" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("CaeEle","DroMel")));maxes=c(maxes,1000);mins=c(800)}
if(("AcyPis" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","DroMel")));maxes=c(maxes,401);mins=c(mins,345)}
if(("ApiMel" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ApiMel","DroMel")));maxes=c(maxes,372);mins=c(mins,317)} ## has fossil
if(("AedAeg" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AedAeg","DroMel")));maxes=c(maxes,206);mins=c(mins,107)} ## has fossil
if(("NilLug" %in% tr$tip.label) & ("AcyPis" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","NilLug")));maxes=c(maxes,346);mins=c(mins,232)}
if(("TetUrt" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("TetUrt","DroMel")));maxes=c(maxes,579);mins=c(mins,539)}
if(("PluXyl" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("PluXyl","BomMor")));maxes=c(maxes,178);mins=c(mins,116)} ## has fossil
    

mycalibration <- makeChronosCalib(tr, node=c(nodes), age.min=mins,age.max=maxes)
mytimetree <- chronos(tr, lambda = 1, model = "correlated", calibration = mycalibration)
write.tree(mytimetree, file='/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/trees/arth_ult_new.tre')
  
test=read.tree('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/OUTPUTS_SEP21_CAFE/outputs/trees/arth_ult_new.tre')

num=test$node.label %>% as.numeric()
cols=c()
for(i in num){
  if(is.na(i)){cols=c(cols,'black')
  }else if(i>90){cols=c(cols,'green4')
  }else if(i>70){cols=c(cols,'blue4')
  }else{cols=c(cols,'red')}
}

ma=max(mytimetree$edge.length)
xma=ma+100
ma.r=seq(0,round(ma,-2),by=100)

diff=ma-round(ma,-2)


gp=ggtree(test)#, mrsd = "2010-01-01")
gp=gp+geom_tiplab(size=6)#,face='bold')
gp=gp+geom_nodepoint(size=2,col=cols)
gp=gp+theme_tree2()
gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
print(gp)

suba=read.tree('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_named_Hemipteran.tre')
a=read.tree('/data2/shane/Documents/SLC_id/CAFE/trees/Hemipteran_tree_ultrametric.tre')
b=fread('/data2/shane/Documents/SLC_id/CAFE/outputs/Hemipteran_SLC_summary.txt_anc.txt')[`Family ID`=='SLC33']

convert='(((PedHum<0>,(((LadFul<2>,CalSpl<4>)<3>,(LocMig<6>,(BlaGer<8>,ZooNev<10>)<9>)<7>)<5>,(((NilLug<12>,HomVit<14>)<13>,(GerBue<16>,((OncFas<18>,HalHal<20>)<19>,(CimLec<22>,RhoPro<24>)<23>)<21>)<17>)<15>,(BemTab<26>,(DiaCit<28>,((AphGly<30>,RhoPad<32>)<31>,(AcyPis<34>,(DiuNox<36>,(MyzCer<38>,MyzPer<40>)<39>)<37>)<35>)<33>)<29>)<27>)<25>)<11>)<1>,DroMel<42>)<41>);'
writeLines(convert,'/data2/shane/Documents/Panos_CAFE_lepi/lab.tree')
lab.tree=read.tree('/data2/shane/Documents/Panos_CAFE_lepi/lab.tree')


c= b %>% select(-matches('[A-z]'))
a$node.label=lab.tree$node.label

d=as.numeric(c)
names(d)=colnames(c)

sort(order(names(a$node.label))[names(d)]) 
sorted=c()
for(i in a$node.label){
  print(i)
  temp=d[names(d)==i]
 sorted=c(sorted,temp)
}
final=c('',as.character(sorted))

u=a
u$node.label=final

gp=ggtree(u)#, mrsd = "2010-01-01")
gp=gp+geom_tiplab(size=6)#,face='bold')
gp=gp+geom_nodepoint(size=2)
gp=gp+geom_nodelab(hjust=0)
gp=gp+theme_tree2()
#gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
print(gp)

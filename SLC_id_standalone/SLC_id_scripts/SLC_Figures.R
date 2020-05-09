#!/usr/bin/env R
library(dplyr)
library(data.table)
library(VennDiagram)
library(svglite)
library(gplots)
library(ggplot2)
library(ggsci)
library(gridExtra)

dir.create('Figures')
setwd('Figures')

co.var=fread('../GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv',header=T) #%>%
  #select(Species_name,abbreviation,Taxonomic_Classification,Phagy,Phagy2,Vory,Diet_category)
full.count=fread('../Final_raw_outputs/TableS7_count_summary.csv')
full.table=fread('../Final_raw_outputs/TableS6_Full_dict_table.csv') ### Need to generate this file
#total.counts=full.table %>% group_by(Species_name) %>% summarize(fam_size=length(Species_name)) %>% data.table() ##  probably can delete



################ FIGURE 2
dir.create('Pipeline_compare')
transporter.db=fread('../GENERAL_REFERENCE/model_SLC_info/Dm_Transporter_DB_manual.csv')
flybase=fread('../GENERAL_REFERENCE/model_SLC_info/DroMel_SLC_table_flybase.csv')

dros.slcs=fread('../real_final_SLC_dicts/DroMel_final_SLC_table.csv')
dros.slcs$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",dros.slcs$code)
#test=fread('./Dm_Database_Generate/DroMel_iterative_search/final_output/total_slc_table.csv')
#test$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",test$code)

l=list(transporter.db$CG,flybase$CG,dros.slcs$CG)
names(l)=c('TransporterDB','Flybase','This Study \n(Denecke et. al 2020)')
ven=venn.diagram(l,filename=NULL,cex=1.5,fill=c('coral','cadetblue','gold'))
ggsave(ven,file='./Figure2_Pipeline_compare.svg',device='svg',width=10,height=10,units='cm')

uni.shane=setdiff(dros.slcs$CG,transporter.db$CG) %>% setdiff(flybase$CG)
uni.tdb=setdiff(transporter.db$CG,dros.slcs$CG) %>% setdiff(flybase$CG)

shane.uni=dros.slcs[CG %in% uni.shane]
fwrite(shane.uni,'./Pipeline_compare/Unique_to_our_study.csv')

tdb.uni=transporter.db[CG %in% uni.tdb]
fwrite(tdb.uni,'./Pipeline_compare/Unique_to_transporter_DB.csv')

#################### FIGURE 3

## Figure 3


#sub.size=total.counts[Species_name==sp.input.mod]$fam_size

gp=ggplot(full.count,aes(x=SLC_total))
gp=gp+geom_histogram(colour="black", fill="grey75",binwidth=20)
gp=gp+geom_density(alpha=.2, fill="#FF6666")
gp=gp+labs(x='\nTotal SLCs Identified in Species',y='Frequency\n')
gp=gp+theme_bw()
gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
            axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
            axis.title=element_text(size=22),axis.text.x=element_text(size=18),
            legend.position = 'none',plot.title = element_text(hjust = 0.5))
#print(gp)

ggsave(gp,file='./FigureS2_histogram.pdf',device='pdf',width=20,height=10,units='cm')

##################### FIGURE 4 #####################
counts.summary=full.count
#counts.summary$SLC_62=NULL
m=merge(counts.summary,co.var,by='abbreviation') %>% filter(abbreviation!='HomSap') %>% data.table()


groups=c('Hymenoptera','Coleoptera','Hemiptera','Lepidoptera','Diptera','Arachnida','Crustacea')
cols=c('firebrick2','blue4','magenta','green3','orange','mediumorchid3','gold3')
#cols=pal_aaas('default')(7)

names(groups)=cols
final.cols=c()
for(i in 1:nrow(m)){
  g=m$Taxonomic_Classification[i]
  
  if(g %in% groups){
  final.cols[i]=names(groups[which(groups==g)])
  }else{
    final.cols[i]='black'
  }
}
counts.matrix=m %>% 
  select(matches("SLC"),-matches('Unsorted'),-matches('Unsorted'),-SLC_total) %>%
  as.matrix() %>% t()
colnames(counts.matrix)=m$Species_name[m$Species_name!='Homo_sapiens']

#par(mar=c(1,4,10,3)) 
pdf('Figure3_heatmap.pdf',width=20,height=10)
heatmap.2(counts.matrix,Rowv=F,Colv=T,scale="row",col=colorpanel(75,'blue','grey','red'),dendrogram = 'column',tracecol=NA,colCol = final.cols,margins = c(10,5),cexRow=.7,
          density.info = 'density',denscol='black')
dev.off()



################# FIGURE 4 ##################33

count.sum=select(full.count,-matches('Unsorted'),-SLC_total)

############### filter out SLC families which have <100 total members in dataset 
v=c()
for(i in colnames(count.sum)[2:length(count.sum)]){
  sub=as.numeric(count.sum[[i]])
  #if(max(sub)<10){# & sub[173]<.5){
  if(sum(sub)<100 | max(sub)<10){
    v=c(v,i)
    count.sum[[i]]=NULL
  }
}


## merge to olympia's data
meta=merge(count.sum,co.var,by="abbreviation")

## set varaibles for loop
comps=c('Taxonomic_Classification','Phagy','Vory')
l=list()

### Perform loop which goes through each SL family and comparison and builds linear model for each
for(i in colnames(meta)[grep('SLC_',colnames(meta))]){
  for(j in comps){
    
    ##Filter for infrequent categories
    good=names(which(table(meta[[j]])>5)) ## filter for occurances which don't occur at least 5 times
    sub=meta[which(meta[[j]] %in% good)] ## subset for comparisons which have at least 5 occurances
    
    if(nrow(sub)>0){ #### catch cases where no comparison available
      if(length(unique(sub[[j]]))>1){
    ##anova
      model=aov(formula=sub[[i]]~sub[[j]])
      pval=summary(model)[[1]][["Pr(>F)"]][1]
      ef=max(abs(model$coefficients[!grepl('Intercept',names(model$coefficients))]))
      l[[paste(i,j,sep='_')]]=data.table(family=i,co_variable=j,pval=pval,max_effect=ef) 
      }
    }
  }
}

summary=rbindlist(l)  
summary$bonf=p.adjust(summary$pval,method='bonferroni') 

#tab.one=summary %>% arrange(bonf) %>% filter(bonf<1e-5)  %>% filter(max_effect>4) %>% data.table()
tab.one=summary %>% arrange(bonf)  %>% filter(max_effect>4) %>% data.table()


fwrite(summary,'./ANOVA_full.csv')
fwrite(tab.one,'./Table1_ANOVA_significant.csv')


dir.create('ANOVA_meta_plots')

plot.table=tab.one %>% distinct(family,co_variable) %>% data.table()


  
for(i in 1:nrow(plot.table)){
  row=plot.table[i]
  co=row$co_variable
  fam=row$family
  red=select(m,fam,co)
  red[[fam]]=as.numeric(red[[fam]])
  
  plot=red[red[[co]] %in% names(which(table(red[[co]])>5))] ### remove elements that aren't present at least 5 times
  
  ylab=paste0('Number of ',gsub('_',' ',fam),' Members','\n')
  tit=paste0(gsub('_',' ',fam),' Family Size vrs. ',gsub('_',' ',co))
  xlab=paste0('\n',gsub('_',' ',co))
  
  gp=ggplot(plot,aes_string(co,y=fam,fill=co))
  gp=gp+geom_boxplot(outlier.size=1)
  gp=gp+labs(x=xlab,y=ylab)
  gp=gp+scale_fill_rickandmorty()
  gp=gp+ggtitle(tit)
  gp=gp+theme_bw()
  gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
              axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
              strip.text=element_text(size=20),strip.background=element_rect("white"),
              axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1),
              legend.position = 'none',plot.title = element_text(hjust = 0.5))
  ggsave(paste0('./ANOVA_meta_plots/',co,'_',fam,'.pdf'),device='pdf',plot=gp)
  
  
  #print(gp)
  
  gp2=gp+labs(x='',y='')
  gp2=gp2+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
                   axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
                   axis.text.x=element_blank(),title=element_text(size=10),
                   legend.position='none',plot.title = element_text(hjust = 0.5))
  gp2=gp2+ggtitle(fam)
  
  assign(paste0(co,'.',fam),gp2)
}

grid.plot=grid.arrange(Taxonomic_Classification.SLC_36,
             Taxonomic_Classification.SLC_22,
             Taxonomic_Classification.SLC_2,
             Taxonomic_Classification.SLC_35,
             Taxonomic_Classification.SLC_60,
             Taxonomic_Classification.SLC_33,
             nrow=2)

ggsave(grid.plot,file='Figure4_ANOVA_comparisons.pdf',device='pdf',width=20,height=15,units='cm')



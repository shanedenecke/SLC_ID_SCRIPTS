#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(ggplot2))
shhh(library(ggsci))
shhh(library(stringi))
shhh(library(seqinr))


######################## functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}

shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}
########################################3



#setwd('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id') #local
setwd('/data2/shane/Documents/SLC_id')
dir.create('TMHMM_filter')
dir.create('Final_raw_outputs')


## read in data
meta.data=fread('./general_reference/Co_variables/Arthropod_species_metadata.csv',header=T)
human.hmm=fread('./general_reference/SLC_info/Human_SLC_HMM.txt')
colnames(human.hmm)=c('code','tm_domains')

comp.score=fread('./genome_score/comp_score.txt') %>% 
  filter(grepl('DROSOPHILA',score.name)) %>% separate(score.name,into=c('dros','abbreviation'),'_') %>% 
  data.table()


dros.hmm=fread('./general_reference/SLC_info/Drosophila_Flybase_SLC_TMHMM.csv')
colnames(dros.hmm)=c('code','tm_domains','family')

slc_fams=readLines('/data2/shane/Documents/SLC_id/general_reference/SLC_info/SLC_families.txt')
slc_fams=sapply(slc_fams,dash.remove)
names(slc_fams)=NULL



### copy DroMel and HomSap databases to final_dicts directory
file.remove('/data2/shane/Documents/SLC_id/final_SLC_dicts/DroMelFinal_SLC_table.csv')
file.remove('/data2/shane/Documents/SLC_id/final_SLC_dicts/HomSapFinal_SLC_table.csv')
file.copy('./Dm_Database_Generate/SLC_source_dict_flybaseXref.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/DroMelFinal_SLC_table.csv')
file.copy('./HomSap_Database/SLC_source_dict.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/HomSapFinal_SLC_table.csv')

file.copy('./general_reference/model_proteomes/DroMel_unigene.faa','./proteomes/')
file.copy('./general_reference/model_proteomes/HomSap_unigene.faa','./proteomes/')




###Generate HMM filter tables
human.hmm$family=sapply(human.hmm$code,dash.remove)
model.hmm=rbindlist(list(human.hmm,dros.hmm))
model.hmm.key=human.hmm %>% group_by(family) %>% summarize(minimum=max(0,min((min(tm_domains)-2),(min(tm_domains)/2)))) %>% data.table()
model.hmm.key=rbindlist(list(model.hmm.key,data.table(family='SLC_Unsorted',minimum=0)),use.names = T)


########## GET PRELIMINARY LIST AND SUBSET FASTA

#i='./final_SLC_dicts/CaeEleFinal_SLC_table.csv'
l=list()
for(i in list.files('./final_SLC_dicts',full.names = T)){
  dict=fread(i)
  abbrev=gsub('./final_SLC_dicts/','',i,fixed=T)
  abbrev=gsub('Final_SLC_table.csv','',abbrev,fixed=T)
  fam=sapply(dict$name,dash.remove)
  dict$abbreviation=abbrev
  dict$family=fam
  dict2=merge(dict,select(meta.data,Species_name,abbreviation),by='abbreviation',use.names=T) %>% mutate(name)
  #dict2=dict
  dict2$name=paste(dict2$Species_name,dict2$abbreviation,dict2$code,dict2$name,sep='___')
  dict2$name=gsub(' ','_',dict2$name)
  l[[i]]=dict2
  ind.rename.dict= dict2 %>% select(name,code)
  writeLines(ind.rename.dict$code,'./TMHMM_filter/slc_unfiltered_codes.txt')
  fwrite(ind.rename.dict,'./TMHMM_filter/temp_rename_dict.csv')
  file.copy(paste0('./proteomes/',abbrev,'_unigene.faa'),'./TMHMM_filter/temp_proteome.faa',overwrite = T)
  
  #system('
        
  #      /data2/shane/Applications/custom/unigene_fa_sub.sh ./TMHMM_filter/temp_proteome.faa ./TMHMM_filter/slc_unfiltered_codes.txt > ./TMHMM_filter/SLC_unfiltered_all_raw.faa 
  #      /data2/shane/Applications/custom/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/temp_rename_dict.csv >> ./TMHMM_filter/Renamed_unfiltered_SLC.faa
  #        ')
}


#system("
#       /home/pioannidis/Programs/tmhmm-2.0c/bin/tmhmm  ./TMHMM_filter/Renamed_unfiltered_SLC.faa | grep 'Number of predicted' | perl -pe 's/..(.+) Number of predicted TMHs:\s+(\S+)/$1\t$2/g' > ./TMHMM_filter/SLC_TMHMM_scores.txt
#       ")

unnamed.fasta=read.fasta('./TMHMM_filter/Renamed_unfiltered_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)


all=rbindlist(l,use.names = T)

arth.hmm=fread('./TMHMM_filter/SLC_TMHMM_scores.txt') 
#arth.hmm=fread('./TMHMM_filter/Cae_SLC_TMHMM_scores.txt',header=F,sep='\t') 
colnames(arth.hmm)=c('name','tm_count')
arth.hmm=arth.hmm[grepl('___',name)]

m=merge(arth.hmm,all,by='name') %>% merge(meta.data,by=c('Species_name','abbreviation'))
#m=merge(arth.hmm,all,by='name') 

g.l=list()
b.l=list()
for(i in 1:nrow(model.hmm.key)){
  row=model.hmm.key[i]
  fam=row$family
  mini=row$minimum
  
  sub=m[family==fam]
  
  good=sub[tm_count>=mini]
  #good=sub[tm_count>=mini | tm_count>4]
  g.l[[i]]=good
  
  bad=sub[tm_count<mini & tm_count<=4]
  b.l[[i]]=bad
}
tmhmm.filtered.full=rbindlist(g.l)
tmhmm.removed=rbindlist(b.l)

count.summary=tmhmm.filtered.full %>% group_by(family,abbreviation) %>% summarize(count=length(family)) %>% spread(key=family,value=count) %>% data.table()
count.summary[is.na(count.summary)]=0
#select(count.summary,matches('SLC')) %>% as.matrix()  

count.summary$SLC_total=rowSums(select(count.summary,matches('SLC')))

writeLines(tmhmm.filtered.full$name,'./TMHMM_filter/Filtered_SLC_codes.txt')

system('
  /data2/shane/Applications/custom/unigene_fa_sub.sh ./TMHMM_filter/Renamed_unfiltered_SLC.faa ./TMHMM_filter/Filtered_SLC_codes.txt > ./TMHMM_filter/SLC_filtered_all_raw.faa 
  python3 /data2/shane/Applications/custom/unigene_fa_sub_update.py ./TMHMM_filter/Renamed_unfiltered_SLC.faa ./TMHMM_filter/Filtered_SLC_codes.txt > ./TMHMM_filter/SLC_filtered_all_raw2.faa 
     
       ')



###### Now filter out those with very bad genome completeness scores


## so you don't have to run the same script over and over again
## count.summary=fread('./Final_raw_outputs/count_summary.csv')
## tmhmm.filtered.full=fread('./Final_raw_outputs/Full_dict_table.csv')

comp.score2=comp.score[!duplicated(abbreviation)] 

comp.score3=merge(comp.score2,count.summary,by='abbreviation')

with(comp.score3,cor.test(comp.score,SLC_total))

pdf('./Figures/FigureS1_length_filter_out.pdf',width=20,height=10)
with(comp.score3,plot(comp.score~SLC_total,xlab='Total Number of SLCs',ylab='Fracion of proteins too short',cex=1.5,cex.lab=1.5)) ### make figure
dev.off()

comp.score4=comp.score3[comp.score<.3]

good_spec=c(comp.score4$abbreviation,'DroMel','HomSap')
bad_spec=comp.score3[comp.score>=.3]$abbreviation
  
count.summary2=count.summary[abbreviation %in% good_spec]
tmhmm.filtered.full2=tmhmm.filtered.full[abbreviation %in% good_spec]

fa=read.fasta('./TMHMM_filter/SLC_filtered_all_raw.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)
fa.final=fa[!grepl(paste0(bad_spec,collapse='|'),names(fa))]


##write totals to file
fwrite(count.summary2,'./Final_raw_outputs/TableS4_count_summary.csv')
fwrite(tmhmm.removed,'./Final_raw_outputs/TMM_filtered_out.csv')
fwrite(tmhmm.filtered.full2,'./Final_raw_outputs/TableS3_Full_dict_table.csv')
write.fasta(fa.final,names=names(fa.final),file.out='./Final_raw_outputs/FileS1_All_SLC_final.faa',nbchar=10000,as.string=T)
#save(fa.final,file='./Final_raw_outputs//All_SLCs_fasta.Robj')
  

a=fread('./Final_raw_outputs/TableS3_Full_dict_table.csv') %>% select(abbreviation,code,family) %>% mutate(name=paste0(abbreviation,'__',code,'__',family,'_')) %>% 
  select(-family) %>% data.table()

lapply(split(a,a$abbreviation),function(x) fwrite(select(x,-abbreviation),file=paste0('./real_final_SLC_tables/',x$abbreviation[1],'_final_SLC_table.csv')))





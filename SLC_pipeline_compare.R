library(dplyr)
library(data.table)
library(VennDiagram)
library(venn)
library(svglite)

setwd('/data2/shane/Documents/SLC_id')
dir.create('Pipeline_compare')
transporter.db=fread('./general_reference/SLC_info/Dm_Transporter_DB_manual.csv')
flybase=fread('./general_reference/SLC_info/DroMel_SLC_table_flybase.csv')

shane=fread('./Dm_Database_Generate/Hs_to_DroMel_Search/final_output/total_slc_table.csv')
shane$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",shane$code)
#test=fread('./Dm_Database_Generate/DroMel_iterative_search/final_output/total_slc_table.csv')
#test$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",test$code)

l=list(transporter.db$CG,flybase$CG,shane$CG)
names(l)=c('TransporterDB','Flybase','This Study \n(Denecke et. al 2020)')
#names(l)=c('TransporterDB','Flybase','shane')
ven=venn.diagram(l,filename=NULL,cex=1.5,fill=c('coral','cadetblue','gold'))
ggsave(ven,file='./Pipeline_compare/Figure2_Pipeline_compare.svg',device='svg',width=10,height=10,units='cm')




####
transporter.db[CG %in% setdiff(l$TransporterDB,l$shane)]
unique.shane=shane[CG %in% (setdiff(setdiff(l$shane,l$TransporterDB),l$Flybase))]
unique.shane.codes=shane[CG %in% (setdiff(setdiff(l$shane,l$TransporterDB),l$Flybase))]$code

fwrite(unique.shane,'./Pipeline_compare/Unique_test_table.csv')
writeLines(unique.shane.codes,'./Pipeline_compare/Unique_test_codes.txt')


transporter.db[CG %in% setdiff(l$TransporterDB,l$initial)] %>% fwrite('initial_missed.csv')
unique.initial=initial[CG %in% (setdiff(setdiff(l$initial,l$TransporterDB),l$Flybase))]





m=length(shane$CG)




length(l$TransporterDB)=m
length(l$Flybase)=m
a=cbind(l$TransporterDB,l$Flybase,l$shane) %>% data.table()
fwrite(a,'Venn_test.csv')
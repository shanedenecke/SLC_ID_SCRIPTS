library(dplyr)
library(data.table)
library(VennDiagram)
library(venn)

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
venn(l,zcolor = c('red','blue','green'))
l2=l[1:3]
venn.diagram(l2,'./Pipeline_compare/test.tiff')


transporter.db[CG %in% setdiff(l$TransporterDB,l$Test)]
unique.test=test[CG %in% (setdiff(setdiff(l$Test,l$TransporterDB),l$Flybase))]
unique.test.codes=test[CG %in% (setdiff(setdiff(l$Test,l$TransporterDB),l$Flybase))]$code

fwrite(unique.test,'./Pipeline_compare/Unique_test_table.csv')
writeLines(unique.test.codes,'./Pipeline_compare/Unique_test_codes.txt')


transporter.db[CG %in% setdiff(l$TransporterDB,l$initial)] %>% fwrite('initial_missed.csv')
unique.initial=initial[CG %in% (setdiff(setdiff(l$initial,l$TransporterDB),l$Flybase))]



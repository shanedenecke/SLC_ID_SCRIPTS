library(dplyr)
library(data.table)
library(VennDiagram)


setwd('/data2/shane/Documents/SLC_id')
dir.create('Pipeline_compare')
transporter.db=fread('./general_reference/SLC_info/Dm_Transporter_DB_manual.csv')
flybase=fread('./general_reference/SLC_info/DroMel_SLC_table_flybase.csv')

test=fread('./Dm_Database_Generate/DroMel_iterative_search/final_output/total_slc_table.csv')
test$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",test$code)

l=list(transporter.db$CG,flybase$CG,test$CG)
names(l)=c('TransporterDB','Flybase','Test')

l2=l[1:3]
venn.diagram(l2,'./Pipeline_compare/test.tiff')


transporter.db[CG %in% setdiff(l$TransporterDB,l$Test)]
unique.test=test[CG %in% (setdiff(setdiff(l$Test,l$TransporterDB),l$Flybase))]
unique.test.codes=test[CG %in% (setdiff(setdiff(l$Test,l$TransporterDB),l$Flybase))]$code

fwrite(unique.test,'./Pipeline_compare/Unique_test_table.csv')
writeLines(unique.test.codes,'./Pipeline_compare/Unique_test_codes.txt')



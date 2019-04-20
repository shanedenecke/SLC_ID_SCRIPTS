library(data.table)

raw=fread('/home/shanedenecke/Documents/omics_data/Bayer_transcriptomes_clean/Helicoverpa_SPKM.csv')

a=readLines('/home/shanedenecke/Downloads/RAxML_bipartitions.SLC_26.tre')

slc.dict=fread('/home/shanedenecke/Documents/SLC_id/final_SLC_dicts/HelArmFinal_SLC_table.csv')
for(file in list.files('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_expression_analysis/trees',pattern='.tre',full.names =T))  
  for(i in 1:length(b$name)){
    raw.tree=readLines(file)
    nam=paste0("HelArm_",slc.dict$name[i])
    cod=slc.dict$code[i]
    new.tree=gsub(nam,cod,raw.tree)
    writeLines(new.tree,file)
  }

writeLines(a,'/home/shanedenecke/Dropbox/wp7_prodrug/SLC_expression_analysis/trees/RAxML_bipartitions.SLC_26___2.tre')

### compare 2 tables for missing values


it=fread('/home/shanedenecke/Documents/SLC_id/Drosophila_Database/Dm_iterative_search/final_output/SLC_final_output.csv')
ori=fread('/home/shanedenecke/Documents/SLC_id/Drosophila_Database/Hs_to_Dm_Database/SLC_source_dict.csv')

m=merge(it,ori,by='code',all=T)
colnames(m)=c('code','iterative','original')

View(m)

m %>% filter(is.na(iterative) | is.na(original))


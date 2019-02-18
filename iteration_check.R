### compare 2 tables for missing values


it=fread('/home/shanedenecke/Documents/SLC_id/iterative_search/iterative_search_BomMor/final_output/SLC_final_output.csv',col.names = c('code','name'))
ori=fread('/home/shanedenecke/Documents/SLC_id/iterative_database/iterative_database_BomMor/SLC_source_dict.csv',col.names = c('code','name'))

it=a
colnames(it)=c('code','name')
m=merge(it,ori,by='code',all=T)
colnames(m)=c('code','iterative','original')
m=m %>% arrange(iterative)
#View(m)

sum=m %>% filter(is.na(iterative) | is.na(original))
sum

View(sum)

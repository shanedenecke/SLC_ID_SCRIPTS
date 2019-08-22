library(venneuler)
library(data.table)
library(dplyr)

a=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/Venn_test.csv')
colnames(a)=c('TDB','Flybase','Shane')

na_replace=function(x){
  v=c()
  for(i in x){
    if(grepl('CG',i)){
      v=c(v,i)
    }
  }
  return(v)
}


all=intersect(a$TDB,a$Flybase) %>% intersect(a$Shane) %>% na_replace() %>% length()

t.only=setdiff(a$TDB,a$Flybase) %>% setdiff(a$Shane) %>% na_replace() %>% length()
s.only=setdiff(a$Shane,a$Flybase) %>% setdiff(a$TDB) %>% na_replace() %>% length()
f.only=setdiff(a$Flybase,a$Shane) %>% setdiff(a$TDB) %>% na_replace() %>% length()

tf.only=intersect(a$TDB,a$Flybase) %>% setdiff(a$Shane) %>% na_replace() %>% length()
sf.only=intersect(a$Shane,a$Flybase) %>% setdiff(a$TDB) %>% na_replace() %>% length()
ts.only=intersect(a$Shane,a$TDB) %>% setdiff(a$Flybase) %>% na_replace() %>% length()



plot(venneuler(c(A=t.only, B=s.only, C=f.only, 'A&B'=ts.only, 'A&C'=tf.only, 'B&C'=sf.only, 'A&B&C'=all)))


plot(venneuler(a))

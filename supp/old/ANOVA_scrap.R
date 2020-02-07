



tukey=TukeyHSD(model)

tab=tukey[[1]] %>% data.table


#minsig=min(tab[["p adj"]])
#effect=max(abs(tukey$diff))
comparison=rownames(tukey[[1]])
#sig=tab$`p adj`
#sig=summary(model)[[1]][["Pr(>F)"]][1]
l[[paste(i,j,sep='_')]]=data.table(tab,family=i,co_variable=j,comp=comparison) 

##lm
#model=lm(formula=sub[[i]]~sub[[j]])
#coef=summary(model)$coefficients 
#comparison=rownames(coef)
#l[[paste(i,j,sep='_')]]= coef %>% data.table(comp=comparison,family=i,co_variable=j) 
}}

## merge all lists into large data table with all associations
#raw.model=rbindlist(l) %>% select(Estimate,`Pr(>|t|)`,comp,family,co_variable) %>% filter(!grepl("Intercept",comp)) %>% data.table() 
#raw.model=rbindlist(l) %>% select(Estimate,`Pr(>|t|)`,family,co_variable) %>% data.table() 
#
#colnames(raw.model)=gsub('Pr(>|t|)','p.value',colnames(raw.model),fixed=T) 
#raw.model$comp=gsub('sub[[j]]','',raw.model$comp,fixed=T)

## Add adjusted estimate and p value corrections
#for(i in 1:nrow(raw.model)){raw.model$adjusted_Estimate[i]=raw.model[i]$Estimate/(size.adjust[fam==raw.model[i]$family])$avg}


raw.model=rbindlist(l)
raw.model$bonf=p.adjust(raw.model[['p adj']],method='bonferroni') 


final=arrange(raw.model,bonf) %>% filter(bonf<1e-10) %>% filter(abs(diff)>5) %>% data.table()

final=arrange(raw.model,bonf) %>% filter(bonf<1e-5) %>% filter(abs(diff)>5) %>% filter(co_variable!='Taxanomic_Classification') %>%
  data.table()



fwrite(raw.model,'./SLC_family_counts/full_correlation_anlaysis.csv')
fwrite(final,'./SLC_family_counts/significant_correlation_anlaysis.csv')


























































################################# COUNTS ANALYSIS

aa=slc.function$AA
sugar=stri_remove_empty(slc.function$Sugar)
drug=stri_remove_empty(slc.function$Drug)
ion=stri_remove_empty(slc.function$Ion)

anal.counts=count.summary %>% select(-slc_total,-SLC_X,-SLC_Unsorted)

anal.counts$SLC_AA=anal.counts %>% select(aa) %>% rowSums()
anal.counts$SLC_sugar=anal.counts %>% select(sugar) %>% rowSums()
anal.counts$SLC_drug=anal.counts %>% select(drug) %>% rowSums()
anal.counts$SLC_ion=anal.counts %>% select(ion) %>% rowSums()




############### filter out SLC families which have <100 total members in dataset 
for(i in colnames(anal.counts)[2:70]){
  sub=as.numeric(anal.counts[[i]])
  #if(max(sub)<10){# & sub[173]<.5){
  if(sum(sub)<100 & max(sub)<10){
    print(i)
    anal.counts[[i]]=NULL
  }
}

mean.fam.size=apply(select(anal.counts,matches('SLC')),2,mean)
size.adjust=data.table(fam=names(mean.fam.size),avg=mean.fam.size)


## merge to olympia's data
full.counts=merge(anal.counts,co.variables,by="abbreviation")


## set varaibles for loop
comps=c('Taxanomic_Classification','Diet_category','Phagy','Vory')
l=list()

### Perform loop which goes through each SL family and comparison and builds linear model for each
for(i in colnames(full.counts)[grep('SLC_',colnames(full.counts))]){
  for(j in comps){
    
    ##Filter for infrequent categories
    good=names(which(table(full.counts[[j]])>5)) ## filter for occurances which don't occur at least 5 times
    sub=full.counts[which(full.counts[[j]] %in% good)] ## subset for comparisons which have at least 5 occurances
    
    ##anova
    model=aov(formula=sub[[i]]~sub[[j]])
    tukey=TukeyHSD(model)
    
    tab=tukey[[1]] %>% data.table
    
    
    #minsig=min(tab[["p adj"]])
    #effect=max(abs(tukey$diff))
    comparison=rownames(tukey[[1]])
    #sig=tab$`p adj`
    #sig=summary(model)[[1]][["Pr(>F)"]][1]
    l[[paste(i,j,sep='_')]]=data.table(tab,family=i,co_variable=j,comp=comparison) 
    
    ##lm
    #model=lm(formula=sub[[i]]~sub[[j]])
    #coef=summary(model)$coefficients 
    #comparison=rownames(coef)
    #l[[paste(i,j,sep='_')]]= coef %>% data.table(comp=comparison,family=i,co_variable=j) 
  }}

## merge all lists into large data table with all associations
#raw.model=rbindlist(l) %>% select(Estimate,`Pr(>|t|)`,comp,family,co_variable) %>% filter(!grepl("Intercept",comp)) %>% data.table() 
#raw.model=rbindlist(l) %>% select(Estimate,`Pr(>|t|)`,family,co_variable) %>% data.table() 
#
#colnames(raw.model)=gsub('Pr(>|t|)','p.value',colnames(raw.model),fixed=T) 
#raw.model$comp=gsub('sub[[j]]','',raw.model$comp,fixed=T)

## Add adjusted estimate and p value corrections
#for(i in 1:nrow(raw.model)){raw.model$adjusted_Estimate[i]=raw.model[i]$Estimate/(size.adjust[fam==raw.model[i]$family])$avg}


raw.model=rbindlist(l)
raw.model$bonf=p.adjust(raw.model[['p adj']],method='bonferroni') 


final=arrange(raw.model,bonf) %>% filter(bonf<1e-10) %>% filter(abs(diff)>5) %>% data.table()

final=arrange(raw.model,bonf) %>% filter(bonf<1e-5) %>% filter(abs(diff)>5) %>% filter(co_variable!='Taxanomic_Classification') %>%
  data.table()


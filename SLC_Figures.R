library(dplyr)
library(data.table)
library(VennDiagram)
library(venn)
library(svglite)
library(gplots)

setwd('/data2/shane/Documents/SLC_id')
dir.create('Figures')
setwd('Figures')




################ FIGURE 2

dir.create('Pipeline_compare')
transporter.db=fread('../general_reference/SLC_info/Dm_Transporter_DB_manual.csv')
flybase=fread('../general_reference/SLC_info/DroMel_SLC_table_flybase.csv')

dros.slcs=fread('./Dm_Database_Generate/Hs_to_DroMel_Search/final_output/total_slc_table.csv')
dros.slcs$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",dros.slcs$code)
#test=fread('./Dm_Database_Generate/DroMel_iterative_search/final_output/total_slc_table.csv')
#test$CG=gsub('FBgn.+_.+_(CG[0-9]+)_.+$',"\\1",test$code)

l=list(transporter.db$CG,flybase$CG,dros.slcs$CG)
names(l)=c('TransporterDB','Flybase','This Study \n(Denecke et. al 2020)')
ven=venn.diagram(l,filename=NULL,cex=1.5,fill=c('coral','cadetblue','gold'))
ggsave(ven,file='./Pipeline_compare/Figure2_Pipeline_compare.svg',device='svg',width=10,height=10,units='cm')


#################### FIGURE 3

full.table=fread('../shiny_prep/Reference_csv_dictionary.csv')

## Figure 3A

total.counts=full.table %>% group_by(Species_name) %>% summarize(fam_size=length(Species_name)) %>% data.table()
sub.size=total.counts[Species_name==sp.input.mod]$fam_size

gp=ggplot(total.counts,aes(x=fam_size))
gp=gp+geom_histogram(colour="black", fill="white",binwidth=20)
gp=gp+geom_density(alpha=.2, fill="#FF6666")
gp=gp+labs(x='Family Size',y='Number of Occurances')
gp=gp+theme_bw()
gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
            axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
            strip.text=element_text(size=20),strip.background=element_rect("white"),
            axis.title=element_text(size=25),axis.text.x=element_text(size=20),
            legend.position = 'none',plot.title = element_text(hjust = 0.5))
print(gp)

ggsave(gp,file='./Figure3A_histogram.svg',device='svg',width=10,height=10,units='cm')

#### FIGURE 3B
counts.summary=fread('/data2/shane/Documents/SLC_id/SLC_family_counts/count_summary.csv') 
  

counts.matrix=counts.summary %>% select(-abbreviation,-SLC_Unsorted,-SLC_X,-slc_total) %>% as.matrix() %>% t()
colnames(counts.matrix)=counts.summary$abbreviation

pdf('heatmap.pdf',width=20)
heatmap.2(counts.matrix,Rowv=F,Colv=T,scale="row",dendrogram = 'column',tracecol=NA)
dev.off()



heatmap.2(matrix,
          #cellnote = matrix, # same data set for cell labels
          main="Heatmap",
          #notecol="black",
          # change font color of cell labels to
          black
          density.info="none", # turns off density plot inside color
          legend
          trace="none",
          # turns off trace lines inside the heat
          map
          margins =c(15,18),
          # widens margins around plot
          col=my_palette,
          # use on color palette defined earlier
          #breaks=col_breaks,
          # enable color transition at specified
          limits
          dendrogram="row",
          # only draw a row dendrogram
          Colv="NA",
          # turn off column clustering
          keysize = 1,
          key.title = NA,
          key.xlab = "Cohen's D",
          #notecex = 0.75,
          srtCol = 90,
          cexCol = 1.5,
          cexRow = 1,labRow = parse(text=rownames(matrix)),
          labCol = parse(text=colnames(matrix)),
          RowSideColors = c(
            # grouping row-variables into different
            rep(geno_colour[1], 3),
            rep(geno_colour[2], 3),
            rep(geno_colour[3], 3),
            rep(geno_colour[4], 3),
            rep(geno_colour[5], 3),
            rep(geno_colour[6], 3),
            rep(geno_colour[7], 3))
)
#par(lend = 1)
legend(.8,.6,
       # location of the legend on the heatmap plot
       legend = parse(text=unique(merge2$gene)), # category labels
       col = geno_colour, # color key
       lty= 1,
       # line style
       lwd = 5,
       # line width
       ncol=1,
       text.width=0.05,
       cex=1,
       bty = "n")
dev.off()
}








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

setwd("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/functional_analysis")
set.seed(123)

library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(siggenes)
library(tidyverse)
library(enrichplot)
library(mygene)
library(pathview)
library(dplyr)

my.cluster<-read.csv("cluster.csv", header = T)[,1]


data.prot<-read.csv("merge_prot.csv", header = T)
data.gene<-read.csv("merge_gene.csv", header = T)
data.metab<-read.csv("merge_metab.csv", header = T)
prot_name<-read.table("prot_name.txt", header = F, sep = "\t", fill = T)
metab_name<-read.table("metab_name.txt", header = F, sep = " ", fill = T)

name.prot<-prot_name[,2]
name.gene<-data.gene[,1]
name.metab<-metab_name[,1]

#data.prot<-data.prot[,-1]
data.gene<-data.gene[,-1]
data.metab<-data.metab[,-1]

data.metab<-apply(data.metab, 2, as.numeric) %>% t()
data.prot<-apply(data.prot, 2, as.numeric) %>% t()

#---------------------PAIRWISE GSEA-------------------------------#

my.data<-data.gene[,my.cluster %in% c("3", "2") %>% which() ]
my.name<-name.gene
na_index<-which(is.na(my.data[1,]))
my.data<-my.data[,-na_index]
my.cluster.sam<-my.cluster[my.cluster %in% c("3","2") %>% which()][-na_index]

my.data<-apply(my.data, 2, as.numeric)

sam.out.all <- sam(my.data, my.cluster.sam, 
                   rand = 123, gene.names = my.name,
                   method = wilc.stat)
sum.sam.out.all <- summary(sam.out.all, findDelta(sam.out.all, 0.01)[1,1])
plot(sam.out.all, findDelta(sam.out.all, 0.01)[1,1])

#write.csv(sum.sam.out.all@mat.sig, "/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/sam1vs2v3.csv", row.names=T)


pval.12<-NULL
for(i in 1:nrow(my.data)){
  w<-wilcox.test(my.data[i,which(my.cluster.sam==3)], 
                 my.data[i,which(my.cluster.sam==2)])
  pval.12<-c(pval.12, w$p.value)
}

#Set FDR cutoff and get genes 
signif <- sum.sam.out.all@mat.sig



#get entrez ID
signif.entrez <- unique(signif %>% rownames())

#signif.entrez <- my.name[which(p.adjust(pval.12) <= 0.01)] %>% unique()

data.gene.gsea<-my.data[which(my.name %in% signif.entrez),]

UA<-data.gene.gsea[,which(my.cluster.sam==3)] %>% rowMeans(na.rm = T)
UB<-data.gene.gsea[,which(my.cluster.sam==2)] %>% rowMeans(na.rm = T)
SigA<-apply(data.gene.gsea[,which(my.cluster.sam==3)], 1, sd, na.rm=T)
SigB<-apply(data.gene.gsea[,which(my.cluster.sam==2)], 1, sd, na.rm=T)
fc<- ((UA-UB)/(SigA+SigB)) #%>% abs()

FC<-data.frame(Gene=name.gene[which(my.name %in% signif.entrez)],
               mean_cg = fc)

FC<-FC[order(FC$mean_cg, decreasing = T),]

##format for gsea 
FC.vec <- FC$mean_cg
names(FC.vec) <- FC$Gene

FC.vec[grep("MTO", names(FC.vec))]

#set score type
scoreType <- "std"

#Get gene set database
C5 <- msigdbr(species = "Homo sapiens", category = "C5")
C6 <- msigdbr(species = "Homo sapiens", category = "C6")
C5.entrez <- dplyr::select(C5, gs_name, gene_symbol)
C6.entrez <- dplyr::select(C6, gs_name, gene_symbol)


C5.ensembl.ls <- C5 %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

C6.ensembl.ls <- C6 %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

names(C5.ensembl.ls)<-gsub("_", " ", names(C5.ensembl.ls))
names(C6.ensembl.ls)<-gsub("_", " ", names(C6.ensembl.ls))

#Run GSEA

#------------GSEA with ClusterProfiler---------------------------------------#

C6_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, gene_symbol)

gsea.C6.cp<-GSEA(FC.vec,
                 exponent = 1,
                 minGSSize = 50,
                 maxGSSize = 500,
                 eps = 0,
                 pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 TERM2GENE = C6_t2g)

#format results
ridgeplot(gsea.C6.cp)
gseaplot2(gsea.C6.cp, title="C6 enrichment",
          geneSetID = 1:3)

dotplot(gsea.C6.cp,orderBy ="Count",
        x="Count", 
        showCategory=20) + ggtitle("dotplot for C6")



C5_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)
head(C5_t2g)
gsea.C5.cp<-GSEA(FC.vec,
                 exponent = 1,
                 minGSSize = 50,
                 maxGSSize = 500,
                 eps = 0,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 TERM2GENE = C5_t2g)


#format results
ridgeplot(gsea.C5.cp)

gseaplot2(gsea.C5.cp, title="C5 enrichment",
          geneSetID = 1:3)
dotplot(gsea.C5.cp, orderBy ="Count",
        x="Count",
        showCategory=20) + ggtitle("dotplot for C5")



#------------------------KEGG--------------------------------------------#
FC.enz<-queryMany(FC.vec %>% names(), scopes="symbol", fields="entrezgene", species="human")
FC.enz<-data.frame(gene_symbol = FC.enz$query,
                   enz_ID = FC.enz$entrezgene)
FC.vec.df<-data.frame(gene_symbol = FC.vec %>% names(),
                      Fc = FC.vec)
FC.vec.enz<-merge(FC.vec.df, FC.enz, "gene_symbol", all=T)
FC.vec.enz<-FC.vec.enz[order(FC.vec.enz$Fc, decreasing = T),]
FC.vec.enz<-FC.vec.enz[-which(FC.vec.enz$enz_ID %>% is.na()),]

KEGG.vec<-FC.vec.enz$Fc
names(KEGG.vec)<-FC.vec.enz$enz_ID

saveRDS(KEGG.vec, "kegg3v2.Rdata")

kk <- gseKEGG(gene         = KEGG.vec,
              organism     = 'hsa',
              minGSSize    = 3,
              maxGSSize    = 500,
              eps          = 0,
              pvalueCutoff = 0.05,
              verbose      = T,
              keyType = "kegg")


hsa04110 <- pathview(gene.data  = KEGG.vec,
                     pathway.id = "hsa04060",
                     species    = "hsa",
                     limit      = list(gene=max(abs(KEGG.vec)), cpd=1))

#----------------------Integrated Pathway-------------------------------------#

p.meta<-NULL
for(i in 1:ncol(data.metab)){
  wm<-wilcox.test(data.metab[my.cluster==1, i] %>% na.omit(),
                  data.metab[my.cluster==2, i] %>% na.omit())
  p.meta<-c(p.meta, (wm$p.value))
}
write.table(name.metab[which(p.meta<0.05)], "metab1vs2.txt", row.names = F)


p.prot<-NULL
for(i in 1:ncol(data.prot)){
  wm<-wilcox.test(data.prot[my.cluster==1, i] %>% na.omit(),
                  data.prot[my.cluster==2, i] %>% na.omit())
  p.prot<-c(p.prot, (wm$p.value))
}
write.table(name.prot[which(p.prot<0.05)], "protein1vs2.txt", row.names = F)

write.table(KEGG.vec %>% names(), "gene1vs2.txt", row.names = F)

#-----------------------------------------------------------------------------#
#-----------------Over Represenation----------------#

my.data<-data.gene
my.name<-name.gene

na_index<-which(is.na(my.data[1,]))
my.data<-my.data[,-na_index]
my.cluster.sam<-my.cluster[-na_index]
my.data<-apply(my.data, 2, as.numeric)

pval.123<-NULL
for(i in 1:nrow(my.data)){
  k<-kruskal.test(my.data[i,], 
                  my.cluster.sam %>% as.factor())
  pval.123<-c(pval.123, k$p.value)
}


#### Enrichment ####

#get entrez ID
#signif.entrez <- unique(signif %>% rownames())
signif.entrez <- my.name[which(p.adjust(pval.123) <= 0.01)] %>% unique()

data.gene.gsea<-my.data[which(my.name %in% signif.entrez),]


#run enrichment
enrich.C5 <- enricher(gene = signif.entrez, 
                      TERM2GENE = C5.entrez, 
                      pvalueCutoff = 0.05,
                      minGSSize = 50,
                      maxGSSize = 500)
enrich.C6 <- enricher(gene = signif.entrez, 
                      TERM2GENE = C6.entrez,
                      pvalueCutoff = 0.05,
                      minGSSize = 50,
                      maxGSSize = 500)


#format results
enrich.C5.df <- enrich.C5@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)

enrich.C6.df <- enrich.C6@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)



#### Visualize results ####
enrich.C5.df %>% 
  filter(p.adjust <= 10e-10) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("C5_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Mtb significant genes (FDR < 10E-16)\nEnriched in C5 gene sets (FDR < 0.01)")


enrich.C6.df %>% 
  filter(p.adjust <= 10e-3) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("C5_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Mtb significant genes (FDR < 10E-16)\nEnriched in C6 gene sets (FDR < 0.01)")


dotplot(enrich.C5, showCategory=20) + ggtitle("dotplot for GSEA")














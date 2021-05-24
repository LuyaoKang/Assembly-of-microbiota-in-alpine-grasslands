# Load libraries

setwd('E:/147asvs/')
library(phyloseq)
library(ape)
library(Biostrings)
library(reshape)
library(ggplot2)


# Construct phyloseq objects
#read in metadata
metadata <- read.delim("./data/metadata.txt", sep = '\t', row.names=1, header = T)
metadata <- sample_data(metadata)

## bacteria
#read in otu table
otu.table <- read.delim('./data/bacteria/otutab_rare.txt',sep="\t", row.names=1, header = T)
otu.table <- as.matrix(otu.table)

#read in taxonomy, seperated by kingdom phylum class order family genus species 
taxonomy <- read.delim('./data/bacteria/taxonomy.txt',sep="\t",row.names=1, header = T)
taxonomy <- as.matrix(taxonomy)

# read in tree
tree <- read_tree('./data/bacteria/otus.nwk')

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./data/bacteria/otus.fa",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)

#merge into one phyloseq object
bac.phylo <- phyloseq(otu.table, taxonomy, metadata, tree, ref_seqs)
bac.phylo
bac.phylo.rel <- microbiome::transform(bac.phylo, "compositional")

## fungi
#read in otu table
otu.table <- read.delim('./data/fungi/otutab_rare.txt',sep="\t", row.names=1, header = T)
otu.table <- as.matrix(otu.table)

#read in taxonomy, seperated by kingdom phylum class order family genus species 
taxonomy <- read.delim('./data/fungi/taxonomy.txt',sep="\t",row.names=1, header = T)
taxonomy <- as.matrix(taxonomy)

# read in tree
tree <- read_tree('./data/fungi/otus.nwk')

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./data/fungi/otus.fa",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)

#merge into one phyloseq object
fungi.phylo <- phyloseq(otu.table, taxonomy, metadata, tree, ref_seqs)
fungi.phylo
fungi.phylo.rel <- microbiome::transform(fungi.phylo, "compositional")

## protist
#read in otu table
otu.table <- read.delim('./data/protist/otutab_rare.txt',sep="\t", row.names=1, header = T)
otu.table <- as.matrix(otu.table)

#read in taxonomy, seperated by kingdom phylum class order family genus species 
taxonomy <- read.delim('./data/protist/taxonomy.txt',sep="\t",row.names=1, header = T)
taxonomy <- as.matrix(taxonomy)

# read in tree
tree <- read_tree('./data/protist/otus.nwk')

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./data/protist/otus.fa",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)

#merge into one phyloseq object
protist.phylo <- phyloseq(otu.table, taxonomy, metadata, tree, ref_seqs)
protist.phylo
protist.phylo.rel <- microbiome::transform(protist.phylo, "compositional")

## animal
#read in otu table
otu.table <- read.delim('./data/animal/otutab_rare.txt',sep="\t", row.names=1, header = T)
otu.table <- as.matrix(otu.table)

#read in taxonomy, seperated by kingdom phylum class order family genus species 
taxonomy <- read.delim('./data/animal/taxonomy.txt',sep="\t",row.names=1, header = T)
taxonomy <- as.matrix(taxonomy)

# read in tree
tree <- read_tree('./data/animal/otus.nwk')

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./data/animal/otus.fa",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)

#merge into one phyloseq object
animal.phylo <- phyloseq(otu.table, taxonomy, metadata, tree, ref_seqs)
animal.phylo
animal.phylo.rel <- microbiome::transform(animal.phylo, "compositional")


# Comunity composition 
# phylum.boxplot, Extract list of top N Taxa, taxrank indicate 1 to kingdom, 2(phylum), 3(class)
new.df <-function(phylo.rel, taxrank, k, N) {
  otu_tab <- otu_table((tax_glom(phylo.rel, taxrank)))
  tax.names <- as.vector(tax_table(tax_glom(phylo.rel,taxrank))[,k])
  rownames(otu_tab) <- tax.names
  rowmean <-sapply(1:nrow(otu_tab),function(x) mean(otu_tab[x,]))
  otu_tab<-otu_tab[order(rowmean,decreasing=TRUE), ]
  taxa_list<-rownames(otu_tab)[1:N]
  new_df<-otu_tab[rownames(otu_tab) %in% taxa_list,]
  new_df <- melt(new_df, id.vars = 'taxa')
  colnames(new_df) <- c('taxa','site','relative_abundance')
  new_df$taxa <- factor(new_df$taxa, ordered = T, levels = rev(taxa_list))
  return(new_df)
}

phylum.plot.fun <- function(new_df) {
  p <- ggplot(new_df, aes(x = taxa, y = relative_abundance *100, fill = taxa)) +
    geom_boxplot(outlier.shape = NA) +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="steelblue", size=1.4, alpha=0.2) +
    theme_classic() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    xlab("")  + ylab('Relative abundance (%)')+
    coord_flip()
  return(p)  
}
new_df_bac <- new.df(bac.phylo.rel, 'Phylum', 2, 10)
new_df_fungi <- new.df(fungi.phylo.rel, 'Phylum', 2, 10)
new_df_protist <- new.df(protist.phylo.rel, 'Phylum', 2, 10)
new_df_animal <- new.df(animal.phylo.rel, 'Class', 3, 10)

bac.phy.boxplot <- phylum.plot.fun(new_df_bac)
fungi.phy.boxplot <- phylum.plot.fun(new_df_fungi)
protist.phy.boxplot <- phylum.plot.fun(new_df_protist)
animal.phy.boxplot <- phylum.plot.fun(new_df_animal)

##determine the class compositions within top 10 phylums##
arrange.tab <- function(phylo, N, taxrank, vect) {
  subphylo <- tax_glom(phylo, taxrank)
  subphylo.rel <- microbiome::transform(subphylo, "compositional")
  ra.tab <- otu_table(subphylo.rel)
  MRA <- rowMeans(ra.tab)
  group <- tax_table(subphylo.rel)[,vect]
  mra.tab <- data.frame(group,MRA)
  colnames(mra.tab) <- c('level1', 'level2', 'MRA')
  #arrange the class table
  library(tidyr)
  mra.tab <- mra.tab %>% spread(level2, MRA)
  mra.tab[is.na(mra.tab)] <- 0
  rownames(mra.tab)<-mra.tab$'level1'
  mra.tab<-as.matrix(t(mra.tab[,-1])*100)
  colsum <-apply(mra.tab,2,sum)
  rowsum<-apply(mra.tab,1,sum)
  top_N_tab<-(mra.tab[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])[,1:N]
  head(top_N_tab)
  top_N_tab<-as.matrix(top_N_tab)
  return(top_N_tab)
}

top10phylum_bcateria <- arrange.tab(bac.phylo, 10, 'Class', c(2,3))
top10phylum_bcateria_tab <-data.frame(phylums =colnames(top10phylum_bcateria),
                                      rela_abun = colSums(top10phylum_bcateria))
top10phylum_bcateria_tab

top10phylum_fungi <- arrange.tab(fungi.phylo, 10, 'Class', c(2,3))
top10phylum_fungi_tab <-data.frame(phylums =colnames(top10phylum_fungi),
                                   rela_abun = colSums(top10phylum_fungi))
top10phylum_fungi_tab

top10phylum_protist <- arrange.tab(protist.phylo, 10, 'Class', c(2,3))
top10phylum_protist_tab <-data.frame(phylums =colnames(top10phylum_protist),
                                     rela_abun = colSums(top10phylum_protist))
top10phylum_protist_tab

#soil animal only has 7 classes
top7class_animal <- arrange.tab(animal.phylo, 7, 'Order', c(3,4))
top7class_animal_tab <-data.frame(class =colnames(top7class_animal),
                                  rela_abun = colSums(top7class_animal))
top7class_animal_tab

# Get the stacked barplot
# create color palette:
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,
          558,43,652,31,610,477,588,99,81,503,562,76,96,495,77,12,90,
          345,255,401,366,276,158,436)

layout(matrix(c(1:8),4,2,byrow = F), c(1,1))
#bacteria
mycol <-colors()[rep(mycol,nrow(top10phylum_bcateria))]
par(mar=c(2,5,2,2))
barplot(top10phylum_bcateria,width=1.8,space=0.4,plot=T,las=2,
        col=mycol[1:nrow(top10phylum_bcateria)],cex.axis=0.8,cex.names=0.7,border=NA,
        xlab = 'Phylum',ylab="Relative abundance(%)",
        offset=0,cex.lab=1)
par(mar=c(2,3.5,3.5,1))
plot.new()
legend("topleft",legend=rownames(top10phylum_bcateria),
       ncol=4,fill=mycol[1:nrow(top10phylum_bcateria)],cex=0.6,bty="n")

#protist
par(mar=c(2,5,2,2))
barplot(top10phylum_protist,width=1.8,space=0.4,plot=T,las=2,
        col=mycol[1:nrow(top10phylum_protist)],cex.axis=0.8,cex.names=0.7,border=NA,
        xlab = 'Phylum',ylab="Relative abundance(%)",
        offset=0,cex.lab=1)
par(mar=c(2,3.5,3.5,1))
plot.new()
legend("topleft",legend=rownames(top10phylum_protist),
       ncol=4,fill=mycol[1:nrow(top10phylum_protist)],cex=0.6,bty="n")

#fungi
par(mar=c(2,5,2,2))
barplot(top10phylum_fungi,width=1.8,space=0.4,plot=T,las=2,
        col=mycol[1:nrow(top10phylum_fungi)],cex.axis=0.8,cex.names=0.7,border=NA,
        xlab = 'Phylum',ylab="Relative abundance(%)",
        offset=0,cex.lab=1)
par(mar=c(2,3.5,3.5,1))
plot.new()
legend("topleft",legend=rownames(top10phylum_fungi),
       ncol=4,fill=mycol[1:nrow(top10phylum_fungi)],cex=0.6,bty="n")

#animal
par(mar=c(2,5,2,2))
barplot(top7class_animal,width=1.8,space=0.4,plot=T,las=2,
        col=mycol[1:nrow(top7class_animal)],cex.axis=0.8,cex.names=0.7,border=NA,
        xlab = 'Phylum',ylab="Relative abundance(%)",
        offset=0,cex.lab=1)
par(mar=c(2,3.5,3.5,1))
plot.new()
legend("topleft",legend=rownames(top7class_animal),
       ncol=4,fill=mycol[1:nrow(top7class_animal)],cex=0.6,bty="n")
par(mfrow=c(1,1))
```
# α diversity
```{r}
alpha_bac <- read.delim('./data/bacteria/result/alpha/vegan.txt',sep="\t", 
                        row.names=1, header = T)
alpha_fungi <- read.delim('./data/fungi/result/alpha/vegan.txt',sep="\t", 
                          row.names=1, header = T)
alpha_protist <- read.delim('./data/protist/result/alpha/vegan.txt',sep="\t", 
                            row.names=1, header = T)
alpha_animal <- read.delim('./data/animal/result/alpha/vegan.txt',sep="\t", 
                           row.names=1, header = T)

alpha_bac <- cbind(tax = rep('Bacteria', 147), alpha_bac[,c(1,4,5)])
alpha_fungi <- cbind(tax = rep('Fungi', 147), alpha_fungi[,c(1,4,5)])
alpha_protist <- cbind(tax = rep('Protist', 147), alpha_protist[,c(1,4,5)])
alpha_animal <- cbind(tax = rep('Animal', 147), alpha_animal[,c(1,4,5)])

diversity <- rbind(alpha_bac, alpha_fungi, alpha_protist, alpha_animal)
melted <- melt(diversity, id.vars = c('tax'))
library(plyr)
ddply(melted, c("tax", "variable"), summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)))

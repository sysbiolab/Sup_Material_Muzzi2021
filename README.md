Script for data analysis
================
João C. D. Muzzi, Mauro A. A. Castro <br>

19 February 2021.


Context
----
Despite progress in understanding the biology of adrenocortical carcinoma (ACC), the treatment options have not changed in the last 3 decades. Our goal was to improve the knowledge on immune pathways that could be druggable targets to enhance the immune response of patients. Our strategy was to revisit the differences based on a molecular classification on High and Low Steroid Phenotypes (HSP and LSP respectively) proposed by Zheng et al. (2016), using The Cancer Genome Atlas (TCGA) ACC database and public data sets published by reference studies of Thorsson et al., (2018) and Sturm et al. (2019) for TCGA samples. This script reproduces all results published in Muzzi et al. (2021) and serves as complementary material.
In Preprocessing.Rmd file, we show how to obtain and preprocess the data sets used in the study, at the end of the section the data sets are saved as RData and can be imported for further analysis. Here we use the preprocessed data and show the steps for data analysis and how to obtain the results and plots shown in Muzzi et al. (2021). 

Package Installation and Loading
----
```r
# Packages Installation
## Find out if packages are already installed, if not proceed with installation
### R packages in CRAN repository
Rpackages <- c("tidyverse","scales","igraph", "ggpubr", "ggforce", "ggrepel",
               "cowplot", "survminer", "survival", "gridExtra")
Rpackages.new <- Rpackages[!(Rpackages %in% installed.packages()[,"Package"])]
if(length(Rpackages.new)) install.packages(Rpackages.new)
### R packages in BioConductor
Rpackages <- c("DESeq2","SummarizedExperiment","tmod", "qvalue","ComplexHeatmap",
               "RedeR")
Rpackages.new <- Rpackages[!(Rpackages %in% installed.packages()[,"Package"])]
if(length(Rpackages.new)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(Rpackages.new)
}
### R packages from other sources
Rpackages <- c("devtools", "ggunchained")
Rpackages.new <- Rpackages[!(Rpackages %in% installed.packages()[,"Package"])]
if(length(Rpackages.new)){
  install.packages("devtools")
  remotes::install_github("JanCoUnchained/ggunchained")
}
```

```r 
# Libraries used
## For Differential Expression (D.E.) analysis
library(SummarizedExperiment)
library(DESeq2)
## For enrichment analysis (CERNO test)
library(tmod)
## Data manipulation and statistics
library(tidyverse)
library(qvalue)
library(scales)
library(ggpubr)
## Plots
library(ggforce)
library(ggrepel)
library(ggunchained)
library(cowplot)
### Table construction
library(gridExtra)
### Heatmap construction
library(ComplexHeatmap)
## Survival analysis
library(survminer)
library(survival)
## Network reconstruction
library(RedeR)
library(igraph)
```

Import Preprocessed data
----

```r
download.file(url = "https://github.com/sysbiolab/Sup_Material_Muzzi2021/blob/main/Pancan_data.RData", destfile = "Pancan_data.RData") 
download.file(url = "https://github.com/sysbiolab/Sup_Material_Muzzi2021/blob/main/tcgaACC_pre_processed.RData", destfile = "tcgaACC_pre_processed.RData")

### Importing pre-processsed data
load("tcgaACC_pre_processed.RData")
load("Pancan_data.RData")
## Retrieving Thorsson et al 2018 preprocessed data set
thorsson_pan <- Pancan_data[["thorsson_pan"]]
## Retrieving Sturm et al 2019 preprocessed data set
sturm_pan <- Pancan_data[["sturm_pan"]]

## Retrieving column annotations from preprocessed SummarizedExperiment object
## This annotation is composed by Zheng et al (2016) original clinical information for 
## this cohort, Thorsson et al. (2018) immune features data set and Sturm et al. (2019)
## xCELL scores for immune infiltrate in TCGA-ACC patients.
clinic.dat <- data.frame(colData(tcgaACC))

```

Reproduction of Figure 1 at Muzzi et al. (2021)
----

```r 

# Boxplot - pan-cancer comparison for Leukocyte Fraction - ACC cohort
dat_imun <- 
  data.frame(cbind(barcode=thorsson_pan$TCGA.Participant.Barcode,
                   TCGA.Study=thorsson_pan$TCGA.Study,
                   Leukocyte.Fraction=thorsson_pan$Leukocyte.Fraction))
dat_imun$Leukocyte.Fraction <- as.numeric(dat_imun$Leukocyte.Fraction)
## Order cohorts by median
t <- aggregate(dat_imun[,"Leukocyte.Fraction"],
               list(dat_imun$TCGA.Study),
               median)
t <- t[order(t$x),]
dat_imun$TCGA.Study <- factor(dat_imun$TCGA.Study,
                              levels=t$Group.1)

## Boxplot
g<- ggplot(dat_imun,
           aes(x=TCGA.Study, y= Leukocyte.Fraction))+
  stat_boxplot(geom='errorbar', lwd=0.2)+
  geom_boxplot(fill="white", outlier.size = 0.1, lwd=0.1)+
  stat_boxplot(data=dat_imun[dat_imun$TCGA.Study=="ACC",],
               aes(x=TCGA.Study, y= Leukocyte.Fraction),
               geom='errorbar', 
               color="green3")+
  geom_boxplot(data=dat_imun[dat_imun$TCGA.Study=="ACC",],
               aes(x=TCGA.Study, y= Leukocyte.Fraction),
               fill="green3", 
               outlier.size = 0.1,
               outlier.colour =  "green3",
               lwd=0.1)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill="grey94"))+
  ylab(label="Leukocyte Fraction")+
  xlab(label="TCGA Cohorts")

# Boxplot - pan-cancer comparison for Leukocyte Fraction - ACC distinction between Steroid Phenotypes
dat_imun$TCGA.Study <- 
  ifelse(
    dat_imun$barcode %in% clinic.dat$patient[clinic.dat$Steroid=="Steroid_High"],
    "ACC HSP",
    ifelse(
      dat_imun$barcode %in% clinic.dat$patient[clinic.dat$Steroid=="Steroid_Low"],
      "ACC LSP",
      as.character(dat_imun$TCGA.Study)))
## Order cohorts by median
t <- aggregate(dat_imun[,"Leukocyte.Fraction"], 
               list(dat_imun$TCGA.Study),
               median)
t<- t[order(t$x),]
dat_imun$TCGA.Study <- factor(dat_imun$TCGA.Study, 
                              levels=t$Group.1)
dat_imun <- dat_imun[dat_imun$TCGA.Study!="ACC",] #only 1 ACC patient without Steroid notation

## Boxplot
g1<- ggplot(dat_imun,
            aes(x=TCGA.Study, y=Leukocyte.Fraction))+
  stat_boxplot(geom='errorbar', lwd=0.2)+
  geom_boxplot(fill="white", outlier.size=0.1, lwd=0.1)+
  stat_boxplot(data= dat_imun[dat_imun$TCGA.Study=="ACC HSP" |
                              dat_imun$TCGA.Study=="ACC LSP",],
               aes(x=TCGA.Study, y= Leukocyte.Fraction), 
               geom='errorbar', 
               color= c("#F8766D","#00BFC4"))+
  geom_boxplot(data= dat_imun[dat_imun$TCGA.Study=="ACC HSP" |
                              dat_imun$TCGA.Study=="ACC LSP",],
               aes(x=TCGA.Study, y= Leukocyte.Fraction),
               fill=c("#F8766D","#00BFC4"),
               outlier.size = 0.1, 
               outlier.color = c("#F8766D"),
               lwd=0.1)+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill="grey94"))+
  ylab(label="Leukocyte Fraction")+
  xlab(label="TCGA Cohorts")
plot_grid(g,g1,labels="AUTO", ncol=1)
```

Reproduction of Figure 2 at Muzzi et al. (2021)
----
```r
# Boxplot - pan-cancer comparison for Immune Infiltrate
## Merging Sturm and Thorsson datasets by barcode
dat <- subset(sturm_pan, 
              sturm_pan$patient_XCELL %in% thorsson_pan$TCGA.Participant.Barcode)
dat_imun <- subset(thorsson_pan,
                   thorsson_pan$TCGA.Participant.Barcode %in% dat$patient_XCELL)
dat_imun<- dat_imun[order(dat_imun$TCGA.Participant.Barcode),]

## Selecting CD8 T cells and NK cells for xCell and Cibersort
dat_imun <- 
  cbind(
    dat_imun[, c("TCGA.Participant.Barcode", "TCGA.Study")], 
    sturm_pan[order(dat$patient_XCELL), c("T.cell.CD8._XCELL", "T.cell.NK_XCELL")],
    dat_imun[, c("T.Cells.CD8", "NK.Cells.Activated")]*dat_imun$Leukocyte.Fraction)
## Cibersort score must be multiplied by Leukocyte Fraction for comparison inter samples

dat_imun$TCGA.Study <- ifelse(
  dat_imun$TCGA.Participant.Barcode %in% clinic.dat$patient[clinic.dat$Steroid=="Steroid_High"],
  "ACC HSP",
  ifelse(
      dat_imun$TCGA.Participant.Barcode %in% clinic.dat$patient[clinic.dat$Steroid=="Steroid_Low"],
      "ACC LSP",
      as.character(dat_imun$TCGA.Study)))

names(dat_imun) <- c( "TCGA.Participant.Barcode","TCGA.Study",
                      "T.cell.CD8\nXCELL", "T.cell.NK\nXCELL",
                      "T.Cells.CD8\nCIBERSORT",
                      "NK.Cells.Activ.\nCIBERSORT")
## Rescale to 0 to 1 range
dat_imun[,3:6] <- lapply(dat_imun[,3:6], rescale)
## Defining maximum y values for each panel
n<- c(0.4,1,0.5,0.4)

## Boxplot
for(i in 3:6){
  dat <- dat_imun
  ### Order cohorts by median
  t <- aggregate(dat[,i], 
                 list(dat$TCGA.Study), 
                 function(a){median(a,na.rm=T)})
  t<- t[order(t$x),]
  dat$TCGA.Study <- factor(dat$TCGA.Study,
                           levels=t$Group.1)
  dat <- dat[dat$TCGA.Study!="ACC",]
  
  ### Names for y axis
  cols <- names(dat)[i]
  ### Plots
  g <- ggplot(dat, aes(x=TCGA.Study, y= dat[,i]))+
    stat_boxplot(geom='errorbar', lwd=0.2)+
    geom_boxplot(fill="white",
                 outlier.shape = NA,
                 lwd=0.1)+
    stat_boxplot(data= dat[dat$TCGA.Study=="ACC HSP" |
                           dat$TCGA.Study=="ACC LSP",],
                 aes(x=TCGA.Study,
                     y= dat[dat$TCGA.Study=="ACC HSP" |
                            dat$TCGA.Study=="ACC LSP",i]),
                 geom='errorbar', 
                 color= c("#F8766D","#00BFC4"))+
    geom_boxplot(data= dat[dat$TCGA.Study=="ACC HSP" |
                           dat$TCGA.Study=="ACC LSP",],
                 aes(x=TCGA.Study, 
                     y= dat[dat$TCGA.Study=="ACC HSP" |
                            dat$TCGA.Study=="ACC LSP",i]),
                 fill=c("#F8766D","#00BFC4"),
                 outlier.shape = NA,
                 lwd=0.1)+
    coord_cartesian(ylim=c(0,n[i-2]))+
    theme(
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 6), 
          axis.text.y = element_text(size=6),
          panel.background = element_rect(fill="grey94"), 
          axis.title = element_text(size=7))+
    ylab(label= paste0(cols," score"))+
    xlab(label= "TCGA Cohorts")
  if(i==3){g1<-plot_grid(g)}
  if(i==4){g2<-plot_grid(g)}
  if(i==5){g3<-plot_grid(g)}
  if(i==6){g4<-plot_grid(g)} 
}

### Plotting panels
 plot_grid(g1,g2,g3,g4, ncol = 1, labels = "AUTO", label_size = 11)
```

Reproduction of Figure 3 at Muzzi et al. (2021)
----
```r

######## Differential expression analysis with DESeq2
dds <- DESeqDataSet(tcgaACC, design = ~ Steroid)
dds <- DESeq(dds)

# Changing rownames for symbols to plot
## Finding duplicated Symbols
(dds@rowRanges@elementMetadata@listData[["external_gene_name"]]
  [duplicated(dds@rowRanges@elementMetadata@listData[["external_gene_name"]])])
#[1] "PINX1" "TMSB15B" "ATF7" "MATR3" 
## Changing duplicated Symbols, adding _ after duplicated gene name
dds@rowRanges@elementMetadata@listData[["external_gene_name"]][
  match(c("PINX1", "TMSB15B", "ATF7", "MATR3"),
         dds@rowRanges@elementMetadata@listData[["external_gene_name"]])] <-
  c("PINX1_", "TMSB15B_", "ATF7_", "MATR3_")
anyDuplicated(dds@rowRanges@elementMetadata@listData[["external_gene_name"]]) #0

# Changing DDS gene names from Ensembl to Symbol
rownames(dds)<-dds@rowRanges@elementMetadata@listData[["external_gene_name"]]

# Changing some genes to friendly names
rownames(dds)[match(c("CD274","PDCD1","PDCD1LG2","AC011462.1",
                      "FLT1", "TNFRSF14", "TNFSF4", "CD276", "TNFSF9"), 
                    rownames(dds))] <- c("PDL1","PD1","PDL2","TGFB1",
                                         "VEGFR1", "HVEM", "OX40L", "B7-H3", "4.1BBL")
## In 2021 new version of TCGABiolinks TGFB1 is named "AC011462.1" - ensembl ENSG00000105329.8

# Presenting D.E. results
res <- results(dds)

## Generating data.frame to export results with q-value and FDR
res.df <- data.frame(res) #19174 genes
res.df <- res.df[res.df$baseMean>0,] #19171 genes
res.df <- res.df[order(res.df$padj),]
### Q-value and FDR calculus
qobj <- qvalue(p= res.df$pvalue)
res.df$qvalue <- qobj$qvalues
res.df$FDR <- qobj$lfdr
#hist(qobj) #p-value density histogram
res.df <- cbind(GENE_SYMBOL=rownames(res.df),
                res.df) 
res.df <- res.df[order(res.df$padj),]
#Sup table 1
#writexl::write_xlsx(res.df, 
#                    path = "/home/muzzi/ACC/Sup_table_1_Dif_exp_genes.xls") 

res.df$'D.E.genes' <- 
  ifelse(res.df$padj<(0.05) & res.df$log2FoldChange >0.5,
         "up",
          ifelse(res.df$padj<(0.05) & res.df$log2FoldChange <(-0.5), 
                 "down", 
                 "NS"))

###### Normalizing counts - vst = variance stabilizing transformation from DESeq2
vsd <- vst(dds, blind = F) # for information on "blind" argument, see DESeq2 vignette

# Extracting normalized gene expression matrix
gen_expt <- as.data.frame(t(assay(vsd)))
# Volcano plot for D.E. results
g<-ggplot(data = res.df[complete.cases(res.df$padj),],
          aes(x=log2FoldChange, y=-log10(padj), color=D.E.genes),
          alpha=0.6)+
  geom_point(size=0.5)+
  scale_color_manual(values=c("blue","grey","red"))+
  geom_hline(yintercept=(-log10(0.05)), linetype=2, color="#800026", size=1)+
  geom_vline(xintercept=0.5, linetype=2, color="#800026")+
  geom_vline(xintercept=-0.5, linetype=2, color="#800026")+
  geom_label_repel(
         data=res.df[1:20,], 
         aes(label=GENE_SYMBOL),
         show.legend = F,
         size=2)+ 
  annotate(
         geom="text",
         y=3, x=-13,
         label="p adj < 0.05", 
         color="#800026", 
         size=3, 
         hjust = 0, 
         fontface="bold")+
  annotate(
         geom="text",
         y=50, x=-7.5,
         label="Down-regulated in LSP\n Up-regulated in HSP", 
         color="blue", 
         size=2.5,
         hjust = "center", 
         fontface="bold")+
  annotate(
         geom="text",
         y=50, x=7.5,
         label="Up-regulated in LSP\nDown-regulated in HSP", 
         color="red",
         size=2.5,
         hjust = "center",
         fontface="bold")+
  annotate(
         geom="text",
         y=-2, x=-1,
         label="LFC < -0.5", 
         color="blue",
         size=3, 
         hjust = "right", 
         fontface="bold")+
  annotate(
         geom="text",
         y=-2, x=1,
         label="LFC > 0.5", 
         color="red",
         size=3,
         hjust = "left", 
         fontface="bold")+
  labs(title="LSP versus HSP")+
  ylab("-Log10 (P-adj)") + 
  xlab("Log2 Fold Change") + 
  coord_fixed(
        ratio = 0.5,
        xlim = c(-12,12),
        ylim = c(-1,50), 
        expand = TRUE, clip = "on")+
  theme(
        plot.title=element_text(hjust=0.5),
        legend.position="bottom",
        legend.margin = margin(-5, -5, 5, -5),
        legend.spacing.y = unit(0, 'cm'), 
        title = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.text.y=element_text(size=7), 
        axis.ticks = element_line(linetype=1, color="grey"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=8 ),
        legend.key.size = unit(0,"mm"),
        legend.background = element_rect(fill=NULL, color=NULL),
        axis.line = element_blank(),
        panel.grid.major.x = element_line(linetype=111, color="grey80", size = 0.4),
        panel.grid.major = element_line(linetype=3, color="grey", size = 0.2),
        panel.background = element_rect(fill = "grey98", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1))+
  guides(color = guide_legend(label.position = "bottom",
                              title.position="bottom",
                              title.hjust=0.5,
                              override.aes=list(size=3,
                                                shape=21,
                                                fill=c("blue","grey","red"))))
g
```

Reproduction of Figure 4 at Muzzi et al. (2021)
----

```r 

######################################################### Beginning - Prof Mauro Castro script
p.threshold <- function(pvals, alpha=0.05, method="BH") {
  pvals <- sort(pvals)
  padj <- p.adjust(pvals, method = method)
  thset <- which(padj <= alpha)
  if(length(thset)>0){
    mx1 <- mx2 <- which.max(thset)
    if(mx2<length(padj)) mx2 <- mx2 + 1
    th <- (pvals[mx1] + min(pvals[mx2],alpha) ) / 2
  } else {
    th <- min(c(alpha,pvals))
  }
  return(th)
}
######################################################### Adapted from Prof Mauro Castro script
my_dot_plot_CERNO <- function(gsea_obj, top=20, abrv=60, fdr_level=0.05, 
                              title="", ylab="Hallmarks (gene set size)") {
  
  # get fdr threshold
  fdr_threshold <- p.threshold(gsea_obj$pval, alpha = fdr_level, method = "BH")
  
  # adjust pval
  gsea_obj$padj <- p.adjust(gsea_obj$pval, method = "BH")
  
  # Subset results
  df <- gsea_obj %>%
    purrr::when(!is.null(top) ~ head(., top), ~ .)
  df <- as.data.frame(df)
  
  # Handle empty dataframes
  if (nrow(df) == 0) return(ggempty())
  
  # Order by significance value
  df <- df[order(df$pval, decreasing = F),]
  
  # Abbreviate labels
  label.abrv <- substr(df$pathway, 1, abrv)
  df$pathway <- label.abrv
  
  at<- 10^(-c(0,10,20,30,40))
  if(all(at>fdr_level))at <- sort(c(max(at),fdr_level))
  
  df$yla <- nrow(df):1*6.5
  df$xla <- -log10(df$padj)
  df$down <- ifelse(df$padj>fdr_threshold,0,df$down)
  df$up <- ifelse(df$padj>fdr_threshold,0,df$up)  
  df$NS <- ifelse(df$padj>fdr_threshold,df$size,df$NS)
  
  df<- gather(df, key=pie, value = cols, "up","NS","down")
  
  g<- ggplot(df,aes(x=xla,y=yla))+
      geom_arc_bar(data=df,
                   mapping = aes(x0=xla, y0=yla, r0=0, r=3, amount=cols, fill=pie),
                   stat="pie", 
                   inherit.aes = F, 
                   alpha=0.7, 
                   col=NA, 
                   show.legend = c(fill=T))+
      coord_fixed(xlim = c(0,40))+
      scale_fill_manual(values = c("blue","gray","red"), aesthetics="fill")+
      geom_vline(xintercept=(-log10(fdr_threshold)), linetype=2, color="#800026") +
      labs(title=title, fill="D.E. genes")+
      scale_x_continuous(breaks = -log10(at),
                         expand = expansion(mult = c(0.07,0.1)), 
                         labels = format(-log10(at), digits = 1)) +
      scale_y_continuous( position = "right", breaks= df$yla,labels= df$pathway)+
      ylab(ylab) + 
      xlab("-Log10(P-value)") + 
      annotate(geom="text",
               y=1.5,
               x=-log10(fdr_threshold-(fdr_threshold/10))+1,
               label=paste0("FDR<",fdr_level), 
               color="#800026", 
               size=3, 
               hjust = 0,
               fontface="bold") +
      theme(
        plot.title=element_text(hjust=0),
        legend.position="bottom",
        legend.margin = margin(-5, -5, 5, -5),
        legend.spacing.y = unit(0, 'cm'), 
        title = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.text.y=element_text(size=7), 
        axis.text.x = element_text(size=10),
        axis.ticks = element_line(linetype=1, color="grey"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=8 ),
        axis.line = element_blank(),
        panel.grid.major.x = element_line(linetype=111, color="grey80", size = 0.4),
        panel.grid.major = element_line(linetype=3, color="grey", size = 0.2),
        panel.background = element_rect(fill = "grey98", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1)) + 
      guides(fill = guide_legend(label.position = "bottom",
                                 title.position="bottom",
                                 title.hjust=0.5,
                                 override.aes=list(size=5,
                                                   shape=21,
                                                   colour="white")))
  return(g)
}
######################################################### End - Prof Mauro Castro script
```


```r 
######## Enrichment analysis with tmod
# Genes list ordered by p-value
l <- rownames(res)[order(res$padj,decreasing = F)]

##### MSig modules to tmod - Hallmarks
# Download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/msigdb_v7.1.xml
download.file(
  url="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/msigdb_v7.1.xml",
  destfile = "msigdb_v7.1.xml")
msig <- tmodImportMSigDB("msigdb_v7.1.xml")
# H = Hallmarks gene sets
sel <- msig$MODULES$Category == "H"
resCERNO <- tmodCERNOtest(l, 
                          mset = msig[sel], 
                          qval = 1)
# Preparing data.frame with results, q-value and FDR to export
# q value and FDR calculus
res.df <- data.frame(resCERNO)
qobj <- qvalue(p= res.df$P.Value, pi0 = 1)
# pi0=1 - Uses Benjamini–Hochberg estimation.
# Little p-values quantity doesn't allow pi estimation.
res.df$qvalue <- qobj$qvalues
res.df$FDR <- qobj$lfdr
res.df.H <- res.df
# Defines for each modules how many genes are up, down or not significant based on D.E. results
# adj pvalue Cutoff = 0.05
# |log2 fold change| Cutoff = 0.5
pie <- tmodDecideTests(rownames(res),
                       lfc = res$log2FoldChange, 
                       pval = res$padj,
                       mset = msig[sel], #Hallmarks
                       labels = "LSP vs HSP")
res.df.H <- cbind(res.df.H,
                  pie[["LSP vs HSP"]][rownames(res.df.H),])
pie[["LSP vs HSP"]] <- subset(pie[["LSP vs HSP"]], 
                              rownames(pie[["LSP vs HSP"]]) %in% rownames(res.df))

# Creating data.frame for input in function to generate plot
df <- cbind(resCERNO,
            pie[["LSP vs HSP"]][rownames(resCERNO),])

df$Title <- gsub("Hallmark ",
                 "",
                 df$Title)
df$Title <- paste0("(",df$N1,") ",toupper(df$Title))
df<- data.frame(cbind(pathway=paste0(" ",df$Title),
                      size=df$N1, 
                      pval=as.numeric(df$P.Value),
                      down=as.numeric(df$down),
                      NS=df$N,
                      up=as.numeric(df$up)))
df[,2:ncol(df)] <- lapply(df[,2:ncol(df)],
                          as.numeric)

# Plot results
g1 <- my_dot_plot_CERNO(df, top=20, title = "     LSP versus HSP")

##### Using Li et al (2016) BTMs modules
resCERNO <- tmodCERNOtest(l, 
                          mset="LI", 
                          qval = 1)
# Preparing data.frame with results, q-value and FDR to export
# q value and FDR calculus
res.df <- data.frame(resCERNO)
qobj <- qvalue(p= res.df$P.Value, pi0 = 1)
# pi0=1 - Uses Benjamini–Hochberg estimation. 
# Little p-values doesn't allow pi estimation.
res.df$qvalue <- qobj$qvalues
res.df$FDR <- qobj$lfdr
res.df.Li <- res.df

# Defines for each modules how many genes are up, down or not significantly based on D.E. results
# adj pvalue Cutoff = 0.05
# log2 fold change = 0.5
pie <- tmodDecideTests(rownames(res),
                       lfc = res$log2FoldChange, 
                       pval = res$padj, 
                       labels = "LSP vs HSP", 
                       mset = "LI")
res.df.Li <- cbind(res.df.Li,
                   pie[["LSP vs HSP"]][rownames(res.df.Li),])
pie[["LSP vs HSP"]] <- subset(pie[["LSP vs HSP"]], 
                              rownames(pie[["LSP vs HSP"]]) %in% rownames(res.df))
# Creating data.frame for input in function to generate plot
df <- cbind(resCERNO,
            pie[["LSP vs HSP"]][rownames(resCERNO),])
df$Title <- paste0("(",df$N1,") ",toupper(df$Title))
## Two gene sets names too large for panel:
# df$Title[4] #"(78) REGULATION OF ANTIGEN PRESENTATION AND IMMUNE RESPONSE"
# df$Title[13] #"(29) MYELOID CELL ENRICHED RECEPTORS AND TRANSPORTERS"     
df$Title[4] <- "(78) REGULATION OF ANTIGEN PRESENTATION"
df$Title[13] <- "(29) MYELOID CELL ENRICHED RECEPTORS"     
df<- data.frame(cbind(pathway=(df$Title),
                      size=df$N1, 
                      pval=as.numeric(df$P.Value),
                      down=as.numeric(df$down),
                      NS=df$N,
                      up=as.numeric(df$up)))
df[,2:ncol(df)] <- lapply(df[,2:ncol(df)],
                          as.numeric)
# Plot results
g2<- my_dot_plot_CERNO(df, top=20, title="     LSP versus HSP",
                       ylab = "Blood Transcriptome Modules (gene set size)")
#Sup table 2
#writexl::write_xlsx(x= list(CERNO.Hallmarks=res.df.H,
#                            CERNO.LI=res.df.Li),
#                    path= "/home/muzzi/ACC/Sup_table_2_CERNO_steroid.xls") 
g1 # Ref fig 4A Muzzi et al. 2021
g2 # Ref fig 4B Muzzi et al. 2021
```

Reproduction of Figure 5 at Muzzi et al. (2021)
----
```r
# Preparing data.frame for heatmap and boxplot - Cibersort (Thorsson et al) & xCell (Sturm et al)

## Retrieving Cibersort and xCell immune infiltrate scores for ACC patients
dat_imun <- cbind(clinic.dat[,c(51:72)],
                  clinic.dat[,80:118])
## Removing columns with 0 counts - Cibersort T.cells.gamma.delta
dat_imun <- dat_imun[,colSums(dat_imun)>0]
## Adjusting column names
names(dat_imun)[58:60] <- c("Immune.Score","Stroma.Score","Microenvironment.Score")

#### Heatmap xCell x Cibersort
# # Z-score calculus for immune infiltrate scores and setting maximum and minimum values
dat_imun[,1:57] <- scale(dat_imun[,1:57])
dat_imun[dat_imun>2] <- 2
dat_imun[dat_imun<(-2)] <- (-2)

# Removing "XCELL" from column names
names(dat_imun) <- gsub("_XCELL","", names(dat_imun))

#### Annotation Data
annot_hm<- clinic.dat[order(rownames(clinic.dat)),
                      c("Steroid","vital_status","Cortisol", "Other.Hormones",
                         "Immune.Subtype", "Leukocyte.Fraction",
                         "Lymphocyte.Infiltration.Signature.Score","Macrophage.Regulation",
                         "Wound.Healing","IFN.gamma.Response", "TGF.beta.Response")]
names(annot_hm)[c(2,7)] <- c("Vital.Status", "Lymphocyte.Inf.Score")

# Clustering patients between Steroid Phenotypes
## Clustering Low Steroid patients
cols <-hclust(dist(dat_imun[rownames(annot_hm)[annot_hm$Steroid=="Steroid_Low"],],
                   method="euclidean"),
              method = "mcquitty")
cols <- cols$labels[rev(cols$order)]
## Clustering High Steroid patients
cols2 <-hclust(dist(dat_imun[rownames(annot_hm)[!annot_hm$Steroid=="Steroid_Low"],],
                    method="euclidean"),
               method = "mcquitty")
cols2 <- cols2$labels[rev(cols2$order)]
# Reorder data_frames with clustered patients 
annot_hm<- annot_hm[c(cols,cols2),]
dat_imun <- dat_imun[rownames(annot_hm),]

# rescale arbitrary values to 0 to 1 range
annot_hm[,7:11] <- lapply(annot_hm[,7:11],rescale)
dat_imun[, 58:60] <- lapply(dat_imun[,58:60], rescale)

### Preparing p-values annotations for Heatmap
### Method Wilcoxon - non-parametric
annot_row <- 
  data.frame(apply(dat_imun,
                   2,
                   function(a){wilcox.test(formula = a~annot_hm$Steroid)[3]}))
names(annot_row) <- names(dat_imun)
annot_row <- -log10(annot_row)

### Heatmap with ComplexHeatmap package
#### Columns Annotation
ha<- columnAnnotation(
                df=annot_hm[,c(1:4)],
                simple_anno_size=unit(4,"mm"),
                height=unit(1,"cm"),
                na_col="white",
                annotation_name_gp = gpar(fontsize=8),
                annotation_name_side = "right", 
                annotation_legend_param = list(
                  Cortisol=list(
                    labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=8, fontface="bold")),
                  Steroid=list(
                    nrow=2, labels=c("HSP","LSP"),
                    labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=8, fontface="bold")),
                  Vital.Status=list(
                    title="Vital Status", 
                    labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=8, fontface="bold")),
                  Other.Hormones=list(
                    title="Other Hormones",
                    labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=8, fontface="bold"))),
                 col=list(
                    Steroid=c("Steroid_High"="#F8766D", "Steroid_Low"="#00BFC4"),
                    Vital.Status=c("Alive"="azure","Dead"="black"),
                    Other.Hormones=c("Sexual"="pink","Mineralcorticoids"="yellow","No"="grey"),
                    Cortisol=c("Cortisol"="navyblue","No"="grey")))
#### Main heatmap for Cibersort results
h1<- Heatmap(matrix = t(dat_imun[,1:21]),
             name= "Cibersort Immune infiltrate scores\nThorsson et al., 2018",
             show_heatmap_legend = T,
             show_column_names = F,
             column_title_gp = gpar(fontsize=10),
             column_title = c("LSP (n=31)","HSP (n=47)"),
             cluster_columns = F, 
             column_split = factor(annot_hm$Steroid, levels=c("Steroid_Low","Steroid_High")),  
             cluster_rows = F,
             row_title = "Cibersort Immune \ninfiltrate scores", 
             row_title_gp = gpar(fontsize=10),
             row_title_side = "left",
             row_names_side = "right",
             row_names_gp = gpar(fontsize=8),
             heatmap_legend_param = list(
               at=c(-2,0,2),
               legend_width= unit(1.7,"cm"),
               labels_gp=gpar(fontsize=8), 
               title_gp=gpar(fontsize=8, fontface="bold"),
               border=T, 
               direction = "horizontal"), 
               col= circlize::colorRamp2(c(-2,0,2), c("darkseagreen","white","orange1")),
             top_annotation = columnAnnotation(
               df=annot_hm[,5:7],
               annotation_name_side = "right",
               annotation_name_gp = gpar(fontsize=8),
               simple_anno_size=unit(3,"mm"),
               annotation_legend_param = list(
                    Immune.Subtype=list(
                      nrow=3, title="Immune Subtype",
                      labels_gp=gpar(fontsize=8), 
                      title_gp=gpar(fontsize=8, fontface="bold")),
                    Leukocyte.Fraction=list(
                      title="Leukocyte\nFraction",
                      border=T, legend_width= unit(1,"cm"),
                      labels_gp=gpar(fontsize=8),
                      title_gp=gpar(fontsize=8, fontface="bold"),
                      direction="horizontal", at=c(0,0.4)),
                     Lymphocyte.Inf.Score=list(
                       title="Lymphocyte\nInf. Score",
                       border=T, legend_width= unit(1,"cm"),
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       direction="horizontal", at=c(0,1))),
                col=list(
                    Immune.Subtype=c("C1"="red","C2"="yellow","C3"="green3",
                                     "C4"="cyan","C5"="blue","C6"="pink"),
                    Leukocyte.Fraction=circlize::colorRamp2(c(0,0.4), c("white","blue")),
                    Lymphocyte.Inf.Score=circlize::colorRamp2(c(0,1), c("white","darkgreen")))),
              right_annotation = rowAnnotation(
                df=t(annot_row[,1:21]), 
                annotation_legend_param = list(
                  at=c(0,2,10),
                  labels=c("",2,10),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  border=T, 
                  legend_width= unit(1.7,"cm"),
                  direction="horizontal",
                  nrow=2),
                annotation_label = "-log(P.Value)", 
                show_annotation_name = F,
                col=list( V1= circlize::colorRamp2(c(0,2,10), c("white","white", "black"))
                                              )))
#### Main heatmap for xCell results
h2<- Heatmap(matrix = t(dat_imun[,22:57]),
             name= "xCell Immune infiltrate scores\nSturm et al., 2019",
             show_heatmap_legend = T,
             show_column_names = T,
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_split = factor(annot_hm$Steroid, levels=c("Steroid_Low","Steroid_High")),
             cluster_columns = F,
             cluster_rows = F,
             row_title_gp = gpar(fontsize=10),
             row_names_side = "right",
             row_names_gp = gpar(fontsize=8),
             row_title_side = "left", 
             row_title = "xCell Immune \ninfiltrate scores",
             heatmap_legend_param = list(
               at=c(-2,0,2),
               legend_width= unit(1.7,"cm"),
               labels_gp=gpar(fontsize=8),
               title_gp=gpar(fontsize=8, fontface="bold"),
               border=T,
               direction = "horizontal"), 
               col= circlize::colorRamp2(c(-2,0,2), c("darkseagreen","white","orange1")),
             top_annotation = columnAnnotation(
               df=dat_imun[,58:60], 
               simple_anno_size=unit(3,"mm"),
               annotation_name_side = "right",
               annotation_name_gp = gpar(fontsize=8), 
               height = unit(1,"cm"),
               annotation_legend_param = list(
                     Immune.Score=list(
                       title="Immune\nScore",
                       border=T,
                       legend_width= unit(1,"cm"),
                       labels_gp=gpar(fontsize=8), 
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       direction="horizontal",
                       at=c(0,1)),
                     Stroma.Score=list(
                       title="Stroma\nScore",
                       border=T, 
                       legend_width= unit(1,"cm"),
                       labels_gp=gpar(fontsize=8), 
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       direction="horizontal", 
                       at=c(0,1)),
                    Microenvironment.Score=list(
                      title="Microenv.\nScore",
                      border=T, 
                      legend_width= unit(1,"cm"),
                      labels_gp=gpar(fontsize=8),
                      title_gp=gpar(fontsize=8, fontface="bold"),
                      direction="horizontal", at=c(0,1))),
                col=list(
                     Immune.Score= circlize::colorRamp2(c(0,1), c("white","mediumvioletred")),
                     Stroma.Score= circlize::colorRamp2(c(0,1), c("white","darkgoldenrod3")),
                     Microenvironment.Score= circlize::colorRamp2(c(0,1), c("white","sienna")))),
              right_annotation = rowAnnotation(
                df=t(annot_row[,22:57]),
                annotation_name_rot=0,
                annotation_legend_param = list(
                  at=c(0,2,10),
                  labels=c("",2,10),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  border=T, 
                  legend_width= unit(1.7,"cm"),
                  direction="horizontal",
                  nrow=2),
                annotation_label = "-log(P.Value)",
                annotation_name_gp = gpar(fontsize=8),
                col=list( V1= circlize::colorRamp2(c(0,2,10), c("white","white", "black"))
                                              )))
#### Concatenate heatmaps                 
h_list<- ha %v% h1 %v% h2
#### Saving object for panel
g<- grid.grabExpr(
  draw(h_list,
       ht_gap=unit(c(1,3),"mm"),
       cluster_columns=F,
       heatmap_legend_side="bottom",
       annotation_legend_side="bottom",
       merge_legends=T),
  height = 15,
  width = 5)

```


```r
# Boxplot xCell & Cibersort
## Selecting Monocytes, T CD8 cells and NK cells from xCell and Cibersort
dat_bp <- cbind(annot_hm[,c("Steroid","Immune.Subtype")],
                dat_imun[rownames(annot_hm),
                         c("Monocyte","T.cell.CD8.",
                           "T.cell.NK", "Monocytes", 
                           "T.Cells.CD8", "NK.Cells.Activated")])
names(dat_bp) <- c( "Steroid","Immune.Subtype", "Monocyte\nXCELL", "T.cell.CD8\nXCELL",
                    "T.cell.NK\nXCELL", "Monocytes\nCIBERSORT", "T.Cells.CD8\nCIBERSORT",
                    "NK.Cells.Activ.\nCIBERSORT")
## Rescale to 0 to 1 range
dat_bp[,3:8] <- lapply(dat_bp[,3:8], rescale)
## Preparing data.frame for boxplot
dat_bp <- pivot_longer(dat_bp,
                       cols = -contains(c("Steroid","Immune.Subtype")),
                       values_to = "value",
                       names_to = "variable")
dat_bp$variable <- factor(dat_bp$variable, 
                          levels = c("Monocyte\nXCELL", "Monocytes\nCIBERSORT",
                                     "T.cell.CD8\nXCELL","T.Cells.CD8\nCIBERSORT",
                                     "T.cell.NK\nXCELL", "NK.Cells.Activ.\nCIBERSORT"))
dat_bp$Steroid <- factor(dat_bp$Steroid,
                         levels = c("Steroid_Low", "Steroid_High"))
## Plotting
g1<- ggplot(dat_bp,
            aes(x=Steroid, y= value, fill=Steroid))+
  facet_wrap(.~variable, ncol=1, scales = "free_y")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+  
  ggunchained::geom_split_violin(
    aes(fill=NULL, col=Immune.Subtype),
    alpha=0.4, 
    show.legend = T)+
  scale_color_manual(values=c("green3","cyan"), 
                     labels=c("C3 (HSP n=7)\n     (LSP n=16)",
                              "C4 (HSP n=37)\n      (LSP n=12)"))+
  labs(y="Immune Cells Infiltrate score", x= NULL) +
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=8),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(5,"mm"),
        legend.text=element_text(size=8), 
        legend.direction = "vertical",
        legend.box = "vertical" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 1.1,
                     label="p.format",
                     size=3)+
  scale_x_discrete(breaks=c("Steroid_Low", "Steroid_High"),
                   labels=c("LSP", "HSP"))+
  scale_fill_manual(values=c("#00BFC4", "#F8766D"), 
                    labels=c("LSP (n=31)","HSP (n=47)"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)))+
  geom_rug(data = dat_bp[dat_bp$Steroid=="Steroid_High",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#F8766D")+
  geom_rug(data = dat_bp[dat_bp$Steroid=="Steroid_Low",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#00BFC4" )+
  guides(fill=guide_legend(nrow = 2), 
         col=guide_legend(nrow=2, title = "Immune\nSubtype"))
# Plotting Heatmap and Boxplot
plot_grid(g,g1,ncol = 2,rel_widths = c(3,1), labels="AUTO")
```

Reproduction of Figure 6 at Muzzi et al. (2021)
----
```r
##### Preparing matrix of immmune modulators (IM) gene expression for Heatmap
#### Immune modulators based on Thorsson et al (2018) list 
#### C10orf54 was not found in Gene Expression Matrix, added CSF2

# Ordering by names inside each category
co_stimulator <- sort(c("CD80","CD28","ICOSLG"))
co_inhibitor <- sort(c("PDL2","PDL1","VTCN1","SLAMF7","BTN3A2","BTN3A1","B7-H3")) #CD276 = B7-H3
ligand <- sort(c("4.1BBL","TNF","OX40L","IL1B","CXCL9","CXCL10","CCL5","VEGFB","CX3CL1",
                 "TGFB1","VEGFA","CD70","CD40LG","IL10", "IFNG", "IL1A","IL12A", "IFNA2",
                 "IFNA1","IL4","IL2","IL13", "CSF2")) #TNFSF4=OX40L "TNFSF9=4.1BBL
receptor <- sort(c("TNFRSF18","TIGIT","PD1","CTLA4","IL2RA","TNFRSF4","CD27",
                   "LAG3","TNFRSF9","ICOS","BTLA","KIR2DL3","KIR2DL1",
                   "HVEM","EDNRB","CD40","ADORA2A","TLR4","HAVCR2")) # TNFRSF14 = HVEM
cell_adhesion <- sort(c("ITGB2","ICAM1","SELP"))
antigen_presentation <-
  sort(c("HLA-DRB5","HLA-DQA1","HLA-DQB1","MICA","MICB","HLA-DQA2","HLA-DQB2","HLA-B",
         "HLA-A","HLA-C","HLA-DRA","HLA-DRB1","HLA-DPB1","HLA-DPA1"))
other <- sort(c("IDO1","GZMA","PRF1","ARG1","HMGB1","ENTPD1"))

## Selecting IM genes from gene expression matrix
dat_imun <-
  gen_expt[,c(co_stimulator,co_inhibitor,ligand,receptor,cell_adhesion,antigen_presentation,other)]
cols <- colnames(dat_imun)
## Ordering by annot_hm clustered patients
dat_imun <- dat_imun[rownames(annot_hm),]

### Preparing p-value annotation for Heatmap
annot_row <- 
  data.frame(res)[c(co_stimulator,co_inhibitor, ligand, receptor,
                    cell_adhesion, antigen_presentation, other),
                  "padj",
                  drop=F]
rownames(annot_row) <- names(dat_imun)
annot_row <- -log10(annot_row)

# Z-score calculus for IM gene expression and setting maximum and minimum values
dat_imun <- data.frame(scale(dat_imun))
dat_imun[dat_imun>2]<-2
dat_imun[dat_imun<(-2)]<-(-2)
colnames(dat_imun) <- cols
# Spliting rows based on genes groups
split<- as.vector(c(rep("Co-\nStimulator",3),
                    rep("Co-\ninhibitor",7),
                    rep("Ligand",23),
                    rep("Receptor",19),
                    rep("Cell\nadhesion",3),
                    rep("Antigen\npresentation",14),
                    rep("Other",6) ))

split<- factor(split,
               levels = c("Co-\nStimulator","Co-\ninhibitor","Ligand","Receptor",
                          "Cell\nadhesion", "Antigen\npresentation","Other"))

# Top annotations
# Adding other molecular classification to annotation
annot_hm$CIMP <- clinic.dat[rownames(annot_hm),]$paper_MethyLevel
annot_hm$C1A.C1B <- clinic.dat[rownames(annot_hm),]$paper_C1A.C1B
annot_hm <- annot_hm[,c("Steroid", "CIMP", "C1A.C1B",
                        "Vital.Status","Cortisol", "Other.Hormones", "Immune.Subtype",
                        "Leukocyte.Fraction", "Lymphocyte.Inf.Score", "Macrophage.Regulation",
                        "Wound.Healing", "IFN.gamma.Response","TGF.beta.Response")]

# Heatmap construction with ComplexHeatmap                   
ha<- columnAnnotation(
              df=annot_hm[,-c(4,6,8)],
              simple_anno_size= unit(3, "mm"),
              height=unit(1,"cm"),
              annotation_name_gp = gpar(fontsize=8),
              annotation_name_side = "right",
              na_col="white", 
              annotation_legend_param = list(
                     Steroid=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       nrow=2,
                       labels=c("HSP","LSP")),
                     CIMP=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       nrow=3, 
                       labels=c("High","Intermediate", "Low")),
                     C1A.C1B=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       nrow=2),
                     Immune.Subtype=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       nrow=3, 
                       title="Immune\nSubtype"),
                     Lymphocyte.Inf.Score=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       title="Lymphocyte\nInf. Score",
                       border=T, 
                       legend_width= unit(1,"cm"),
                       direction="horizontal",
                       at=c(0,1)),
                     Macrophage.Regulation=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       title="Macrophage\nRegulation",
                       border=T,
                       legend_width= unit(1,"cm"),
                       direction="horizontal",
                       at=c(0,1)),
                     Wound.Healing=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       title="Wound\nHealing",
                       border=T,
                       legend_width= unit(1,"cm"),
                       direction="horizontal",
                       at=c(0,1)),
                     IFN.gamma.Response=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       title="IFN gamma\nResponse",
                       border=T,
                       legend_width= unit(1,"cm"),
                       direction="horizontal", 
                       at=c(0,1)),
                     TGF.beta.Response=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"),
                       title="TGF beta\nResponse",
                       border=T, 
                       legend_width= unit(1,"cm"),
                       direction="horizontal", 
                       at=c(0,1)),
                     Cortisol=list(
                       labels_gp=gpar(fontsize=8),
                       title_gp=gpar(fontsize=8, fontface="bold"))),
              col=list(
                    CIMP=c("CIMP-high"="brown2",
                             "CIMP-intermediate"="chocolate1",
                             "CIMP-low"="burlywood1"),
                    C1A.C1B=c("C1A"="plum1","C1B"="darkolivegreen2"),
                    Steroid=c("Steroid_High"="#F8766D", "Steroid_Low"="#00BFC4"),
                    Immune.Subtype=c("C1"="red","C2"="yellow","C3"="green3",
                                       "C4"="cyan","C5"="blue","C6"="pink"),
                    Cortisol=c("Cortisol"="navyblue","No"="grey"),
                    Lymphocyte.Inf.Score=circlize::colorRamp2(c(0,1), c("white","darkgreen")),
                    Macrophage.Regulation =circlize::colorRamp2(c(0,1), c("white","darkblue")),
                    Wound.Healing=circlize::colorRamp2(c(0,1), c("white","firebrick2")),
                    IFN.gamma.Response=circlize::colorRamp2(c(0,1), c("white","gold3")),
                    TGF.beta.Response=circlize::colorRamp2(c(0,1), c("white","deeppink"))
                      ))    

h1<- Heatmap(matrix = t(dat_imun),
             name = "IM Expression\nmatrix",
             show_column_names = T, 
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_title_gp = gpar(fontsize=10), 
             column_title = c("LSP (n=31)","HSP (n=47)"),
             column_split = factor(annot_hm$Steroid, levels=c("Steroid_Low","Steroid_High")),
             cluster_columns = F, 
             cluster_rows = F, 
             row_split = split, 
             row_title_rot = 0, 
             row_title_gp = gpar(fontsize=8),
             row_names_gp = gpar(fontsize=6),
             row_title_side ="left",
             heatmap_legend_param = list(
               labels_gp=gpar(fontsize=8),
               at=c(-2,0,2),
               title_gp=gpar(fontsize=8, fontface="bold"),
               direction = "horizontal",
               border=T,
               legend_width= unit(1.7,"cm")),
             right_annotation = rowAnnotation(
               df=(annot_row),
               annotation_name_rot=0,
               annotation_name_side="bottom",
               annotation_legend_param = list(
                 at=c(0,2,10),
                 labels=c("",2,">10"),
                 labels_gp=gpar(fontsize=8),
                 title_gp=gpar(fontsize=8, fontface="bold"),
                 border=T,
                 legend_width= unit(1.7,"cm"),
                 direction="horizontal",
                 nrow=2),
               annotation_label = "-log(P.Value)",
               annotation_name_gp = gpar(fontsize=8),
               col=list(padj= circlize::colorRamp2(c(0,2,10), c("white","white", "black"))
                                              )))
#### Concatenate heatmaps                 
h_list<- ha %v% h1
#### Saving object for panel construction
g<- grid.grabExpr(
  draw(h_list,
       cluster_columns=F,
       heatmap_legend_side="bottom",
       merge_legends=T),
  height = 15,
  width = 5)
```


```r 
# Boxplot Thorsson immune characteristics
## Adjustating data.frame for boxplot
dat_bp <- pivot_longer(data = annot_hm[,-c(2:6)], 
                       cols = -contains(c("Steroid","Immune.Subtype")), 
                       values_to = "value", 
                       names_to = "variable")
dat_bp$variable <- factor(dat_bp$variable,
                          levels = c("Leukocyte.Fraction", "Lymphocyte.Inf.Score", 
                                     "Macrophage.Regulation", "Wound.Healing",
                                     "IFN.gamma.Response", "TGF.beta.Response"))
dat_bp$Steroid <- factor(dat_bp$Steroid, 
                         levels = c("Steroid_Low", "Steroid_High"))
## Boxplot
g1<- ggplot(dat_bp, aes(x=Steroid, y= value))+
  facet_wrap(.~variable, 
             ncol=1, 
             scales = "free_y", 
             labeller=labeller(
               variable=c("Leukocyte.Fraction" = "Leukocyte\nFraction",
                          "Lymphocyte.Inf.Score"="Lymphocyte\nInfiltrate Score",
                          "Macrophage.Regulation"= "Macrophage\nRegulation",
                          "Wound.Healing"="Wound\nHealing", 
                          "IFN.gamma.Response"="IFN gamma\nResponse",
                          "TGF.beta.Response"="TGF beta\nResponse")))+
  geom_boxplot(aes(fill=Steroid), width=0.4, outlier.size = 0.8, lwd=0.2)+
  ggunchained::geom_split_violin(
    aes(fill=NULL, col=Immune.Subtype), 
    alpha=0.4, 
    show.legend = T)+
  scale_color_manual(values=c("green3","cyan"), 
                     labels=c("C3 (HSP n=7)\n     (LSP n=16)",
                              "C4 (HSP n=37)\n      (LSP n=12)"))+
  labs(y="Scores for immune characteristics", x= NULL) +
  theme(legend.position = "bottom",
        legend.title = element_text(face="bold", size=8),
        legend.key=element_rect(size=10, color=NA),
        legend.key.size=unit(5,"mm"),
        legend.text=element_text(size=8), 
        legend.direction = "vertical",
        legend.box = "vertical" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3),
                     label.y.npc = c(1),
                     label="p.format",
                     size=3)+
  scale_fill_manual(values=c("#00BFC4", "#F8766D"), 
                    labels=c("LSP (n=31)","HSP (n=47)"))+
  scale_x_discrete(breaks=c("Steroid_Low", "Steroid_High"),
                   labels=c("LSP", "HSP"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)))+
  geom_rug(data = dat_bp[dat_bp$Steroid=="Steroid_High",],
           aes(x=NULL), 
           sides=  "r",
           colour= "#F8766D")+
  geom_rug(data = dat_bp[dat_bp$Steroid=="Steroid_Low",],
           aes(x=NULL), 
           sides=  "l",
           colour="#00BFC4" )+
  guides(fill= guide_legend(nrow = 2),
         col= guide_legend(nrow=2, title="Immune\nSubtype"))
## Plotting Panel 
plot_grid(g, g1, ncol = 2, rel_widths = c(3,1), labels="AUTO") 

```

Statistical correlation between Leukocyte Fraction with driver gene mutations
----
```r
# Mutations and Leukocyte fraction correlation
dat <- clinic.dat[,c("Steroid", "Immune.Subtype", "Leukocyte.Fraction")] 

# cBioPortal - TCGA ACC project, mRNA patients (78) barcode for mutations in TP53, CTNNB1,MEN1 
tp53<- c("A5J2","A5J5","A5J8","A5JA","A5JB","A5JG","A5JJ","A5JL","A5JP","A5K4",
         "A5KO","A5KY","A5KZ","A5LE","A5LJ","A5HB")#Patients with TP53 mutation
# Retrieve barcodes for TP53 mutated patients
tp53 <- rownames(dat[grep(paste(tp53,collapse="|"),rownames(dat)),])

ctnnb1<- c("A5J2","A5JJ","A5JM","A5JS","A5JY","A5K5","A5K6","A5K9","A5KU","A5L3",
           "A5L8","A5LE","A5LJ")#Patients with CTNNB1 mutation
# Retrieve barcodes for CTNNB1 mutated patients
ctnnb1 <- rownames(dat[grep(paste(ctnnb1,collapse="|"),rownames(dat)),])

men1 <- c("A5JA","A5K6","A5KW","A5KY","A5LJ","A5LT","A5OF") #Patients with MEN1 mutation
# Retrieve barcodes for MEN1 mutated patients
men1 <- rownames(dat[grep(paste(men1,collapse="|"),rownames(dat)),])

#Creating columns with markers for patients with mutations
dat$TP53 <- ifelse(rownames(dat) %in% tp53,
                   "TP53",
                   "no")
dat$CTNNB1 <- ifelse(rownames(dat) %in% ctnnb1,
                     "CTNNB1",
                     "no")
dat$MEN1 <- ifelse(rownames(dat) %in% men1,
                   "MEN1",
                   "no")
dat$Leukocyte.Fraction<- as.numeric(dat$Leukocyte.Fraction)
# Statistical correlation for Leukocyte Fraction with presence of mutations 
# Wilcoxon test for Leukocyte Fraction and driver genes mutations
##TP53
wilcox.test(dat$Leukocyte.Fraction~dat$TP53)#p-value=0.18
## MEN1
wilcox.test(dat$Leukocyte.Fraction~dat$MEN1)#p-value=0.48
## CTNNB1
wilcox.test(dat$Leukocyte.Fraction~dat$CTNNB1)#p-value=0.038
## CTNNB1 only with HSP patients
wilcox.test(dat[dat$Steroid=="Steroid_High","Leukocyte.Fraction"]~
              dat[dat$Steroid=="Steroid_High","CTNNB1"])#p-value=0.99

```
 
Statistical correlation between molecular classifications, Leukocyte Fraction and Cortisol
----
```r 
# Statistical correlations between Leukocyte Fraction, molecular classifications and Cortisol
## Leukocyte Fraction and Steroid Phenotype
wilcox.test(annot_hm$Leukocyte.Fraction~annot_hm$Steroid) # p-value = 5.604e-10
## Leukocyte Fraction and C1A/C1B
wilcox.test(annot_hm$Leukocyte.Fraction~annot_hm$C1A.C1B) # p-value = 1.064e-06
## Leukocyte Fraction and Clinical Cortisol Excess
wilcox.test(annot_hm$Leukocyte.Fraction~annot_hm$Cortisol) # p-value = 0.01449
## Leukocyte Fraction and CIMP High versus CIMP Intermediate
wilcox.test(annot_hm[annot_hm$CIMP=="CIMP-high"|
                       annot_hm$CIMP=="CIMP-intermediate",]$Leukocyte.Fraction~
              annot_hm[annot_hm$CIMP=="CIMP-high"|
                         annot_hm$CIMP=="CIMP-intermediate",]$CIMP) # p-value = 0.0205
## Leukocyte Fraction and CIMP High versus CIMP Low
wilcox.test(annot_hm[annot_hm$CIMP=="CIMP-high"|
                       annot_hm$CIMP=="CIMP-low",]$Leukocyte.Fraction~
              annot_hm[annot_hm$CIMP=="CIMP-high"|
                         annot_hm$CIMP=="CIMP-low",]$CIMP) # p-value = 3.641e-06
## Leukocyte Fraction and CIMP Intermediate versus CIMP Low
wilcox.test(annot_hm[annot_hm$CIMP=="CIMP-intermediate"|
                       annot_hm$CIMP=="CIMP-low",]$Leukocyte.Fraction~
              annot_hm[annot_hm$CIMP=="CIMP-intermediate"|
                         annot_hm$CIMP=="CIMP-low",]$CIMP) # p-value = 0.002241
## Kruskal Wallis test for Leukocyte Fraction versus CIMP classification
kruskal.test(annot_hm$Leukocyte.Fraction~annot_hm$CIMP) # p-value = 1.397e-05
## Correlation between Steroid Phenotype and clinical cortisol
fisher.test(annot_hm$Steroid, annot_hm$Cortisol) # p-value = 0.0002549
## Correlation between Steroid Phenotype and C1A/C1B
fisher.test(annot_hm$Steroid, annot_hm$C1A.C1B) # p-value = 1.874e-11
## Correlation between Steroid Phenotype and CIMP
fisher.test(annot_hm$Steroid, annot_hm$CIMP) # p-value = 1.396e-06
```

Reproduction of Figure 7A at Muzzi et al. (2021)
----
```r
##### iAtlas (Thorsson et al 2018) Extracellular communication network reconstruction
##### Importing data from CRI Atlas for TCGA ACC C3 and C4 patients
##### Abundance > 0.5 and Concordance > 2.5
##### https://isb-cgc.shinyapps.io/shiny-iatlas/ Extracellular Networks

download.file(url = "https://github.com/sysbiolab/Sup_Material_Muzzi2021/blob/main/networks.RData", destfile = "networks.RData")

# Loading pre-processed network 
load("networks.RData")
redeC3 <- networks[["redeC3"]]
redeC4 <- networks[["redeC4"]]

######## Setting attributes to C3 network
# Atributes for node Type (cell, ligand, receptor)
redeC3 <- att.setv(g=redeC3,
                   from="Type",
                   to="nodeColor",
                   cols = c("darkseagreen2","khaki1","rosybrown1"))
# Atributes for node shape (for Types: cell, ligand, receptor)
redeC3 <- att.setv(g=redeC3,
                   from="Type",
                   to="nodeShape",
                   shapes =  c("ROUNDED_RECTANGLE","ELLIPSE","ELLIPSE"))
# Defining node size for each Type
redeC3 <- att.setv(g=redeC3,
                   from="Type",
                   to="nodeSize", 
                   xlim=c(50,25,5))
# Defining node line width proportional to Abundance
redeC3 <- att.setv(g=redeC3,
                   from="Abundance",
                   to="nodeLineWidth",
                   xlim=c(1,5,25))
# Defining node line color to grey
redeC3 <- att.setv(g=redeC3,
                   from="Abundance",
                   to="nodeLineColor", 
                   cols=c("grey"))
# Defining edge width proportional to Concordance
redeC3 <- att.sete(g=redeC3,
                   from="Concordance",
                   to="edgeWidth", 
                   xlim=c(1,7,10), 
                   breaks = c(2.5,5,7.5,10,12.5,15))
# Edge color = green for C3 network
V(redeC3)$edgeColor <- c("green3")
redeC3 <- att.sete(g=redeC3,
                   from="Group",
                   to="edgeColor", 
                   cols=c("green3"))

######## Setting attributes to C4 network
# Atributes for node Type (cell, ligand, receptor)
redeC4 <- att.setv(g=redeC4,
                   from="Type",
                   to="nodeColor", 
                   cols = c("darkseagreen2","khaki1","rosybrown1"))
# Atributes for node shape (for Types: cell, ligand, receptor)
redeC4 <- att.setv(g=redeC4,
                   from="Type",
                   to="nodeShape", 
                   shapes =  c("ROUNDED_RECTANGLE","ELLIPSE","ELLIPSE"))
# Defining node size for each Type
redeC4 <- att.setv(g=redeC4,
                   from="Type",
                   to="nodeSize",
                   xlim=c(50,25,5))
# Defining node line width proportional to Abundance
redeC4 <- att.setv(g=redeC4,
                   from="Abundance",
                   to="nodeLineWidth",
                   xlim=c(1,5,25))
# Defining node line color to grey
redeC4 <- att.setv(g=redeC4,
                   from="Abundance",
                   to="nodeLineColor", 
                   cols=c("grey"))
# Defining edge width proportional to Concordance
redeC4 <- att.sete(g=redeC4,
                   from="Concordance",
                   to="edgeWidth", 
                   xlim=c(1,7,10),
                   breaks = c(2.5,5,7.5,10,12.5,15))
# Edge color = green for C4 network
redeC4 <- att.sete(g=redeC4,
                   from="Group",
                   to="edgeColor",
                   cols=c("cyan"))

## The code below open a new window in RedeR to construct the network
## and then it can be exported as pdf file. To run the code below
## remove the # before the lines
#rdp <- RedPort()
#calld(rdp)
#addGraph(rdp, redeC3, layout.auto(redeC3))
#addGraph(rdp,redeC4, layout.auto(redeC4)) #changing layout to hierarchical

### Add legend based on node Type
#scl <- redeC3$legNodeColor$scale
#leg <- redeC3$legNodeColor$legend
#addLegend.color(rdp, colvec=scl,labvec=leg,title="Node Color - Type")

#resetd(rdp)
#RedeR::exitd(rdp)
```

Reproduction of Figure 7B at Muzzi et al. (2021)
----
```r 

# Heatmap of Network nodes gene expression
# Data preparation for Network Heatmap
## Retrieving pan-cancer gene expression matrix for 25 nodes from Extracellular Network
gen_exp_pan <- Pancan_data[["gen_exp_pan"]]

## Ranking by each gene expression values - lowest = 1 and highest = 9361
gen_exp_pan[,2:26] <- apply(gen_exp_pan[,2:26],
                            2,
                            rank) # 9361 patients
## Changing to friendly gene names
names(gen_exp_pan)[
  match(c("PDCD1","PDCD1LG2", "FLT1", "TNFRSF14", "TNFSF4", "CD276", "FLT1", "TNFSF9"), 
        names(gen_exp_pan))] <- 
        c("PD1","PDL2","VEGFR1", "HVEM", "OX40L", "B7-H3", "VEGFR1", "4.1BBL") 
## Retrieving total number of patients for tertiles calculus
m <- nrow(gen_exp_pan) 
## Filtering for 78 ACC patients by barcodes and ordering samples as in clinic.dat
gen_exp_pan <- gen_exp_pan[clinic.dat$patient,-1]
## Changing rownames as complete barcodes to merge with other data.frames previously constructed
rownames(gen_exp_pan) <- rownames(clinic.dat)
## Ordering data.frame by annotation clustered patients and for genes names in alphabetical order
gen_exp_pan <- gen_exp_pan[rownames(annot_hm),order(names(gen_exp_pan))]
```

```r  
# Heatmap construction with ComplexHeatmap
## Annotation columns
ha<- columnAnnotation(df=annot_hm[,c("Steroid", "Immune.Subtype")],
                      annotation_name_gp = gpar(fontsize=10),
                      height=unit(1,"cm"),
                      na_col="white",
                      annotation_name_side = "right",
                      annotation_legend_param = list(
                       Steroid=list(
                          labels_gp=gpar(fontsize=8),
                          title_gp=gpar(fontsize=8, fontface="bold"),
                          nrow=2, 
                          labels=c("HSP","LSP")),
                       Immune.Subtype=list(
                          labels_gp=gpar(fontsize=8),
                          title_gp=gpar(fontsize=8, fontface="bold"),
                          nrow=3, 
                          title="Immune\nSubtype")),
                     col=list(
                       Steroid=c("Steroid_High"="#F8766D", "Steroid_Low"="#00BFC4"),
                       Immune.Subtype=c("C1"="red","C2"="yellow","C3"="green3",
                                        "C4"="cyan","C5"="blue","C6"="pink")))    
## Main heatmap with gene expression level
h1<- Heatmap(matrix = t(gen_exp_pan),
             name = "Gene Expression\nLevel",
             show_column_names = F,
             column_split = factor(annot_hm$Steroid, levels=c("Steroid_Low","Steroid_High")),
             column_title_gp = gpar(fontsize=10), 
             column_title = c("LSP (n=31)","HSP (n=47)"),
             cluster_columns = F, 
             cluster_rows = F, 
             row_title_rot = 90, 
             row_title_gp = gpar(fontsize=8),
             row_names_gp = gpar(fontsize=8),
             row_title_side ="left",
             heatmap_legend_param = list(
               labels_gp=gpar(fontsize=8), 
               at=c(1,m/2,m),
               labels=c("Low","Interm.","High"),
               title_gp=gpar(fontsize=8, fontface="bold"),
               direction = "horizontal",
               border=T, 
               legend_width= unit(1.7,"cm")),
             col= circlize::colorRamp2(c(1,round(m/3),round(m/3)+1,round(m/3)*2,round(m/3)*2+1,m),
                                       c("gray55","gray55","khaki3","khaki1", "green3","green")))
# Concatenate heatmaps
h_list<- ha %v% h1
# Grid object
g<- grid.grabExpr(
  draw(h_list,
       cluster_columns=F, 
       heatmap_legend_side="bottom",
       merge_legends=T))
  
plot_grid(g)
```

Reproduction of Figure 8 at Muzzi et al. (2021)
----
```r
# Scatter plot CD8B X PDL1 and PDL2
## Preparing data.frame with expression values, steroid phenotype, vital status and cortisol
dat_imun <- gen_expt[,c("CD8B", "PDL1", "PDL2")]
dat_imun$Steroid <- as.factor(annot_hm[rownames(dat_imun),"Steroid"])
dat_imun$Vital_Status <- as.factor(annot_hm[rownames(dat_imun),"Vital.Status"])
dat_imun$Cortisol <- annot_hm[rownames(dat_imun),"Cortisol"]

## Creating a blank plot to compose panel
blankPlot <- ggplot()+
  geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(), 
   panel.border = element_blank(),
   panel.background = element_blank(),
   axis.title.x = element_blank(),
   axis.title.y = element_blank(),
   axis.text.x = element_blank(), 
   axis.text.y = element_blank(),
   axis.ticks = element_blank()
     )
## Retrieve adjusted p value (Wald's test) from D.E. results for CD8B
p <- format(res[rownames(res)=="CD8B",]$padj, digits=2, scientific = T)
p<- paste0("p =",p)
## Density plot for CD8B expression
xdensity <- 
  ggplot(dat_imun, aes(CD8B, fill=Steroid))+
  geom_density(alpha=.5)+
  coord_cartesian(expand=F)+
  theme(legend.position="none",
        text= element_text(size=8),
        plot.margin = margin(1,19.7,0,19.7, unit='mm'),
        axis.title.x = element_blank())+
  annotate("text",label=p, x=7.5, y=1, size=2.5)+
  scale_x_continuous(limits = c(5,9.5))
## Scatter plot for Cd8B X PDL1 expression
scatterPlot2 <- 
  ggplot(dat_imun, aes(x=CD8B, y=PDL1 ))+
  geom_point(size=1.5, aes(color=Steroid, shape=Vital_Status))+
  theme(legend.position = "none",
        text= element_text(size=8),
        plot.margin = margin(0.1,0.1,0.1,0.1,unit="pt"))+
  labs(y = "log2 PDL1 gene expression", x = element_blank())+
  geom_smooth(method="lm", formula=y~x, color="black")+
  coord_fixed(ratio = 1, expand = F, clip = "on")+
  stat_cor(method="pearson", size = 2.5, label.x = 5.5, label.y = 9.5)+
  scale_x_continuous(limits = c(5,9.5))+
  scale_y_continuous(limits = c(5,10))+
  geom_rug(aes(color=Cortisol, y=PDL1, x=NULL), show.legend = F)+
  scale_color_manual(
      values=c("Steroid_High"="#F8766D",
               "Steroid_Low"="#00BFC4",
               "Cortisol"="navyblue",
               "No"="grey"),
      na.value="white")
## Retrieve adjusted p value (Wald's test) from D.E. results for PDL1
p <- format(res[rownames(res)=="PDL1",]$padj, digits=2, scientific = T)
p<- paste0("p = ",p)
## Density plot for PDL1 expression
ydensity2 <-
  ggplot(dat_imun, aes(PDL1,fill=Steroid)) +
  geom_density(alpha=.5)+
  theme(legend.position ="none", 
        text = element_text(size=8), 
        axis.title = element_blank())+
  coord_flip(xlim=c(5,10),ylim=c(0,1.5), expand = F)+
  annotate("text",label=p, x=9.5, y=0.75, colour ="black", size=2.5)

## Scatter plot for Cd8B X PDL2 expression
scatterPlot3 <- 
  ggplot(dat_imun,aes(x=CD8B, y=PDL2))+
  geom_point(size=1.5, aes(color=Steroid, shape=Vital_Status))+
  theme(legend.position = "none",
        text= element_text(size=8), 
        plot.margin = margin(0.1,0,0.1,0,unit="pt"))+
  labs(x = "log2 CD8B gene expression", y= "log2 PDL2 gene expression")+
  geom_smooth(method="lm", formula=y~x, color="black")+
  stat_cor(method="pearson", size = 2.5, label.x = 5.5, label.y = 9.5)+
  scale_color_hue(labels=c("Steroid High (n=47)","Steroid Low (n=31)"))+
  scale_shape(labels=c("Alive (n=51)","Dead (n=27)"))+
  coord_fixed(ratio = 1, expand = F, clip = "on")+
  scale_x_continuous(limits = c(5,9.5))+
  scale_y_continuous(limits = c(5,10))
## Retrieve scatterPlot3 legend for Steroid and Vital Status
leg <- get_legend(
  scatterPlot3+ 
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold"),
        legend.direction = "vertical",
        legend.box = "horizontal", 
        legend.text= element_text(size=7),
        legend.spacing.y= unit(0.2,"mm")))
## Create legend for Cortisol levels
leg1<- get_legend(
  ggplot(dat_imun)+
    geom_rug(aes(color=Cortisol), show.legend = T)+
    scale_color_manual(values=c("Cortisol"="navyblue","No"="grey"),
                       na.value="white",
                       labels=c("Cortisol (n=32)", "No (n=41)", "NA (n=5)"))+
    theme(legend.position = "bottom", 
          text= element_text(size=8), 
          legend.title = element_text(face="bold"), 
          legend.direction = "vertical",
          legend.box = "horizontal", 
          legend.text= element_text(size=7),
          legend.spacing.y= unit(0.2,"mm")))
## Adding rug plot in scatterPlot3
scatterPlot3 <- 
  scatterPlot3+
  geom_rug(aes(color=Cortisol), show.legend = F)+
  scale_color_manual(
    values=c("Steroid_High"="#F8766D", 
             "Steroid_Low"="#00BFC4",
             "Cortisol"="navyblue",
             "No"="grey"),
    na.value="white")
## Retrieve adjusted p value (Wald's test) from D.E. results for PDL2
p <- format(res[rownames(res)=="PDL2",]$padj, digits=2, scientific = T)
p<- paste0("p = ",p)
## Density plot for PDL2 expression
ydensity3 <- 
  ggplot(dat_imun, aes(PDL2,fill=Steroid)) +
  geom_density(alpha=.5)+
  theme(legend.position="none",
        text= element_text(size=8), 
        axis.title.y  = element_blank())+
  coord_flip(xlim=c(5,10),ylim=c(0,1.5), expand = F)+
  annotate("text",label=p, x=9.5, y=0.75, colour ="black", size=2.5)

## Aligning plots
row0<- plot_grid(xdensity, blankPlot, ncol=2, rel_widths = c(1.2,0.3))
row1<- plot_grid(scatterPlot2, ydensity2 ,ncol=2, align = "h", rel_widths = c(1.2,0.3))
row2<- plot_grid(scatterPlot3, ydensity3 ,ncol=2, align = "h", rel_widths = c(1.2,0.3))
row3<- plot_grid(leg,leg1, ncol=2, align = "hv", rel_widths = c(2,2))
## Plotting panel
plot_grid(row0,row1, row2, row3, ncol=1, rel_heights = c(.5,2,2,.7) , align = "hv")
```

Reproduction of Figure 9 at Muzzi et al. (2021)
----
```r
#### Survival analysis
## Retrieving survival data from Thorsson et al (2018) data set
dat_imun <- clinic.dat[,c("Steroid", "Immune.Subtype", "PFI", "PFI.Time", "OS", "OS.Time")]

## Kaplan-Meier for Steroid Phenotypes
### Overall Survival (OS)
fit <- survfit(Surv(time=OS.Time, event = OS)~Steroid, data=dat_imun)
g1 <- ggsurvplot(fit,
                 data=dat_imun,
                 title="Overall Survival",
                 xlab="Time in days",
                 legend.title="",
                 legend.labs=c("HSP (n=47)", "LSP (n=31)"),
                 conf.int = T, 
                 conf.int.style= "ribbon",
                 conf.int.alpha=0.15,
                 pval = T,
                 pval.size=3, 
                 risk.table = TRUE,
                 risk.table.y.text=F, 
                 risk.table.fontsize=3, 
                 risk.table.title="",
                 ggtheme = theme_grey() )
g1$plot <- ggpar(g1$plot, font.x = 0, font.title = 10)
g1$table <- ggpar(g1$table, font.y = 0)
### Progression Free Interval (PFI)
fit <- survfit(Surv(time=PFI.Time, event = PFI)~Steroid, data=dat_imun)
g2<-ggsurvplot(fit,
               data=dat_imun,
               title="Progression Free Interval", 
               xlab="Time in days",
               legend.labs=c("HSP (n=47)", "LSP (n=31)"), 
               legend.title="",
               conf.int = T, 
               conf.int.style= "ribbon",
               conf.int.alpha=0.15, 
               pval = T,
               pval.size=3, 
               risk.table = TRUE,
               risk.table.y.text=F,
               risk.table.fontsize=3,
               risk.table.title="",
               ggtheme = theme_grey() )
g2$plot <- ggpar(g2$plot, font.x = 0, font.title = 10)
g2$table <- ggpar(g2$table, font.y = 0)   
### Saving plots for Steroid Phenotype Survival analysis
row1<- plot_grid(g1$plot, 
                 g2$plot,
                 g1$table, 
                 g2$table,
                 ncol = 2, 
                 rel_heights = c(1,0.5), 
                 align="v")

## Kaplan-Meier for Immune Subtypes
### Selecting C3 and C4 subtypes for comparison (too few participants in the other groups)
dat_imun <- dat_imun[dat_imun$Immune.Subtype=="C3" | dat_imun$Immune.Subtype=="C4",]
### Overall Survival (OS)
fit <- survfit(Surv(time=OS.Time, event = OS)~Immune.Subtype, data=dat_imun)
g1<-ggsurvplot(fit,
               data=dat_imun,
               title="Overall Survival", 
               xlab="Time in days",
               legend.title="",
               legend.labs=c("C3 (n=23)" , "C4 (n=49)"),
               palette = c("green3","cyan"), 
               conf.int = T,
               conf.int.style= "ribbon", 
               conf.int.alpha=0.15,
               pval = T,
               pval.size=3,
               risk.table = T,
               risk.table.y.text=F,
               risk.table.fontsize=3, 
               risk.table.title="",
               ggtheme = theme_grey())
g1$plot <- ggpar(g1$plot, font.x = 0, font.title = 10)
g1$table <- ggpar(g1$table, font.y = 0)

### Progression Free Interval (PFI)
fit <- survfit(Surv(time=PFI.Time, event = PFI)~Immune.Subtype, data=dat_imun)
g2<-ggsurvplot(fit,
               data=dat_imun,
               title="Progression Free Interval",
               xlab="Time in days",
               legend.title="", 
               legend.labs=c("C3 (n=23)" , "C4 (n=49)"),
               palette = c("green3","cyan"), 
               conf.int = T,
               conf.int.style= "ribbon",
               conf.int.alpha=0.15,
               pval = T,
               pval.size=3,
               risk.table = TRUE,
               risk.table.y.text=F,
               risk.table.fontsize=3, 
               risk.table.title="",
               ggtheme = theme_grey() )
g2$plot <- ggpar(g2$plot, font.x = 0, font.title = 10)
g2$table <- ggpar(g2$table, font.y = 0)  
### Saving plots for Immune Subtypes Survival analysis
row2<- plot_grid(g1$plot,
                 g2$plot, 
                 g1$table, 
                 g2$table, 
                 ncol = 2, 
                 rel_heights = c(1,0.5),
                 align="v")

## Kaplan-Meier for LSP C3 X Others
### Preparing data for analysis
dat_imun <- clinic.dat[,c("Steroid", "Immune.Subtype", "PFI", "PFI.Time", "OS", "OS.Time")]
### Merging C3 with LSP groups
dat_imun$Immune.Subtype <- 
  factor(if_else(dat_imun$Immune.Subtype == "C3" & dat_imun$Steroid == "Steroid_Low", 
                 "LSP.C3",
                 "Other" ),
         levels = c("LSP.C3","Other"))
### Overall Survival (OS)
fit <- survfit(Surv(time=OS.Time, event = OS)~Immune.Subtype, data=dat_imun)
g1 <- ggsurvplot(fit,
                 data=dat_imun,
                 title="Overall Survival", 
                 xlab="Time in days",
                 legend.title="", 
                 legend.labs=c("LSP C3 (n=16)", "Others (n=62)"),
                 palette = c("darkgreen","chocolate1"),
                 conf.int = T, 
                 conf.int.style= "ribbon", 
                 conf.int.alpha=0.15,
                 pval = T,
                 pval.size=3,
                 risk.table = TRUE, 
                 risk.table.y.text=F, 
                 risk.table.fontsize=3, 
                 risk.table.title="", 
                 ggtheme = theme_grey())
g1$plot <- ggpar(g1$plot, font.x = 0, font.title = 10)
g1$table <- ggpar(g1$table, font.y = 0) 

### Progression Free Interval (PFI)
fit <- survfit(Surv(time=PFI.Time, event = PFI)~Immune.Subtype, data=dat_imun)
g2<-ggsurvplot(fit,
               data=dat_imun,
               title="Progression Free Interval",
               xlab="Time in days",
               legend.title="",
               legend.labs=c("LSP C3 (n=16)", "Others (n=62)"),
               palette = c("darkgreen","chocolate1"), 
               conf.int = T, 
               conf.int.style= "ribbon", 
               conf.int.alpha=0.15, 
               pval = T, 
               pval.size=3,
               risk.table = TRUE, 
               risk.table.y.text=F,
               risk.table.fontsize=3, 
               risk.table.title="", 
               ggtheme = theme_grey())
g2$plot <- ggpar(g2$plot, font.x = 0, font.title = 10)
g2$table <- ggpar(g2$table, font.y = 0)   
### Saving plots for LSP C3 Survival analysis
row3<- plot_grid(g1$plot,
                 g2$plot,
                 g1$table, 
                 g2$table, 
                 ncol = 2, 
                 rel_heights = c(1,0.5),
                 align="v")
## Plotting all panels
plot_grid(row1,row2,row3,ncol=1, labels = "AUTO")
```

Reproduction of Table 1 at Muzzi et al. (2021)
----
```r
# Clinical table
tab<- clinic.dat[order(clinic.dat$patient),
                 c("patient", "age_at_index", "gender", "tumor_stage", "vital_status", 
                   "Steroid","Proliferation_mRNA","Cortisol", "Other.Hormones", "Immune.Subtype")]
colnames(tab)[2] <- "Age"
## Group by Steroid and Vital Status and count patients for each category
tab <- tab %>% group_by(Steroid,vital_status)

tab1 <- tab %>% dplyr::select(-contains(c("patient","Age"), ignore.case = F))
tab1 <- tab1 %>% pivot_longer(cols=-c("Steroid","vital_status"),
                              names_to="name",
                              values_to="value" )

tab1 <- 
  tab1%>% 
  dplyr::count(value, .drop = F) %>% 
  pivot_wider(names_from = value, 
              values_from=n,
              values_fill=0 )
# Create a TOTAL count for Steroid and Vital.Status
tab2 <- tab %>% summarise(TOTAL=n())

# Unify total counts with counts for each category
tab1<- cbind(tab1[,c("Steroid", "vital_status")], 
             tab2[,"TOTAL"],
             tab1[,c("male","female", "Proliferation_High",
                     "stage i","stage ii","stage iii","stage iv",
                     "Cortisol", "Mineralcorticoids","Sexual",
                     "C1", "C2", "C3", "C4","C5", "C6")])
# Sum dead and alive and adjusting table               
row1 <- paste0( (tab1[1,3:ncol(tab1)]+tab1[2,3:ncol(tab1)]),
                " (",
                tab1[2,3:ncol(tab1)], ")")   
row2 <- paste0( (tab1[3,3:ncol(tab1)]+tab1[4,3:ncol(tab1)]),
                " (",
                tab1[4,3:ncol(tab1)], ")")   
tab1 <- data.frame(tab1)
tab1[1,3:ncol(tab1)] <- row1             
tab1[3,3:ncol(tab1)] <- row2              
tab1 <- tab1[c(1,3),]
tab1 <- rbind(names(tab1),tab1)
# Adding mean age with standard deviation
tab1[,2] <- c("Mean Age \u00b1SD",
               paste0(round(mean(tab[tab$Steroid!="Steroid_Low",]$Age),digits=1),
                      "\u00b1",
                      round(sd(tab[tab$Steroid!="Steroid_Low",]$Age),digits=1)),
               paste0(round(mean(tab[tab$Steroid=="Steroid_Low",]$Age),digits=1),
                      "\u00b1",
                      round(sd(tab[tab$Steroid=="Steroid_Low",]$Age),digits=1)) )
tab1[,1] <- c(" ", "HSP","LSP")
tab1[,4] <- c("Males\nMean age\u00b1SD",
              paste0(tab1[2,4],"\n",
                     round(mean(tab[tab$Steroid!="Steroid_Low" &
                                    tab$gender=="male",]$Age), digits=1),
                     "\u00b1",
                     round(sd(tab[tab$Steroid!="Steroid_Low" &
                                  tab$gender=="male",]$Age), digits=1)),
              paste0(tab1[3,4],"\n",
                     round(mean(tab[tab$Steroid=="Steroid_Low" &
                                    tab$gender=="male",]$Age), digits=1),
                     "\u00b1",
                     round(sd(tab[tab$Steroid=="Steroid_Low" &
                                  tab$gender=="male",]$Age), digits=1)))
tab1[,5] <- c("Females\nMean age\u00b1SD",
              paste0(tab1[2,5],"\n",
                     round(mean(tab[tab$Steroid!="Steroid_Low" &
                                    tab$gender=="female",]$Age),digits=1),
                     "\u00b1",
                     round(sd(tab[tab$Steroid!="Steroid_Low" &
                                  tab$gender=="female",]$Age),digits=1)),
              paste0(tab1[3,5],"\n",
                     round(mean(tab[tab$Steroid=="Steroid_Low" &
                                    tab$gender=="female",]$Age),digits=1),
                     "\u00b1",
                     round(sd(tab[tab$Steroid=="Steroid_Low" &
                                  tab$gender=="female",]$Age),digits=1)))
# Define column names and order
tab1[1,]<- c("Steroid","Mean Age \u00b1SD","Total (Dead)",
             "Males\nMean age\u00b1SD","Females\nMean age\u00b1SD",
             "High Proliferation",
             "Stage I","Stage II","Stage III","Stage IV",
             "Cortisol","Mineralcorticoids","Sexual",
             "C1", "C2", "C3", "C4","C5", "C6")
tab1 <- cbind(tab1[,1:6], c("Tumor Stage","",""), 
              tab1[,7:10],
              c("Hormone Excess","",""),
              tab1[,11:13],
              c("Immune Subtype","",""),
              tab1[,14:19])
# Separate HSP and LSP cells for first row of final table
col <- tableGrob(t(tab1[2:3,1]),
                 rows=NULL,
                 cols=c("HSP","LSP"),
                 theme = ttheme_default(core=list(fg_params=list(fontsize=8)),
                                        colhead=list(fg_params=list(fontsize=8))))
# Preparing table in tableGrob format with clinical characteristics
tab1<- tableGrob(t(tab1[,-1]),
                 rows = NULL,
                 cols=NULL,
                 theme = ttheme_default(core=list(fg_params=list(fontsize=8))))

# Adding HSP and LSP cells in the first row of the table
tab <- gtable_combine(col[1,], tab1, along=2)
# Defining cells width and height
tab$widths <- rep(max(tab$widths), length(tab$widths))
tab$height <- tab$height*0.8
# Adjusting first row
tab$layout[1:4,c("l","r")] <- list(c(2,3),c(2,3))

#### Changing fontface 

#Discover index of cells in table (t=row and l=column)
which(tab$layout$t==21 & tab$layout$l==1 & tab$layout$name=="core-fg") #row 2:22 column 1 
# index 5:25
## Changing fontface to bold
for(i in 5:25){
  tab$grobs[i][[1]][["gp"]] <- gpar(fontface="bold", fontsize=8)
}
which(tab$layout$t==3 & tab$layout$l==3 & tab$layout$name=="core-fg") #row 3 column 2 and 3 
#index 27 and 48
tab$grobs[27][[1]][["gp"]] <- gpar(fontface="bold", fontsize=8)
tab$grobs[48][[1]][["gp"]] <- gpar(fontface="bold", fontsize=8)

### Changing background colour
which(tab$layout$t==22 & tab$layout$l==1 & tab$layout$name=="core-bg") #row 2:22 column 1 
#68:88
for(i in 68:88){
  tab$grobs[i][[1]][["gp"]] <- gpar(fill="grey90", lwd=0, col="white")
}

which(tab$layout$t==7 & tab$layout$l==2 & tab$layout$name=="core-bg") #row 7 columns 1:3
#Tumor Stage # index 73, 94, 115
tab$grobs[73][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")
tab$grobs[94][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")
tab$grobs[115][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")

which(tab$layout$t==12 & tab$layout$l==1 & tab$layout$name=="core-bg") #row 12 columns 1:3
#Hormone Excess #index 78,99,120
tab$grobs[78][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")
tab$grobs[99][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")
tab$grobs[120][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")

which(tab$layout$t==16 & tab$layout$l==3 & tab$layout$name=="core-bg") #row 16 columns 1:3
#Immune Subtype #index 82, 103, 124
tab$grobs[82][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")
tab$grobs[103][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")
tab$grobs[124][[1]][["gp"]] <- gpar(fill="grey77", lwd=0, col="grey77")

plot_grid(tab)
```

R session
----
```r
sessionInfo()
```

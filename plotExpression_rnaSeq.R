# Script to analyze mValidation gene distribution + plot
## Use CRAN package librarian to install / update / attach R packages
#pkg <- "librarian"[!("librarian" %in% installed.packages()[,"Package"])]
#if(length(pkg)) install.packages(pkg)
#library(librarian)
#shelf(reshape2, plyr, dplyr, fgsea, tidyr, ggplot2, ggpubr, xlsx, quiet=TRUE)
cran_pkgs <- c("reshape2", "plyr", "dplyr", "tidyr", "ggplot2", "ggpubr", "xlsx")
new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[,"Package"])]
if(length(new_cran)) install.packages(new_cran, repos="https://cran.r-project.org")

bioc_pkgs <- "fgsea"
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc)

suppressMessages({
  library(reshape2)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(xlsx)
  library(fgsea)
})

## USER INPUT ##
dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forMichael_rnaSeq", dataDir)

# Output directory
dt <- format(Sys.Date(), "20%y%m%d")
outDir  <- sprintf("%s/res/%s_out_forMichael_rnaSeq", dataDir, dt)
if (!file.exists(outDir)) dir.create(outDir)

# Paths to mouse pathway annotation files (interested in autophagy)
pathGOF <- sprintf("%s/anno/Human_GOBP_AllPathways_no_GO_iea_October_01_2018_symbol.gmt", dataDir)
pathGO <- gmtPathways(pathGOF)

# Manually curated pathways
pathCurate <- read.xlsx(sprintf("%s/HumanCyc_LipidMetabolism.xlsx", geneDir), sheetIndex=1)
pathCurate <- split(pathCurate, pathCurate$Pathway, drop=TRUE)

## START ANALYSIS ##
datF <- list.files(pattern="RNAseq.*csv", path=geneDir, full.names=T)
dat <- read.csv(datF)

plotStuff <- function(set, anno, pathID) {

  # Get data of interest
  if (set == "logFPKM") {
    dat2 <- dat[,c(1,23:28)] # logFPKM columns
    dat_plot <- melt(dat2) # melt to table
    dat_plot2 <- separate(dat_plot, variable, c("KO_Screen", "Experiment"), sep="_")
  }
  if (set == "logFC") {
    dat2 <- dat[,c(1,29:31)] # logFC columns
    dat_plot <- melt(dat2)  # melt to table
    dat_plot2 <- separate(dat_plot, variable, c("Experiment", "KO_Screen"), sep="_")
  }

  dat_plot2 <- na.omit(dat_plot2) # remove missing values
  colnames(dat_plot2)[1] <- "Gene" # rename column

  # Annotating genes in metabolic pathways
  if (anno == "GOBP") {
    path <- pathGO[grep(pathID, names(pathGO))] # get pathway corresponding to ID
    path2 <- as.character(unlist(path))         # vector of pathway genes
    plot_name <- sprintf("%s/gene_%s_%s_%s.pdf", outDir, set, anno, pathID) # plot file name
    title <- sub("\\%.*", "", names(path))    # plot title
  }

  if (anno == "manual") {
    path <- pathCurate[which(names(pathCurate) %in% pathID)] # get pathway corresponding to path name
    path <- melt(path)                                       # melt to table
    path2 <- as.character(path$Gene)                         # vector of pathway genes
    path_name <- gsub(" ", "_", as.character(unique(path$Pathway))) # replace spaces with underscores in pathway name
    path_source <- as.character(unique(path$Source))              # source of pathway
    plot_name <- sprintf("%s/gene_%s_%s_%s.pdf", outDir, set, path_name, path_source) # plot file name
    title <- as.character(unique(path$Pathway))                                   # plot title
  }

  dat_plot2$Pathway <- NA # indicate which genes are in pathway
  dat_plot2$Pathway <- ifelse(dat_plot2$Gene %in% path2, "Metabolic pathway", NA)
  dat_plot2 <- na.omit(dat_plot2) # remove genes not in pathway

  # Control ordering of boxplots + facet order
  dat_plot2$KO_Screen <- factor(dat_plot2$KO_Screen, levels=c("WT", "C12orf49", "SREBF2"))

  if (set == "logFPKM") {
    dat_plot2$Experiment <- factor(dat_plot2$Experiment, levels=c("plus", "minus"))
  }

  # Plot lines
  p <- ggplot(dat_plot2, aes(x=KO_Screen, y=value)) +
          geom_hline(yintercept=0, linetype="dashed", size=0.5) + # dashed line at 0
          geom_boxplot(aes(fill=KO_Screen)) +                     # color boxplots by gene
          labs(y=set, x=NULL, title=title) +                      # labels
          theme_bw() +                                            # black / white theme
          scale_fill_brewer(palette="Set2") +                     # colour palette
          theme(text=element_text(family="sans", size=14),        # text family / size
                panel.grid=element_blank(),                       # remove border lines
                panel.border=element_blank(),
                axis.text.x=element_blank(),                      # remove x axis text
                axis.ticks.x=element_blank(),
                axis.line=element_line(colour="black"),
                strip.text=element_text(face="bold"),             # bold facet labels
                legend.position="bottom",                         # legend at bottom
                legend.title=element_blank())                     # remove legend title

  if (set == "logFPKM") {
    p <- p + facet_grid(.~Experiment)  # facet plots by 'minus' / 'plus' when plotting logFPKM
  }

  # Calculate p value (using WT as reference group)
  p <- p + stat_compare_means(label="p.signif", method="t.test", ref.group="WT")

  # Extract ggplot data
  pg <- ggplot_build(p)
  pg_p <- pg$data[[3]]
  write.table(pg_p, sprintf("%s_pvals.txt", plot_name), col=T, row=F, quote=F, sep="\t")

  ggsave(plot_name, p, width=3.5, height=5.5, dpi=300)
}

# Plot pathway boxplots
## GO IDs
set_list <- c(rep("logFPKM", 4), rep("logFC", 4))
pathID_list <- rep(c("0006695", "0045540", "0006629", "0008610"), 2)
mapply(plotStuff, set_list, "GOBP", pathID_list)

## Manually curated pathways (10 pathways)
set_list <- c(rep("logFPKM", 10), rep("logFC", 10))
pathID_list <- rep(names(pathCurate), 2)
mapply(plotStuff, set_list, "manual", pathID_list)

## Alternatively, can call function like this for each iteration, without mapply
plotStuff(set="logFPKM", anno="GOBP", pathID="0006695")




## --- Set working directory --- #
setwd("~/Desktop/Master_VHIR/TFM/")

## --- Load libraries, if necessary --- ##
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(RColorBrewer)

## ---------------------------------- DOMAINS ---------------------------------- ##
## --- MOSQUITO --- ##
## --- Load data --- ##
culex.d <- read_xlsx("Domains_culex.xlsx")

## compute cpm and relative abundance

culex.d$cpm <- (culex.d$`Number of reads`/culex.d$`Total reads`)* 10^6
culex.d$abundance <- (culex.d$`Number of reads`/culex.d$`Total reads`)* 100
culex.d <- culex.d[culex.d$abundance > 0.1,]

ggplot(culex.d, aes(x = as.factor(Sample), y = abundance, fill =Domain, col=Domain)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Relative abundance (%)") +
  theme(legend.title = element_text(), legend.position = "top", axis.text=element_text(size=15), axis.title = element_text(size=15), legend.text = element_text(size=12)) +
 labs(fill="Domains") + 
  scale_fill_manual(values = c("#6698FF", "#153E7E", "#98AFC7", "#C0C0C0")) + 
  scale_color_manual(values = c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) + guides(col=FALSE)

## --- FISH --- ##
## --- Load data --- ##
fish.d <- read_xlsx("domains_gambusia.xlsx")

## compute cpm and relative abundance

fish.d$cpm <- (fish.d$`Number of reads`/fish.d$`Total reads`)* 10^6
fish.d$abundance <- (fish.d$`Number of reads`/fish.d$`Total reads`)* 100
fish.d <- fish.d[fish.d$abundance > 0.1,]

ggplot(fish.d, aes(x = as.factor(Sample), y = abundance, fill =Domain, col=Domain)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Relative abundance (%)") +
  theme(legend.title = element_text(), legend.position = "top", axis.text=element_text(size=15), axis.title = element_text(size=15), legend.text = element_text(size=12)) +
  labs(fill="Domains")  + scale_fill_manual(values = c("#FF0000","#6698FF", "#153E7E", "#98AFC7", "#C0C0C0")) + 
  scale_color_manual(values = c("#FF0000","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF", "#FFFFFF")) + guides(col=FALSE)

## ---------------------------------- VIROME ---------------------------------- ##
## --- MOSQUITO --- ##
## --- Load data --- ##
culex.v <- read_xlsx("culex_virome.xlsx")

## compute cpm and relative abundance

culex.v$cpm <- (culex.v$Reads/culex.v$`Total reads`)* 10^6
culex.v$abundance <- (culex.v$Reads/culex.v$`Total reads`)* 100
culex.v <- culex.v[culex.v$abundance > 0.1,]

# Define the number of colors you want
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(culex.v, aes(x = as.factor(Sample), y = abundance, fill =Family)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Relative abundance (%)") +
  theme(legend.title = element_text(), legend.position = "top", axis.text=element_text(size=15), axis.title = element_text(size=15), legend.text = element_text(size=12)) + 
 labs(fill="Family") +  scale_fill_manual(values = mycolors)


## --- MOSQUITO --- ##
## --- Load data --- ##
fish.v <- read_xlsx("virome_fish.xlsx")

## compute cpm and relative abundance

fish.v$cpm <- (fish.v$Reads/fish.v$`Total reads`)* 10^6
fish.v$abundance <- (fish.v$Reads/fish.v$`Total reads`)* 100
fish.v <- fish.v[fish.v$abundance > 0.1,]

# Define the number of colors you want
nb.cols <- 30
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(fish.v, aes(x = as.factor(Sample), y = abundance, fill =Family)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Relative abundance (%)") +
  theme(legend.title = element_text(), legend.position = "top", axis.text=element_text(size=15), axis.title = element_text(size=15), legend.text = element_text(size=12)) + 
  labs(fill="Family") +  scale_fill_manual(values = mycolors)



## ------------ SPECIES VIROME GRAPHS ------------ ##
bu.all <- read_xlsx("Culex_species.xlsx")

bu.all$cpm <- (bu.all$Reads/bu.all$`Total reads`)* 10^6
bu.all$abundance <- (bu.all$Reads/bu.all$`Total reads`)* 100
bu.all <- bu.all[bu.all$abundance > 0.1,]

nb.cols <- 6
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(bu.all, aes(x = as.factor(Sample), y = abundance, fill =Species, col=Species)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Relative abundance (%)") +
  theme(legend.title = element_text(), legend.position = "top", axis.text=element_text(size=15), axis.title = element_text(size=15), legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(fill="Family") +  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF", "#CCCCCC")) + guides(col=FALSE)



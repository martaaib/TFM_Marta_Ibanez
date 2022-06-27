#########################################
##  R version 4.0.4 (2021-02-15)
##
## Bioinformatics pipeline graphs 
##
## Copyright: Marta Ibañez Lligoña (marta.ibanez@vhir.org)
##
#########################################


## --- set working directory --- ##
setwd("~/Desktop/Master_VHIR/TFM")

## --- Load packages --- ##
library(ggplot2)
library(readxl)
library(tidyr)
library(ggbreak)
library(scales)
library(stringr)

## --- REMOVAL OF HOST GENOMES --- ##
## --- Load data --- ##
culex <- read_xlsx("host_genome_mosquito.xlsx")
affinis <- read_xlsx("host_genome_fish.xlsx")


## --- Transform to long format --- ##
## --- Culex sp. --- ##
culex$CP <- as.numeric(culex$CP)
culex$CQ <- as.numeric(culex$CQ)
culex$HS <- as.numeric(culex$HS)
culex$unknown <- 100-(culex$CP+ culex$HS + culex$CQ)
colnames(culex) <- c("Sample", "C. Pipiens", "C. Quinquefasciatus", "Human", "Others")
culex.long <- gather(culex, genome, percentage , `C. Pipiens`:Others, factor_key=TRUE)
culex.long <- culex.long[order(culex.long$Sample),]

## --- PLOT --- ##
ggplot(culex.long, aes(x = as.factor(Sample), y = percentage, fill=genome, col=genome)) + geom_col(position = "dodge") + theme_bw()  + xlab("Samples") + ylab("Percentage of reads") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="") + 
  scale_fill_manual(values = c("mediumseagreen", "lightseagreen", "turquoise3", "navyblue")) +scale_color_manual(values =c("seagreen4", "seagreen", "turquoise4", "navyblue")) +
  guides(col=FALSE)

## --- G. affinis --- ##
affinis$GA <- as.numeric(affinis$GA)
affinis$unknown <- 100-affinis$GA
colnames(affinis) <- c("Sample", "G. affinis", "Others")
affinis.long <- gather(affinis, genome, percentage , `G. affinis`:Others, factor_key=TRUE)
affinis.long <- affinis.long[order(affinis.long$Sample),]

## --- PLOT --- ##
ggplot(affinis.long, aes(x = as.factor(Sample), y = percentage, fill=genome, col=genome)) + geom_col(position = "dodge") + theme_bw()  + xlab("Samples") + ylab("Percentage of reads") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="") + scale_fill_manual(values = c("lightsteelblue",  "navyblue")) +scale_color_manual(values = c("lightblue4",  "navyblue")) +
  guides(col=FALSE)


## --- Sensitivity study --- ##
## --- Load data --- ##
culex.s <- read_xlsx("host_genome_mosquito_SENS.xlsx")
affinis.s <- read_xlsx("host_genome_fish_SENS.xlsx")


## --- Transform to long format --- ##
## --- Culex sp. --- ##
culex.s$CP <- as.numeric(culex.s$CP)
culex.s$CQ <- as.numeric(culex.s$CQ)
culex.s$HS <- as.numeric(culex.s$HS)
culex.s$unknown <- 100-(culex.s$CP+ culex.s$HS + culex.s$CQ)
colnames(culex.s) <- c("Sample", "C. Pipiens", "C. Quinquefasciatus", "Human", "Others")
culex.long2 <- gather(culex.s, genome, percentage , `C. Pipiens`:Others, factor_key=TRUE)
culex.long2 <- culex.long2[order(culex.long2$Sample),]

## --- PLOT --- ##
ggplot(culex.long2, aes(x = as.factor(Sample), y = percentage, fill=genome, col=genome)) + geom_col(position = "dodge") + theme_bw()  + xlab("Samples") + ylab("Percentage of reads") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="") + 
  scale_fill_manual(values = c("mediumseagreen", "lightseagreen", "turquoise3", "navyblue")) +scale_color_manual(values =c("seagreen4", "seagreen", "turquoise4", "navyblue")) +
  guides(col=FALSE)

## --- G. affinis --- ##
affinis.s$GA <- as.numeric(affinis.s$GA)
affinis.s$unknown <- 100-affinis.s$GA
colnames(affinis.s) <- c("Sample", "G. affinis", "Others")
affinis.long2 <- gather(affinis.s, genome, percentage , `G. affinis`:Others, factor_key=TRUE)
affinis.long2 <- affinis.long2[order(affinis.long2$Sample),]

## --- PLOT --- ##
ggplot(affinis.long2, aes(x = as.factor(Sample), y = percentage, fill=genome, col=genome)) + geom_col(position = "dodge") + theme_bw()  + xlab("Samples") + ylab("Percentage of reads") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="") + scale_fill_manual(values = c("lightsteelblue",  "navyblue")) +scale_color_manual(values = c("lightblue4",  "navyblue")) +
  guides(col=FALSE)


### ---- REMOVAL OF DUPLICATED AND LOW-QUALITY READS --- ###
## --- MOSQUITO SAMPLES --- ##
culex.reads <- read_xlsx("removal_reads_mosquito.xlsx")
head(culex.reads)
culex.reads_long <- gather(culex.reads, Reads, number , Initial_reads:Classified_reads, factor_key=TRUE)
culex.reads_long$Reads <- str_replace_all(culex.reads_long$Reads, "Initial_reads","Initial number of reads" ) 
culex.reads_long$Reads <- str_replace_all(culex.reads_long$Reads, "after_trim", paste0( "Number of reads after removal", "\n",  "of duplications and low quality reads"))
culex.reads_long$Reads <- str_replace_all(culex.reads_long$Reads, "after_mapping","Number of reads after removal of host genome" ) 
culex.reads_long$Reads <- str_replace_all(culex.reads_long$Reads, "Classified_reads","Number of reads classified as taxa" ) 

ggplot(culex.reads_long, aes(x = as.factor(Sample), y = number, fill=Reads, col=Reads)) + geom_col(position = "dodge") + theme_bw() + scale_y_break(c(4533209, 25562418), scales=0.2) +
  scale_y_continuous(name ="Number of reads", breaks=c(0,200000,500000,1e+06, 2e+06,3e+06, 4e+06, 3e+07, 5e+07, 6e+07), position = "left") +xlab("Samples") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=11), axis.title = element_text(size=11), legend.text = element_text(size=11), axis.text.y.right  = element_blank(), axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="")  + guides(col=FALSE, fill= guide_legend(ncol = 2, nrow = 2, byrow = TRUE)) +
  scale_color_manual(values =c("skyblue4", "darkgoldenrod", "forestgreen", "tomato3")) + scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#D55E00"))
  
## --- FISH SAMPLES --- ##
## Read and prepare data
fish.reads <- read_xlsx("removal_reads_fish.xlsx")
head(fish.reads)
fish.reads_long <- gather(fish.reads, Reads, number , Initial_reads:Classified_reads, factor_key=TRUE)
fish.reads_long$Reads <- str_replace_all(fish.reads_long$Reads, "Initial_reads","Initial number of reads" ) 
fish.reads_long$Reads <- str_replace_all(fish.reads_long$Reads, "after_trim", paste0( "Number of reads after removal", "\n",  "of duplications and low quality reads"))
fish.reads_long$Reads <- str_replace_all(fish.reads_long$Reads, "after_mapping","Number of reads after removal of host genome" ) 
fish.reads_long$Reads <- str_replace_all(fish.reads_long$Reads, "Classified_reads","Number of reads classified as taxa" ) 

## PLOT
ggplot(fish.reads_long, aes(x = as.factor(Sample), y = number, fill=Reads, col=Reads)) + geom_col(position = "dodge") + theme_bw() + scale_y_break(c(4533209, 25562418), scales=0.2) +
  scale_y_continuous(name ="Number of reads", breaks=c(0,200000,500000,1e+06, 2e+06,3e+06, 4e+06, 3e+07, 5e+07, 6e+07), position = "left") +xlab("Samples") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=11), axis.title = element_text(size=11), legend.text = element_text(size=11), axis.text.y.right  = element_blank(), axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="")  + guides(col=FALSE, fill= guide_legend(ncol = 2, nrow = 2, byrow = TRUE)) +
  scale_color_manual(values =c("skyblue4", "darkgoldenrod", "forestgreen", "tomato3")) + scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#D55E00"))

## --- SENSITIVITY STUDY: MOSQUITO SAMPLES --- ##
culex_SENS.reads <- read_xlsx("removal_reads_sens_M.xlsx")
head(culex_SENS.reads)
culex_SENS.reads_long <- gather(culex_SENS.reads, Reads, number , Initial_reads:Classified_reads, factor_key=TRUE)
culex_SENS.reads_long$Reads <- str_replace_all(culex_SENS.reads_long$Reads, "Initial_reads","Initial number of reads" ) 
culex_SENS.reads_long$Reads <- str_replace_all(culex_SENS.reads_long$Reads, "after_trim", paste0( "Number of reads after removal", "\n",  "of duplications and low quality reads"))
culex_SENS.reads_long$Reads <- str_replace_all(culex_SENS.reads_long$Reads, "after_mapping","Number of reads after removal of host genome" ) 
culex_SENS.reads_long$Reads <- str_replace_all(culex_SENS.reads_long$Reads, "Classified_reads","Number of reads classified as taxa" ) 


## prepare labels
molecules <-c(3000000, 1.5e+06, 300000, 30000, 3000)
my.labels <- paste0(culex_SENS.reads_long$Sample, "\n", format(molecules, scientific = TRUE), " molecules")

ggplot(culex_SENS.reads_long, aes(x = as.factor(Sample), y = number, fill=Reads, col=Reads)) + geom_col(position = "dodge") + theme_bw() + scale_y_break(c(4533209, 25562418), scales=0.2) +
  scale_y_continuous(name ="Number of reads", breaks=c(0,200000,500000,1e+06, 2e+06,3e+06, 4e+06, 3e+07, 5e+07, 6e+07), position = "left") +xlab("Samples") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=11), axis.title = element_text(size=11), legend.text = element_text(size=11), axis.text.y.right  = element_blank(), axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="")  + guides(col=FALSE, fill= guide_legend(ncol = 2, nrow = 2, byrow = TRUE)) +
  scale_color_manual(values =c("skyblue4", "darkgoldenrod", "forestgreen", "tomato3")) + scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#D55E00")) + 
  scale_x_discrete(labels=my.labels)


## --- SENSITIVITY STUDY: MOSQUITO SAMPLES --- ##
fish_SENS.reads <- read_xlsx("removal_reads_Fish_sens.xlsx")
head(fish_SENS.reads)
fish_SENS.reads_long <- gather(fish_SENS.reads, Reads, number , Initial_reads:Classified_reads, factor_key=TRUE)
fish_SENS.reads_long$Reads <- str_replace_all(fish_SENS.reads_long$Reads, "Initial_reads","Initial number of reads" ) 
fish_SENS.reads_long$Reads <- str_replace_all(fish_SENS.reads_long$Reads, "after_trim", paste0( "Number of reads after removal", "\n",  "of duplications and low quality reads"))
fish_SENS.reads_long$Reads <- str_replace_all(fish_SENS.reads_long$Reads, "after_mapping","Number of reads after removal of host genome" ) 
fish_SENS.reads_long$Reads <- str_replace_all(fish_SENS.reads_long$Reads, "Classified_reads","Number of reads classified as taxa" ) 


## prepare labels
molecules <-c(3000000, 1.5e+06, 300000, 30000, 3000)
my.labels2 <- paste0(fish_SENS.reads_long$Sample, "\n", format(molecules, scientific = TRUE), " molecules")

ggplot(fish_SENS.reads_long, aes(x = as.factor(Sample), y = number, fill=Reads, col=Reads)) + geom_col(position = "dodge") + theme_bw() + scale_y_break(c(4533209, 25562418), scales=0.2) +
  scale_y_continuous(name ="Number of reads", breaks=c(0,200000,500000,1e+06, 2e+06,3e+06, 4e+06, 3e+07,6e+07, 9e+07), position = "left") +xlab("Samples") +
  theme(legend.title = element_text(), legend.position = "top", axis.text =element_text(size=11), axis.title = element_text(size=11), legend.text = element_text(size=11), axis.text.y.right  = element_blank(), axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill="")  + guides(col=FALSE, fill= guide_legend(ncol = 2, nrow = 2, byrow = TRUE)) +
  scale_color_manual(values =c("skyblue4", "darkgoldenrod", "forestgreen", "tomato3")) + scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#D55E00")) + 
  scale_x_discrete(labels=my.labels2)






## R version 4.0.4 (2022-05-31)
##
## Analysis of kaiju results 
##
##
##  Copyright: Marta Ibañez Lligoña (marta.ibanez@vhir.org)
##

## --- Set working directory --- ##
setwd("~/Desktop/Master_VHIR/TFM/Metagenomics_fastq2")

## Load packages
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(edgeR)
library(ggbreak)

### --- MOSQUITO --- ##
## --- Load data --- ##
## --- SAMPLE: MA --- ##
MA <- read_xlsx("MA/kaiju.taxonpaths_MA.xlsx")

## --- SAMPLE: MB --- ##
MB <- read_xlsx("MB/kaiju.taxonpaths_MB.xlsx")

## --- SAMPLE: MC --- ##
MC <- read_xlsx("MC/kaiju.taxonpaths_MC.xlsx")

## --- SAMPLE: MD --- ##
MD <- read_xlsx("MD/kaiju.taxonpaths_MD.xlsx")

## --- SAMPLE: ME --- ##
ME <- read_xlsx("ME/kaiju.taxonpaths_ME.xlsx")

## --- Keep HCV classification only --- ##
## Selecting rows that have been detectet as hepacivirus
## --- MA --- ##
MA.HCV <- MA %>% filter(Species == "Hepacivirus C")
MA.reads_HCV <- sum(MA.HCV$`Number of reads`)
MA.all <- sum(MA$`Number of reads`)
MA.classified <- MA %>% filter(Root == "root")
MA.classified <- sum(MA.classified$`Number of reads`)

## --- MB --- ##
MB.HCV <- MB %>% filter(Species == "Hepacivirus C")
MB.reads_HCV <- sum(MB.HCV$`Number of reads`)
MB.all <- sum(MB$`Number of reads`)
MB.classified <- MB %>% filter(Root == "root")
MB.classified <- sum(MB.classified$`Number of reads`)

## --- MC --- ##
MC.HCV <- MC %>% filter(Species == "Hepacivirus C")
MC.reads_HCV <- sum(MC.HCV$`Number of reads`)
MC.all <- sum(MC$`Number of reads`)
MC.classified <- MC %>% filter(Root == "root")
MC.classified <- sum(MC.classified$`Number of reads`)

## --- MD --- ##
MD.HCV <- MD %>% filter(Species == "Hepacivirus C")
MD.reads_HCV <- sum(MD.HCV$`Number of reads`)
MD.all <- sum(MD$`Number of reads`)
MD.classified <- MD %>% filter(Root == "root")
MD.classified <- sum(MD.classified$`Number of reads`)

## --- ME --- ##
ME.HCV <- ME %>% filter(Species == "Hepacivirus C")
ME.reads_HCV <- sum(ME.HCV$`Number of reads`)
ME.all <- sum(ME$`Number of reads`)
ME.classified <- ME %>% filter(Root == "root")
ME.classified <- sum(ME.classified$`Number of reads`)

## --- Create table with values --- ##
M.values <- t(data.frame(MA.all, MB.all, MC.all, MD.all, ME.all))
M.classified <- t(data.frame(MA.classified, MB.classified, MC.classified, MD.classified, ME.classified))
M.classified_per <- data.frame(round(M.classified*100/M.values))
M.HCV <- t(data.frame(MA.reads_HCV, MB.reads_HCV, MC.reads_HCV, MD.reads_HCV, ME.reads_HCV))
M.HCV_per <- data.frame((M.HCV*100/M.values))
M.QA <- data.frame(M.values, M.classified, M.HCV)
M.QA$sample <- c("MA", "MB", "MC", "MD", "ME")
colnames(M.QA) <- c("All_reads", "Classified_reads", "HCV", "Sample")

## --- Save table --- ##
write.csv2(M.QA, "M.QA.csv")

## --- Transform table to long format --- ##
M.QA_long <- gather(M.QA, condition, reads , All_reads:HCV, factor_key=TRUE)

## --- Plot --- ##
ggplot(M.QA_long, aes(x = Sample, y = reads, group=condition, color=condition)) + geom_line() + geom_point() + theme_bw() +
  facet_grid(condition ~., scales="free_y") + ylab("Number of reads")

### --- FISH --- ##
## --- Load data --- ##
## --- SAMPLE: PA --- ##
PA <- read_xlsx("PA/kaiju.taxonpaths_PA.xlsx")

## --- SAMPLE: PB --- ##
PB <- read_xlsx("PB/kaiju.taxonpaths_PB.xlsx")

## --- SAMPLE: PC --- ##
PC <- read_xlsx("PC/kaiju.taxonpaths_PC.xlsx")

## --- SAMPLE: PD --- ##
PD <- read_xlsx("PD/kaiju.taxonpaths_PD.xlsx")

## --- SAMPLE: PE --- ##
PE <- read_xlsx("PE/kaiju.taxonpaths_PE.xlsx")

## --- Keep HCV classification only --- ##
## Selecting rows that have been detectet as hepacivirus
## --- PA --- ##
PA.HCV <- PA %>% filter(Species == "Hepacivirus C")
PA.reads_HCV <- sum(PA.HCV$`Number of reads`)
PA.all <- sum(PA$`Number of reads`)
PA.classified <- PA %>% filter(Root == "root")
PA.classified <- sum(PA.classified$`Number of reads`)
## --- PB --- ##
PB.HCV <- PB %>% filter(Species == "Hepacivirus C")
PB.reads_HCV <- sum(PB.HCV$`Number of reads`)
PB.all <- sum(PB$`Number of reads`)
PB.classified <- PB %>% filter(Root == "root")
PB.classified <- sum(PB.classified$`Number of reads`)

## --- PC --- ##
PC.HCV <- PC %>% filter(Species == "Hepacivirus C")
PC.reads_HCV <- sum(PC.HCV$`Number of reads`)
PC.all <- sum(PC$`Number of reads`)
PC.classified <- PC %>% filter(Root == "root")
PC.classified <- sum(PC.classified$`Number of reads`)

## --- PD --- ##
PD.HCV <- PD %>% filter(Species == "Hepacivirus C")
PD.reads_HCV <- sum(PD.HCV$`Number of reads`)
PD.all <- sum(PD$`Number of reads`)
PD.classified <- PD %>% filter(Root == "root")
PD.classified <- sum(PD.classified$`Number of reads`)

## --- PE --- ##
PE.HCV <- PE %>% filter(Species == "Hepacivirus C")
PE.reads_HCV <- sum(PE.HCV$`Number of reads`)
PE.all <- sum(PE$`Number of reads`)
PE.classified <- PE %>% filter(Root == "root")
PE.classified <- sum(PE.classified$`Number of reads`)

## --- Create table with values --- ##
P.values <- t(data.frame(PA.all, PB.all, PC.all, PD.all, PE.all))
P.classified <- t(data.frame(PA.classified, PB.classified, PC.classified, PD.classified, PE.classified))
P.HCV <- t(data.frame(PA.reads_HCV, PB.reads_HCV, PC.reads_HCV, PD.reads_HCV, PE.reads_HCV))
P.QA <- data.frame(P.values, P.classified, P.HCV)
P.QA$sample <- c("PA", "PB", "PC", "PD", "PE")
colnames(P.QA) <- c("All_reads", "Classified_reads", "HCV", "Sample")

## --- Save table --- ##
write.csv2(P.QA, "P.QA.csv")

## --- Transform table to long format --- ##
P.QA_long <- gather(P.QA, condition, reads , All_reads:HCV, factor_key=TRUE)

## --- Plot --- ##
ggplot(P.QA_long, aes(x = Sample, y = reads, group=condition, color=condition)) + geom_line() + geom_point() + theme_bw() +
  facet_grid(condition ~., scales="free_y") + ylab("Number of reads")

## --- CPM: IN-DEPTH ANALYSIS --- ##
## --- Mosquito --- ##
cpm.M.HCV <- (t(M.HCV)/t(M.values))* 10^6
#cpm.M.class <- (t(M.classified)/t(M.values))* 10^6
colnames(cpm.M.HCV) <- c("MA", "MB", "MC", "MD", "ME")
#colnames(cpm.M.class) <- c("MA", "MB", "MC", "MD", "ME")
M.valuescpm <- (t(M.values)/t(M.values))* 10^6
colnames(M.valuescpm) <- c("MA", "MB", "MC", "MD", "ME")
cpm.all <- rbind(M.valuescpm, cpm.M.HCV)
cpm.all_long <- t(cpm.all)
colnames(cpm.all_long) <- c("Total", "HCV")
cpm.all_long <- as.data.frame(cpm.all_long)
cpm.all_long$sample <- rownames(cpm.all_long)

## --- Transform table to long format --- ##
cpm.all_long2 <- gather(cpm.all_long, condition, reads , Total:HCV, factor_key=TRUE)
cpm.reads <- cpm.all_long2[cpm.all_long2$condition =="HCV",]


## --- Plot: counts per million --- ##
raw_class <- read_xlsx("raw_classification.xlsx")
raw_class.culex <- raw_class[raw_class$Tipus == "Mosquit",]
raw_class.culex$cpm <- (raw_class.culex$`HCV reads`/raw_class.culex$`Nº reads inicial`)* 10^6

## --- With total number of reads --- ##
ggplot(cpm.all_long2, aes(x = sample, y = reads, group=condition, color=condition)) + geom_line() + geom_point() + theme_bw() +
  facet_grid(condition ~., scales="free_y") + ylab("Counts per million")

## --- With raw number of reads classified as HCV --- ##
cpm.reads_mosquito <- data.frame(cpm.reads$sample, cpm.reads$reads)
colnames(cpm.reads_mosquito) <- c("sample", "cpm")
cpm.reads_mosquito$group <- "Filtered reads"
cpm.reads_mosquito2 <- data.frame(cpm.reads$sample, raw_class.culex$cpm)
colnames(cpm.reads_mosquito2) <- c("sample", "cpm")
cpm.reads_mosquito2$group <- "Raw reads"
cpm.reads_mosquito <- rbind(cpm.reads_mosquito, cpm.reads_mosquito2)
cpm.reads_mosquito$number_molecules<- c(3000000, 1500000, 300000, 30000, 3000, 3000000, 1500000, 300000, 30000, 3000)
molecules <-c(3000000, 1.5e+06, 300000, 30000, 3000)
my.labels <- paste0(cpm.reads_mosquito2$sample, "\n", format(molecules, scientific = TRUE), " molecules")
## plot
mosquito <- ggplot(cpm.reads_mosquito, aes(x = sample, y = cpm, group=group, color=group)) + geom_line() + geom_point() + theme_bw() + ylab("Counts per million") + xlab("Samples")
mosquito + labs(color='Type of reads') + theme(legend.position = "top" , axis.text =element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), legend.title = element_text(size=12)) + scale_x_discrete(labels=my.labels) + scale_fill_brewer()

## --- FISH --- ##
cpm.P.HCV <- (t(P.HCV)/t(P.values))* 10^6
colnames(cpm.M.HCV) <- c("PA", "PB", "PC", "PD", "PE")
P.valuescpm <- (t(P.values)/t(P.values))* 10^6
colnames(P.valuescpm) <- c("PA", "PB", "PC", "PD", "PE")
Pcpm.all <- rbind(P.valuescpm, cpm.P.HCV)
Pcpm.all_long <- t(Pcpm.all)
colnames(Pcpm.all_long) <- c("All", "HCV")
Pcpm.all_long <- as.data.frame(Pcpm.all_long)
Pcpm.all_long$sample <- rownames(Pcpm.all_long)

## --- Transform table to long format --- ##
Pcpm.all_long2 <- gather(Pcpm.all_long, condition, reads , All:HCV, factor_key=TRUE)

## --- Plot: counts per million --- ##

ggplot(Pcpm.all_long2, aes(x = sample, y = reads, group=condition, color=condition)) + geom_line() + geom_point() + theme_bw() +
  facet_grid(condition ~., scales="free_y") + ylab("Counts per million")

## --- With raw number of reads classified as HCV --- ##
cpm.reads_fish <- data.frame(Pcpm.all_long$sample, t(cpm.P.HCV))
colnames(cpm.reads_fish) <- c("sample", "cpm")
cpm.reads_fish$group <- "Filtered reads"
## select G. affinis HCV reads
raw_class.fish <- raw_class[raw_class$Tipus == "Peix",]
raw_class.fish$cpm <- (raw_class.fish$`HCV reads`/raw_class.fish$`Nº reads inicial`)* 10^6
## create data framw with cpm values
cpm.reads_fish2 <- data.frame(cpm.reads_fish$sample, raw_class.fish$cpm)
colnames(cpm.reads_fish2) <- c("sample", "cpm")
cpm.reads_fish2$group <- "Raw reads"
cpm.reads_fish <- rbind(cpm.reads_fish, cpm.reads_fish2)
molecules <-c(3000000, 1.5e+06, 300000, 30000, 3000)
my.labels <- paste0(cpm.reads_fish2$sample, "\n", format(molecules, scientific = TRUE), " molecules")
## plot
fish <- ggplot(cpm.reads_fish, aes(x = sample, y = cpm, group=group, color=group)) + geom_line() + geom_point() + theme_bw() + ylab("Counts per million") + xlab("Samples")
fish + labs(color='Type of reads') + theme(legend.position = "top", axis.text.x=element_text(size=12), axis.text =element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), legend.title = element_text(size=12)) + scale_x_discrete(labels=my.labels) + scale_fill_brewer()


theme(axis.text.x=element_text(size=15))


### --- MOSQUITO --- ###
## --- Correlation: Number of reads vs molecules --- ##
cpm.all_molecules <- cpm.all 
cpm.all_molecules  <- as.data.frame(cpm.all_molecules)
cpm.all_molecules[3,] <- c(3000000, 1500000, 300000, 30000, 3000)
cpm.all_molecules <- t(cpm.all_molecules)
colnames(cpm.all_molecules) <- c("Reads", "HCV", "Molecules")
cpm.all_molecules <- as.data.frame(cpm.all_molecules)
cpm.all_molecules$Sample <- rownames(cpm.all_molecules)


## --- Transform table to long format --- ##
cpm.all_molecules2 <- gather(as.data.frame(cpm.all_molecules), condition, reads , Reads:Molecules, factor_key=TRUE)

## --- PLOT --- ##
total_counts_single <- ggplot(cpm.all_molecules,
                              aes(x = HCV , y = Molecules)) +
  geom_point(aes(col = as.factor(Sample))) +
  #geom_smooth(method = "lm") +
  labs(x = "HCV reads cpm",
       y = "Number of molecules",
       title = "Effect of read depth on HCV") + theme_bw()

cor.test(cpm.all_molecules$HCV, cpm.all_molecules$Molecules)

cpm.all_molecules_cv <- cpm.all_molecules
cpm.all_molecules_cv$raw_reads <- M.HCV


## Conversion reads to molecules 
cpm.all_molecules_cv$lost_per <- (cpm.all_molecules_cv$Molecules-cpm.all_molecules_cv$raw_reads)/cpm.all_molecules_cv$Molecules

## plot
ggplot(data = cpm.all_molecules_cv, aes(x=Sample, y=conv)) + geom_point()


### --- FISH --- ##
Pcpm.all_molecules <- Pcpm.all 
Pcpm.all_molecules  <- as.data.frame(Pcpm.all_molecules)
Pcpm.all_molecules[3,] <- c(3000000, 1500000, 300000, 30000, 3000)
Pcpm.all_molecules <- t(Pcpm.all_molecules)
colnames(Pcpm.all_molecules) <- c("Reads", "HCV", "Molecules")
Pcpm.all_molecules <- as.data.frame(Pcpm.all_molecules)
Pcpm.all_molecules$Sample <- rownames(Pcpm.all_molecules)
Pcpm.all_molecules$raw_reads <- P.HCV
Pcpm.all_molecules$lost_per <- (Pcpm.all_molecules$Molecules-Pcpm.all_molecules$raw_reads)/Pcpm.all_molecules$Molecules





## ---- METHODS: SUMMARY PLOT (clean reads) ---- ##
summary_M <- read_xlsx("TFM_reads.xlsx") 
colnames(summary_M) <- c("Type", "Sample", "HCV_mol", "Initial_reads", "After_quality_trim", "After_host_genome_mapping", "Classified_reads", "Per_classified")

## --- Transform table to long format --- ##
summary_M.long <- gather(summary_M, condition, reads , Initial_reads:Classified_reads, factor_key=TRUE)

## --- Plot: Evolution of reads --- ##
m <-ggplot(summary_M.long, aes(x = Sample, y = reads, fill=condition)) + geom_col(position = "dodge") +  theme_bw() +
  theme(legend.position = "top") + ylab("Number of reads")  + guides(fill=guide_legend(title=element_blank())) + 
  scale_fill_brewer(name = "condition", labels = c("Initial reads", "Reads after quality trim", "Reads after host genome mapping", "Classified reads"),palette = "Dark2") 



## ---- KAIJU RESULTS FROM NON-PROCESSED DATA ---- ##
## --- MOSQUITO --- ##
detected_hepacivirus.M <- data.frame(753,15,23,5,21)
colnames(detected_hepacivirus.M) <- c("MA", "MB","MC", "MD", "ME")
detected_hepacivirus.M <- t(detected_hepacivirus.M)
colnames(detected_hepacivirus.M) <- "HCV_reads"
detected_hepacivirus.M <- as.data.frame(detected_hepacivirus.M)
detected_hepacivirus.M$all_reads <- M.QA$All_reads

## --- Add columns to data frame --- ##
detected_hepacivirus.M$per_hcv <- (detected_hepacivirus.M$HCV_reads/detected_hepacivirus.M$all_reads)*100
detected_hepacivirus.M$cpm <- (detected_hepacivirus.M$HCV_reads/detected_hepacivirus.M$all_reads)* 10^6
detected_hepacivirus.M$sample <-  rownames(detected_hepacivirus.M)

## --- plot --- ##
detected_hepacivirus.M_long <- gather(detected_hepacivirus.M, condition, reads , HCV_reads:cpm, factor_key=TRUE)
ggplot(detected_hepacivirus.M_long, aes(x = sample, y = reads, group=condition, color=condition)) + geom_line() + geom_point() + theme_bw() + ylab("Counts per million") + 
  geom_line() + facet_grid(condition ~., scales="free_y")

## --- FISH --- ##
detected_hepacivirus.P <- data.frame(440,21,12,14,0)
colnames(detected_hepacivirus.P) <- c("PA", "PB","PC", "PD", "PE")
detected_hepacivirus.P <- t(detected_hepacivirus.P)
colnames(detected_hepacivirus.P) <- "HCV_reads"
detected_hepacivirus.P <- as.data.frame(detected_hepacivirus.P)
detected_hepacivirus.P$all_reads <- P.QA$All_reads

## --- Add columns to data frame --- ##
detected_hepacivirus.P$per_hcv <- (detected_hepacivirus.P$HCV_reads/detected_hepacivirus.P$all_reads)*100
detected_hepacivirus.P$cpm <- (detected_hepacivirus.P$HCV_reads/detected_hepacivirus.P$all_reads)* 10^6
detected_hepacivirus.P$sample <-  rownames(detected_hepacivirus.P)

## --- plot --- ##
detected_hepacivirus.P_long <- gather(detected_hepacivirus.P, condition, reads , HCV_reads:cpm, factor_key=TRUE)
ggplot(detected_hepacivirus.P_long, aes(x = sample, y = reads, group=condition, color=condition)) + geom_line() + geom_point() + theme_bw() + ylab("Counts per million") + 
  geom_line() + facet_grid(condition ~., scales="free_y")


## --- JOIN NO-FILTERED AND FILTERED DATA -- ##
## --- G. affinis --- ##
## --- Plot: all reads, filtered vs non-filtered -- ##
all.reads.P <- data.frame(detected_hepacivirus.P$HCV_reads, detected_hepacivirus.P$all_reads)
colnames(all.reads.P) <- c("No_filter_HCV.reads", "Filtered.all_reads" )

all.reads.P$Filtered.HCV_reads <- P.QA$HCV
all.reads.P$No_filter.all_reads <- c(63396643,
                                     42902068,
                                     100826653,
                                     36098688,
                                     18796935)
all.reads.P$sample <- c("PA", "PB","PC", "PD", "PE")
## Transform to long format
all.reads.P_long <- gather(all.reads.P, condition, reads , No_filter_HCV.reads:No_filter.all_reads, factor_key=TRUE)

all.reads.P_long$color <- "Filtered"
all.reads.P_long$color[grep("No", all.reads.P_long$condition)] <- "Non-filtered"
all.reads.P_long$data <- "HCV"
all.reads.P_long$data[grep("all", all.reads.P_long$condition)] <- "All"

ggplot(all.reads.P_long, aes(x = sample, y = reads, group=condition, color=color, shape=data)) + geom_line() + geom_point() + theme_bw() + ylab("Number of reads") + xlab("Samples") +
  geom_line()  ##facet_grid(data ~., scales="free_y") +
  ##scale_fill_brewer(name = "condition", labels = c("Initial reads", "Reads after quality trim", "Reads after host genome mapping", "Classified reads"),palette = "Dark2") 
  
## --- Only hcv reads cpm --- ##
## Calculate normalization and percentages of HCV reads for both filtered and unfiltered data
all.reads.P$cpm_HCV.no_filtered <- (all.reads.P$No_filter_HCV.reads/all.reads.P$No_filter.all_reads)* 10^6
all.reads.P$per_HCV.no_filtered <- (all.reads.P$No_filter_HCV.reads/all.reads.P$No_filter.all_reads)*100
all.reads.P$cpm_HCV.filtered <- (all.reads.P$Filtered.HCV_reads/all.reads.P$Filtered.all_reads)* 10^6
all.reads.P$per_HCV.filtered <- (all.reads.P$Filtered.HCV_reads/all.reads.P$Filtered.all_reads)*100

## Transform to long format
all.reads.P_long2 <- gather(all.reads.P, condition, reads , No_filter_HCV.reads:No_filter.all_reads, factor_key=TRUE)

## --- PLOT: cpm values for HCV (filtered and unfiltered) --- ##
ggplot(all.reads.P_long, aes(x = sample, y = reads, group=condition, color=color, shape=data)) + geom_line() + geom_point() + theme_bw() + ylab("Number of reads") + xlab("Samples") +
  geom_line()

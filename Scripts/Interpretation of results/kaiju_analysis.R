#########################################
###### Kaiju results interpretation #####
#### Copyright: Marta Ibañez Lligoña ####
#########################################

## --- Set working directory --- #
setwd("~/Desktop/Master_VHIR/TFM/Metagenomics_fastq1/kaiju_analysis_processed/")

## --- Load libraries, if necessary --- ##
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggvenn)
library(plyr)
library(ggrepel)

## --- Culex sp. --- ##
data1_L1 <- read_excel("kaiju.taxonpaths_1_L1.xlsx")
data2_L1 <- read_excel("kaiju.taxonpaths_2_L1.xlsx")
data3_L1 <- read_excel("kaiju.taxonpaths_3_L1.xlsx")
data4_L1 <- read_excel("kaiju.taxonpaths_4_L1.xlsx")
data5_L1 <- read_excel("kaiju.taxonpaths_5_L1.xlsx")
data6_L1 <- read_excel("kaiju.taxonpaths_6_L1.xlsx")
class(data1_L1$`Number of reads`)

## --- Sort rows by number of reads assigned to each taxa --- ##
data1_L1 <- data1_L1[order(data1_L1$`Number of reads`, decreasing = TRUE),]
data2_L1 <- data2_L1[order(data2_L1$`Number of reads`, decreasing = TRUE),]
data3_L1 <- data3_L1[order(data3_L1$`Number of reads`, decreasing = TRUE),]
data4_L1 <- data4_L1[order(data4_L1$`Number of reads`, decreasing = TRUE),]
data5_L1 <- data5_L1[order(data5_L1$`Number of reads`, decreasing = TRUE),]
data6_L1 <- data6_L1[order(data6_L1$`Number of reads`, decreasing = TRUE),]

## --- Domains --- ##
## environmental samples --> taxa that cannot be cultured under laboratory methods
## --- SAMPLE 1 --- ##
unique(data1_L1$Domain)
data1_L1[1,4] <- "Unclassified"
data1_L1.longD <- data1_L1[,c(1,4)]
data1_L1.longD <- aggregate(.~Domain,data=data1_L1.longD,FUN=sum)

## --- SAMPLE 2 --- ##
unique(data2_L1$Domain)
data2_L1[1,4] <- "Unclassified"
data2_L1.longD <- data2_L1[,c(1,4)]
data2_L1.longD <- aggregate(.~Domain,data=data2_L1.longD,FUN=sum)

## --- SAMPLE 3 --- ##
unique(data3_L1$Domain)
data3_L1[1,4] <- "Unclassified"
data3_L1.longD <- data3_L1[,c(1,4)]
data3_L1.longD <- aggregate(.~Domain,data=data3_L1.longD,FUN=sum)

## --- SAMPLE 4 --- ##
unique(data4_L1$Domain)
data4_L1[1,4] <- "Unclassified"
data4_L1.longD <- data4_L1[,c(1,4)]
data4_L1.longD <- aggregate(.~Domain,data=data4_L1.longD,FUN=sum)

## --- SAMPLE 5 --- ##
unique(data5_L1$Domain)
data5_L1[1,4] <- "Unclassified"
data5_L1.longD <- data5_L1[,c(1,4)]
data5_L1.longD <- aggregate(.~Domain,data=data5_L1.longD,FUN=sum)

## --- SAMPLE 6 --- ##
unique(data6_L1$Domain)
data6_L1[1,4] <- "Unclassified"
data6_L1.longD <- data6_L1[,c(1,4)]
data6_L1.longD <- aggregate(.~Domain,data=data6_L1.longD,FUN=sum)


## --- VIRUS --- ##

## --- Select viruses from table --- ##
## --- SAMPLE 1 --- ##
all_viruses1 <- data1_L1[data1_L1$Type == "Viruses",]
all_viruses1 <- all_viruses1[is.na(all_viruses1$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses1$Family)
vir1.fam <- data.frame(all_viruses1$`Number of reads`, all_viruses1$Family)
colnames(vir1.fam) <- c("Reads", "Family")

## NA values can be specified with a higher classification
vir1.fam_long <- aggregate(.~Family,data=vir1.fam,FUN=sum)


## --- SAMPLE 2 --- ##
all_viruses2 <- data2_L1[data2_L1$Type == "Viruses",]
all_viruses2 <- all_viruses2[is.na(all_viruses2$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses2$Family)
vir2.fam <- data.frame(all_viruses2$`Number of reads`, all_viruses2$Family)
colnames(vir2.fam) <- c("Reads", "Family")
## NA values can be specified with a higher classification
vir2.fam_long <- aggregate(.~Family,data=vir2.fam,FUN=sum)

## --- SAMPLE 3 --- ##
all_viruses3 <- data3_L1[data3_L1$Type == "Viruses",]
all_viruses3 <- all_viruses3[is.na(all_viruses3$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses3$Family)
vir3.fam <- data.frame(all_viruses3$`Number of reads`, all_viruses3$Family)
colnames(vir3.fam) <- c("Reads", "Family")
## NA values can be specified with a higher classification
vir3.fam_long <- aggregate(.~Family,data=vir3.fam,FUN=sum)

## --- SAMPLE 4 --- ##
all_viruses4 <- data4_L1[data4_L1$Type == "Viruses",]
all_viruses4 <- all_viruses4[is.na(all_viruses4$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses4$Family)
vir4.fam <- data.frame(all_viruses4$`Number of reads`, all_viruses4$Family)
colnames(vir4.fam) <- c("Reads", "Family")
## NA values can be specified with a higher classification
vir4.fam_long <- aggregate(.~Family,data=vir4.fam,FUN=sum)

## --- SAMPLE 5 --- ##
all_viruses5 <- data5_L1[data5_L1$Type == "Viruses",]
all_viruses5 <- all_viruses5[is.na(all_viruses5$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses5$Family)
vir5.fam <- data.frame(all_viruses5$`Number of reads`, all_viruses5$Family)
colnames(vir5.fam) <- c("Reads", "Family")
## NA values can be specified with a higher classification
vir5.fam_long <- aggregate(.~Family,data=vir5.fam,FUN=sum)

## --- SAMPLE 6 --- ##
all_viruses6 <- data6_L1[data6_L1$Type == "Viruses",]
all_viruses6 <- all_viruses6[is.na(all_viruses6$`Number of reads`)==FALSE,]
unique(all_viruses6$Family)
vir6.fam <- data.frame(all_viruses6$`Number of reads`, all_viruses6$Family)
colnames(vir6.fam) <- c("Reads", "Family")
## NA values can be specified with a higher classification
vir6.fam_long <- aggregate(.~Family,data=vir6.fam,FUN=sum)

## --- G. affinis --- ##

## Load data
data7_L1 <- read_excel("kaiju.taxonpaths_7_L1.xlsx")
data8_L1 <- read_excel("kaiju.taxonpaths_8_L1.xlsx")
data9_L1 <- read_excel("kaiju.taxonpaths_9_L1.xlsx")
data10_L1 <- read_excel("kaiju.taxonpaths_10_L1.xlsx")
data11_L1 <- read_excel("kaiju.taxonpaths_11_L1.xlsx")
data12_L1 <- read_excel("kaiju.taxonpaths_12_L1.xlsx")

## --- Sort rows by number of reads assigned to each taxa --- ##
data7_L1 <- data7_L1[order(data7_L1$`Number of reads`, decreasing = TRUE),]
data8_L1 <- data8_L1[order(data8_L1$`Number of reads`, decreasing = TRUE),]
data9_L1 <- data9_L1[order(data9_L1$`Number of reads`, decreasing = TRUE),]
data10_L1 <- data10_L1[order(data10_L1$`Number of reads`, decreasing = TRUE),]
data11_L1 <- data11_L1[order(data11_L1$`Number of reads`, decreasing = TRUE),]
data12_L1 <- data12_L1[order(data12_L1$`Number of reads`, decreasing = TRUE),]

## --- Domains --- ##
## environmental samples --> taxa that cannot be cultured under laboratory methods
## --- SAMPLE 7 --- ##
unique(data7_L1$Domain)
data7_L1[1,4] <- "Unclassified"
data7_L1.longD <- data7_L1[,c(1,4)]
data7_L1.longD <- aggregate(.~Domain,data=data7_L1.longD,FUN=sum)

## --- SAMPLE 8 --- ##
unique(data8_L1$Domain)
data8_L1[1,4] <- "Unclassified"
data8_L1.longD <- data8_L1[,c(1,4)]
data8_L1.longD <- aggregate(.~Domain,data=data8_L1.longD,FUN=sum)

## --- SAMPLE 9 --- ##
unique(data9_L1$Domain)
data9_L1[1,4] <- "Unclassified"
data9_L1.longD <- data9_L1[,c(1,4)]
data9_L1.longD <- aggregate(.~Domain,data=data9_L1.longD,FUN=sum)

## --- SAMPLE 10 --- ##
unique(data10_L1$Domain)
data10_L1[1,4] <- "Unclassified"
data10_L1.longD <- data10_L1[,c(1,4)]
data10_L1.longD <- aggregate(.~Domain,data=data10_L1.longD,FUN=sum)

## --- SAMPLE 11 --- ##
unique(data11_L1$Domain)
data11_L1[1,4] <- "Unclassified"
data11_L1.longD <- data11_L1[,c(1,4)]
data11_L1.longD <- aggregate(.~Domain,data=data11_L1.longD,FUN=sum)

## --- SAMPLE 12 --- ##
unique(data12_L1$Domain)
data12_L1[1,4] <- "Unclassified"
data12_L1.longD <- data12_L1[,c(1,4)]
data12_L1.longD <- aggregate(.~Domain,data=data12_L1.longD,FUN=sum)


## --- VIRUS --- ##

## --- Select viruses from table --- ##
## --- SAMPLE 7 --- ##
all_viruses7 <- data7_L1[data7_L1$Type == "Viruses",]
all_viruses7 <- all_viruses7[is.na(all_viruses7$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses7$Family)
vir7.fam <- data.frame(all_viruses7$`Number of reads`, all_viruses7$Family)
colnames(vir7.fam) <- c("Reads", "Family")
vir7.fam_long <- aggregate(.~Family,data=vir7.fam,FUN=sum)

## --- SAMPLE 8 --- ##
all_viruses8 <- data8_L1[data8_L1$Type == "Viruses",]
all_viruses8 <- all_viruses8[is.na(all_viruses8$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses8$Family)
vir8.fam <- data.frame(all_viruses8$`Number of reads`, all_viruses8$Family)
colnames(vir8.fam) <- c("Reads", "Family")
vir8.fam_long <- aggregate(.~Family,data=vir8.fam,FUN=sum)

## --- SAMPLE 9 --- ##
all_viruses9 <- data9_L1[data9_L1$Type == "Viruses",]
all_viruses9 <- all_viruses9[is.na(all_viruses9$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses9$Family)
vir9.fam <- data.frame(all_viruses9$`Number of reads`, all_viruses9$Family)
colnames(vir9.fam) <- c("Reads", "Family")
vir9.fam_long <- aggregate(.~Family,data=vir9.fam,FUN=sum)

## --- SAMPLE 10 --- ##
all_viruses10 <- data10_L1[data10_L1$Type == "Viruses",]
all_viruses10 <- all_viruses10[is.na(all_viruses10$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses10$Family)
vir10.fam <- data.frame(all_viruses10$`Number of reads`, all_viruses10$Family)
colnames(vir10.fam) <- c("Reads", "Family")
vir10.fam_long <- aggregate(.~Family,data=vir10.fam,FUN=sum)

## --- SAMPLE 11 --- ##
all_viruses11 <- data11_L1[data11_L1$Type == "Viruses",]
all_viruses11 <- all_viruses11[is.na(all_viruses11$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses11$Family)
vir11.fam <- data.frame(all_viruses11$`Number of reads`, all_viruses11$Family)
colnames(vir11.fam) <- c("Reads", "Family")
vir11.fam_long <- aggregate(.~Family,data=vir11.fam,FUN=sum)

## --- SAMPLE 12 --- ##
all_viruses12 <- data12_L1[data12_L1$Type == "Viruses",]
all_viruses12 <- all_viruses12[is.na(all_viruses12$`Number of reads`)==FALSE,]
## Which families are present?
unique(all_viruses12$Family)
vir12.fam <- data.frame(all_viruses12$`Number of reads`, all_viruses12$Family)
colnames(vir12.fam) <- c("Reads", "Family")
vir12.fam_long <- aggregate(.~Family,data=vir12.fam,FUN=sum)

## --- Save tables to files --- ##
write.csv(all_viruses12, "virome_12.csv")


## --- Culex sp: Domains --- ##

## --- SAMPLE 1 --- ##
data1_L1.longD$sample <- 1
data1_L1.longD$cpm <- (data1_L1.longD$`Number of reads`/sum(data1_L1.longD$`Number of reads`))* 10^6
data1_L1.longD$per <- (data1_L1.longD$`Number of reads`/sum(data1_L1.longD$`Number of reads`))* 100

## --- SAMPLE 2 --- ##
data2_L1.longD$sample <- 2
data2_L1.longD$cpm <- (data2_L1.longD$`Number of reads`/sum(data2_L1.longD$`Number of reads`))* 10^6
data2_L1.longD$per <- (data2_L1.longD$`Number of reads`/sum(data2_L1.longD$`Number of reads`))* 100

## --- SAMPLE 3 --- ##
data3_L1.longD$sample <- 3
data3_L1.longD$cpm <- (data3_L1.longD$`Number of reads`/sum(data3_L1.longD$`Number of reads`))* 10^6
data3_L1.longD$per <- (data3_L1.longD$`Number of reads`/sum(data3_L1.longD$`Number of reads`))* 100

## --- SAMPLE 4 --- ##
data4_L1.longD$sample <- 4
data4_L1.longD$cpm <- (data4_L1.longD$`Number of reads`/sum(data4_L1.longD$`Number of reads`))* 10^6
data4_L1.longD$per <- (data4_L1.longD$`Number of reads`/sum(data4_L1.longD$`Number of reads`))* 100

## --- SAMPLE 5 --- ##
data5_L1.longD$sample <- 5
data5_L1.longD$cpm <- (data5_L1.longD$`Number of reads`/sum(data5_L1.longD$`Number of reads`))* 10^6
data5_L1.longD$per <- (data5_L1.longD$`Number of reads`/sum(data5_L1.longD$`Number of reads`))* 100

## --- SAMPLE 6 --- ##
data6_L1.longD$sample <- 6
data6_L1.longD$cpm <- (data6_L1.longD$`Number of reads`/sum(data6_L1.longD$`Number of reads`))* 10^6
data6_L1.longD$per <- (data6_L1.longD$`Number of reads`/sum(data6_L1.longD$`Number of reads`))* 100

## --- Culex sp: Join all data frames --- ##

Culex_sp.domains <- rbind(data1_L1.longD,data2_L1.longD,data3_L1.longD,data4_L1.longD,data5_L1.longD,data6_L1.longD)

## Choose colors
Culex_sp.domains$color <- "#999999"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Archaea"] <- "#CC99CC"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Eukaryota"] <- "#66CCCC"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Bacteria"] <- "#0C7BDC"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Baculoviridae"] <- "#DC3220"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Duplodnaviria"] <- "#E66100"
Culex_sp.domains$color[Culex_sp.domains$Domain == "environmental samples"] <- "#E1BE6A"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Monodnaviria"] <- "#D35FB7"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Nudiviridae"] <- "#CCCCFF"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Polydnaviridae"] <- "#000066"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Riboviria"] <- "#339900"
Culex_sp.domains$color[Culex_sp.domains$Domain == "unclassified viruses"] <- "#00CCCC"
Culex_sp.domains$color[Culex_sp.domains$Domain == "Varidnaviria"] <- "#663399"

## --- Plot of domains --- ##
## --- cpm --- ##
ggplot(Culex_sp.domains, aes(x = as.character(sample), y = cpm, fill =Domain)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Counts per million") +
  theme(legend.title = element_text(), legend.position = "top") + labs(fill="Domains") +  scale_fill_manual(values = Culex_sp.domains$color)

## --- percentage --- ##
ggplot(Culex_sp.domains, aes(x = as.character(sample), y = per, fill =Domain)) + geom_col(aes()) + theme_bw()  + xlab("Samples") + ylab("Percentage of reads") +
  theme(legend.title = element_text(), legend.position = "top") + labs(fill="Domain") + scale_fill_manual(values = Culex_sp.domains$color)

## --- G. affinis: Domains --- ##

## --- SAMPLE 7 --- ##
data7_L1.longD$sample <- 7
data7_L1.longD$cpm <- (data7_L1.longD$`Number of reads`/sum(data7_L1.longD$`Number of reads`))* 10^6
data7_L1.longD$per <- (data7_L1.longD$`Number of reads`/sum(data7_L1.longD$`Number of reads`))* 100

## --- SAMPLE 8 --- ##
data8_L1.longD$sample <- 8
data8_L1.longD$cpm <- (data8_L1.longD$`Number of reads`/sum(data8_L1.longD$`Number of reads`))* 10^6
data8_L1.longD$per <- (data8_L1.longD$`Number of reads`/sum(data8_L1.longD$`Number of reads`))* 100

## --- SAMPLE 9 --- ##
data9_L1.longD$sample <- 9
data9_L1.longD$cpm <- (data9_L1.longD$`Number of reads`/sum(data9_L1.longD$`Number of reads`))* 10^6
data9_L1.longD$per <- (data9_L1.longD$`Number of reads`/sum(data9_L1.longD$`Number of reads`))* 100

## --- SAMPLE 10 --- ##
data10_L1.longD$sample <- 10
data10_L1.longD$cpm <- (data10_L1.longD$`Number of reads`/sum(data10_L1.longD$`Number of reads`))* 10^6
data10_L1.longD$per <- (data10_L1.longD$`Number of reads`/sum(data10_L1.longD$`Number of reads`))* 100

## --- SAMPLE 11 --- ##
data11_L1.longD$sample <- 11
data11_L1.longD$cpm <- (data11_L1.longD$`Number of reads`/sum(data11_L1.longD$`Number of reads`))* 10^6
data11_L1.longD$per <- (data11_L1.longD$`Number of reads`/sum(data11_L1.longD$`Number of reads`))* 100

## --- SAMPLE 12 --- ##
data12_L1.longD$sample <- 12
data12_L1.longD$cpm <- (data12_L1.longD$`Number of reads`/sum(data12_L1.longD$`Number of reads`))* 10^6
data12_L1.longD$per <- (data12_L1.longD$`Number of reads`/sum(data12_L1.longD$`Number of reads`))* 100

## --- G. affinis: Join all data frames --- ##

Affinis.domains <- rbind(data7_L1.longD,data8_L1.longD,data9_L1.longD,data10_L1.longD,data11_L1.longD,data12_L1.longD)

## Choose colors
Affinis.domains$color[Affinis.domains$Domain == "Unclassified"] <- "#999999"
Affinis.domains$color[Affinis.domains$Domain == "Archaea"] <- "#CC99CC"
Affinis.domains$color[Affinis.domains$Domain == "Eukaryota"] <- "#66CCCC"
Affinis.domains$color[Affinis.domains$Domain == "Bacteria"] <- "#0C7BDC"
Affinis.domains$color[Affinis.domains$Domain == "Baculoviridae"] <- "#DC3220"
Affinis.domains$color[Affinis.domains$Domain == "Duplodnaviria"] <- "#E66100"
Affinis.domains$color[Affinis.domains$Domain == "environmental samples"] <- "#E1BE6A"
Affinis.domains$color[Affinis.domains$Domain == "Monodnaviria"] <- "#D35FB7"
Affinis.domains$color[Affinis.domains$Domain == "unclassified viruses"] <- "#CCCCFF"
Affinis.domains$color[Affinis.domains$Domain == "Polydnaviridae"] <- "#000066"
Affinis.domains$color[Affinis.domains$Domain == "Riboviria"] <- "#339900"
Affinis.domains$color[Affinis.domains$Domain == "unclassified bacterial viruses"] <- "#00CCCC"
Affinis.domains$color[Affinis.domains$Domain == "Varidnaviria"] <- "#663399"

## --- Plot of domains --- ##
## --- cpm --- ##
ggplot(Affinis.domains, aes(x = as.character(sample), y = cpm, fill =Domain)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Samples") + ylab("Counts per million") +
  theme(legend.title = element_text(), legend.position = "top") + labs(fill="Domains") #+  scale_fill_manual(values = Affinis.domains$color)

## --- percentage --- ##
ggplot(Affinis.domains, aes(x = as.character(sample), y = per, fill =Domain)) + geom_col(aes()) + theme_bw()  + xlab("Samples") + ylab("Percentage of reads") +
  theme(legend.title = element_text(), legend.position = "top") + labs(fill="Domain") + scale_fill_manual(values = Affinis.domains$color)


## ---- VIRUSES: IN-DEPTH --- ##

## --- Culex sp. --- ##
## --- SAMPLE 1 --- ##
bu.1 <- all_viruses1[all_viruses1$Family == "Bunyavirales",]
bu.1_reads <- data.frame(bu.1$`Number of reads`, bu.1$Species)
colnames(bu.1_reads) <- c("Reads","Species")
bu.1_reads <- na.omit(bu.1_reads)
bu.1_reads <- aggregate(.~Species,data=bu.1_reads,FUN=sum)
bu.1 <- data.frame(unique(bu.1$Species))
colnames(bu.1) <- "Species"


## --- SAMPLE 2 --- ##
bu.2 <- all_viruses2[all_viruses2$Family == "Bunyavirales",]
bu.2_reads <- data.frame(bu.2$`Number of reads`, bu.2$Species)
colnames(bu.2_reads) <- c("Reads","Species")
bu.2_reads <- na.omit(bu.2_reads)
bu.2_reads <- aggregate(.~Species,data=bu.2_reads,FUN=sum)
bu.2 <- data.frame(unique(bu.2$Species))
colnames(bu.2) <- "Species"

## --- SAMPLE 3 --- ##
bu.3 <- all_viruses3[all_viruses3$Family == "Bunyavirales",]
bu.3_reads <- data.frame(bu.3$`Number of reads`, bu.3$Species)
colnames(bu.3_reads) <- c("Reads","Species")
bu.3_reads <- na.omit(bu.3_reads)
bu.3_reads <- aggregate(.~Species,data=bu.3_reads,FUN=sum)
bu.3 <- data.frame(unique(bu.3$Species))
colnames(bu.3) <- "Species"

## --- SAMPLE 4 --- ##
bu.4 <- all_viruses4[all_viruses4$Family == "Bunyavirales",]
bu.4_reads <- data.frame(bu.4$`Number of reads`, bu.4$Species)
colnames(bu.4_reads) <- c("Reads","Species")
bu.4_reads <- na.omit(bu.4_reads)
bu.4_reads <- aggregate(.~Species,data=bu.4_reads,FUN=sum)
bu.4 <- data.frame(unique(bu.4$Species))
colnames(bu.4) <- "Species"

## --- SAMPLE 5 --- ##
bu.5<- all_viruses5[all_viruses5$Family == "Bunyavirales",]
bu.5_reads <- data.frame(bu.5$`Number of reads`, bu.5$Species)
colnames(bu.5_reads) <- c("Reads","Species")
bu.5_reads <- na.omit(bu.5_reads)
bu.5_reads <- aggregate(.~Species,data=bu.5_reads,FUN=sum)
bu.5 <- data.frame(unique(bu.5$Species))
colnames(bu.5) <- "Species"

## --- SAMPLE 6 --- ##
bu.6 <- all_viruses6[all_viruses6$Family == "Bunyavirales",]
bu.6_reads <- data.frame(bu.6$`Number of reads`, bu.6$Species)
colnames(bu.6_reads) <- c("Reads","Species")
bu.6_reads <- na.omit(bu.6_reads)
bu.6_reads <- aggregate(.~Species,data=bu.6_reads,FUN=sum)
bu.6 <- data.frame(unique(bu.6$Species))
colnames(bu.6) <- "Species"

## --- Join species --- ##
bu.1_2 <- rbind(bu.1, bu.2)
bu.all <- rbind(bu.3, bu.4, bu.5, bu.6)
bu.all_unique <- data.frame(unique(bu.all))

#### -------- VENN DIAGRAM -------- ##

## prepare data for venn diagram
mos <- list(
  M3 = bu.3_reads$Species, 
  M4 = bu.4_reads$Species, 
  M5 = bu.5_reads$Species,
  M6 = bu.6_reads$Species
)

## make venn diagram
ggvenn(
  mos, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

## Compute intersections between data frames to find out common viruses
common.vir2 <- inner_join( bu.3_reads,bu.4_reads, bu.5_reads, bu.6_reads )
x <-data.frame(bu.3_reads$Species)
colnames(x) <- "Species"
y <- data.frame(bu.5_reads$Species )
colnames(y) <- "Species"
common.vir3 <- inner_join(x , y)

## compute intersection between the two groups of mosquitoes
x <-data.frame(bu.1_reads$Species)
colnames(x) <- "Species"
y <- data.frame(rbind(bu.3_reads,bu.4_reads, bu.5_reads, bu.6_reads))
y <- data.frame(unique(y$Species))
colnames(y) <- "Species"
common.vir4 <- inner_join(x , y)

mos2 <- list(
  M1.2 = bu.1_reads$Species, 
  M3.4.5.6 = c(y$Species)
)
## make venn diagram
ggvenn(
  mos2, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 2
)


## ------- IN DEPTH ANALYSIS OF G. AFFINIS' VIROME ------ ##
setwd("~/Desktop/Master_VHIR/TFM/")
## --- By host --- ##
host <- read_xlsx("host_virome_fish.xlsx")
total_virome <- nrow(unique(host))

freq <- data.frame(table(host$Host))
write.csv(file = "Frequency_host.csv", x = freq)
freq$percentage <- (freq$Freq/total_virome )*100
colnames(freq)[1] <- "Host"

nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

freq$labels <- paste(round(freq$percentage),"%", sep = "")
df2 <- freq %>% 
  mutate(
    cs = rev(cumsum(rev(percentage))), 
    pos = percentage/2 + lead(cs, 1),
    pos = if_else(is.na(pos), percentage/2, pos))

ggplot(freq, aes(x="", y=percentage, fill=Host)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_manual(values = mycolors) + theme(legend.text = element_text(size=12), legend.title = element_text(size=14))



 # geom_label_repel(data = df2,
                   aes(y = pos, label = labels),
                   size = 4.5, nudge_x = 1, show.legend = FALSE)

## --- By disease --- ##
disease <- read_xlsx("Host_disease_fish.xlsx")
total_virome <- length(unique(disease$Disease))

freq <- data.frame(table(disease$Disease))
freq$percentage <- (freq$Freq/total_virome )*100
colnames(freq)[1] <- "Disease"
write.csv(file ="disease_affinis.csv", x=freq)

nb.cols <- 26
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

freq$labels <- paste(round(freq$percentage),"%", sep = "")


ggplot(freq, aes(x="", y=percentage, fill=Disease)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_manual(values = mycolors)  + guides(fill= guide_legend(ncol = 1, byrow = FALSE)) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=14))



# geom_label_repel(data = df2,
aes(y = pos, label = labels),
size = 4.5, nudge_x = 1, show.legend = FALSE)


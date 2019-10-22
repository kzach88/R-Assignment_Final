---
title: "Final R"
author: "Zach"
date: "October 21, 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##Install needed Libraries

```{r}

library(tidyverse)
library(reshape2)
```
##Import the original data files

```{r}
fang <- read.table("https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep="\t", stringsAsFactors = FALSE)
snps <- read.table("https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/snp_position.txt",  header = TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE)
```
#File inspection using the folowing commands
```{r}
str(fang)
str(snps)
nrow(fang)
nrow(fang)
ncol(fang)
nrow(snps)
ncol(snps)
summary(fang)
summary(snps)
typeof(fang)
typeof(snps)
class(fang)
class(snps)
colnames(fang)
colnames(snps)
str(fang)
str(snps)
head(fang)
head(snps)

```
##After inspecting these files, i observed the following

fang_et_al_genotypes.txt

Size: 11.05MB

This file has 986 colums and 2782 rowsum;

This file has 16 Groups in the Group column:

Based on the head command, genotypes has missing values coded with '?

snp_position.txt

This file is 82.76KB 

It has 986 rows and 15 columns

This file has 339 candidates and 644 random SNPS. 

SNP Position has column names for SNP ID, marker ID, Chromosome, Position, alternative and multiple positions, amplicon, feature name, gene'


##PART 2

##Data Processing

First, I will subset the fang dataframe into two objects, one for maize (ZMMIL, ZMMLR, and ZMMMR) and one for teosinte (ZMPBA, ZMPIL, and ZMPJA) ans then select Maize or Teosinte Groups based on Group column*

```{r}

maizegroup <- c("ZMMIL", "ZMMLR", "ZMMMR")
teosintegroup <- c("ZMPBA", "ZMPIL", "ZMPJA")
maizegenotypes <-fang[fang$Group %in% maizegroup, ]
row.names(maizegenotypes) <- maizegenotypes[,1]
teosintegenotypes <- fang[fang$Group %in% teosintegroup, ]
row.names(teosintegenotypes) <- teosintegenotypes[,1]
```
##Remove unwanted columns in maizefang and teosintefang dataframes (Sample_ID, JG_OTU, and Group). leaveonly the names of the SNPs and the data.
```{r}
cutmaizegenotypes <- maizegenotypes[,-1:-3]
cutteosintegenotypes <- teosintegenotypes[,-1:-3]
```
##To remove unwanted columns in the snps dataframe i will use cut command (cdv marker ID and everything after position).
```{r}
cutposition <-snps[, c("SNP_ID","Chromosome","Position")]
```
##Transpose the genotypes dataframes, making sure that the final result is a dataframe**  

```{r}

transposedmaize <- as.data.frame(t(cutmaizegenotypes))

transposedteosinte <- as.data.frame(t(cutteosintegenotypes))

is.data.frame(transposedmaize)

is.data.frame(transposedteosinte)
```

#Sort the genotype and SNP position dataframes by the SNP name, the column that the two files have in common
```{r}
sortedposition <- cutposition[order(cutposition$SNP_ID),] #sort by SNP_ID

SNPstransposedmaize <- cbind(SNP_ID = rownames(transposedmaize), transposedmaize) 
```
##For the transposed genotype files, we need to create a new column with the SNP IDs and 
```{r}
rownames(SNPstransposedmaize) <- NULL
```
##Delete the rownames created from transposing the files. That way, we can sort the new column of SNP IDs
```{r}
SNPstransposedteosinte <- cbind(SNP_ID = rownames(transposedteosinte), transposedteosinte)

rownames(SNPstransposedteosinte) <- NULL

sortedmaize <- SNPstransposedmaize[order(SNPstransposedmaize$SNP_ID),] #sort by SNP_ID

sortedteosinte <- SNPstransposedteosinte[order(SNPstransposedteosinte$SNP_ID),] #sort by SNP_ID
```
##Join the sortedposition dataframe to each of the genotype dataframes using SNP_ID
```{r}
table(sortedposition$SNP_ID %in% sortedmaize$SNP_ID)
table(sortedposition$SNP_ID %in% sortedteosinte$SNP_ID)
joinedmaize <- merge(sortedposition, sortedmaize, by.x="SNP_ID", by.y="SNP_ID") 
```
##now merge by the SNP_ID column
```{r}
joinedteosinte <- merge(sortedposition, sortedteosinte, by.x="SNP_ID", by.y="SNP_ID")
```
##Next, I will isolate each chromosome and sort by position. For the first 10 files, we need 1 for each chromosome with SNPs ordered by increasing position values and missing data shown by '?'. To do this, first order the datasets by increasing position
```{r}
library(gtools)
orderedmaizeincrease <- joinedmaize[mixedorder(joinedmaize$snps),] 
```
##order the dataset by increasing position
```{r}
orderedteosinteincrease <- joinedteosinte[mixedorder(joinedteosinte$Position),]
```
##It is also necessary to write a function to pull out each chromosome and write them to a file using the package dplyr This was already loaded in tydiverse
```{r}
library(dplyr)
maize_no_ambig <- subset(orderedmaizeincrease, orderedmaizeincrease$Chromosome!="unknown" & orderedmaizeincrease$Chromosome!="multiple") 
```
##I removed character strings from the Chromosome column of the dataset, leaving only the numeric chromosome values for the function
```{r}
MaizeChromosomeQ <- tbl_df(maize_no_ambig)
```
##need to translate the dataframe into dplyr-readable format
```{r}
filewriteMQ  = function(DF) {

write.table(DF,file = paste0("Maize_Chromosome_Q",unique(DF$Chromosome),".txt"), sep = "\t", row.names = FALSE)

return(DF)

} #my function to write files with names based on chromosome


MaizeChromosomeQ %>% 

group_by(Chromosome) %>% 

do(filewriteMQ(.)) #this is the dplyr piping of the dataframe to a grouping function, and then the files are written for each group
```
##I did the same for teosinte, writing a modified function to change the filenames
```{r}
teosinte_no_ambig <- subset(orderedteosinteincrease, orderedteosinteincrease$Chromosome!="unknown" & orderedteosinteincrease$Chromosome!="multiple")

TeosinteChromosomeQ <- tbl_df(teosinte_no_ambig) 



filewriteTQ  = function(DF) {

write.table(DF,file = paste0("Teosinte_Chromosome_Q",unique(DF$Chromosome),".txt"), sep = "\t", row.names = FALSE)

return(DF)

} 
TeosinteChromosomeQ %>% 

group_by(Chromosome) %>% 

do(filewriteTQ(.)) 

````
##For the next ten files of maize and teosinte with '-' denoting missing values:
```{r}
missingvaluemaize <- as.data.frame(lapply(joinedmaize, function(x) {gsub("\\?","-", x)})) #replace missing values with -

missingvalueteosinte <- as.data.frame(lapply(joinedteosinte, function(x) {gsub("\\?","-", x)}))
```
##Then, I sorted the files by decreasing order:
```{r}
orderedmaizedecrease <- missingvaluemaize[mixedorder(as.character(missingvaluemaize$Position), decreasing=TRUE),] 
```
#order the dataset by decreasing position
```{r}
orderedteosintedecrease <- missingvalueteosinte[mixedorder(as.character(missingvalueteosinte$Position), decreasing=TRUE),]
```
##Next, using the function from before, I pulled out each chromosome and created a file for maize chromosomes
```{r}
maize_no_ambig2 <- subset(orderedmaizedecrease, orderedmaizedecrease$Chromosome!="unknown" & orderedmaizedecrease$Chromosome!="multiple")

MaizeChromosomeD <- tbl_df(maize_no_ambig2) 



filewriteMD  = function(DF) {

write.table(DF,file = paste0("Maize_Chromosome_D",unique(DF$Chromosome),".txt"), sep = "\t", row.names = FALSE)

return(DF)

} 
MaizeChromosomeD %>% 

group_by(Chromosome) %>% 

do(filewriteMD(.))

teosinte_no_ambig2 <- subset(orderedteosintedecrease, orderedteosintedecrease$Chromosome!="unknown" & orderedteosintedecrease$Chromosome!="multiple")

TeosinteChromosomeD <- tbl_df(teosinte_no_ambig2) 



filewriteTD  = function(DF) {

write.table(DF,file = paste0("Teosinte_Chromosome_D",unique(DF$Chromosome),".txt"), sep = "\t", row.names = FALSE)

return(DF)

} 

TeosinteChromosomeD %>% 

group_by(Chromosome) %>% 

do(filewriteTD(.))
```
#Part TWO
##Question 1: SNPs per Chromosome
###Graph
#In this part i am to visualize the data in different ways. To start iaam retaking the joint dataframes without removing the missing data and NA values, so as to crecreate the joint files, using a previous command
```{r}
transposed_fang <- as.data.frame(t(fang[,-1]))
colnames(transposed_fang) <- fang$Sample_ID 
snp_id_chr_pos <-snps[,c(1,3,4)]
fang_joint <- merge(x = snp_id_chr_pos, y = transposed_fang, by.x = "SNP_ID", by.y ="row.names", all.y = TRUE)
```
##load packages 
```{r}
library(reshape2)
library(ggplot2)
library(plyr)
```

##Now run the codes 
```{r}
fang_joint$Chromosome <- factor(fang_joint$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown", "NA"))



ggplot(genotypes_joint, aes((Chromosome))) + geom_bar()
```
##Then, we can plot the SNPs according to the groups they belong, using the genotypes dataframe. That will tell us What groups contribute most of these SNPs?
```{r}
ggplot(fang, aes(Group)) +

geom_bar()
"The ZMMIL, ZMMLR and ZMPBA contribute the most to the SNP count, see graph below"
```
##to analyze the amount of heterozygicity we melt the datasets to make them tidy, using a vector with the names of all the SNPs
```{r}
headers_names<- colnames(fang)[-c(1:3)]

genotypes_melted <- melt(fang, measure.vars = headers_names)
```
###Then create a new column to indicate whether a particular site is homozygous, that is A/A, C/C, G/G, T/T) or heterozygous (otherwise)). First we assign all missing values as NA:
```{r}
genotypes_melted[ genotypes_melted == "?/?" ] = NA
genotypes_melted$isHomozygous <- (genotypes_melted$value=="A/A" | genotypes_melted$value=="C/C" | genotypes_melted$value=="G/G" | genotypes_melted$value=="T/T")
```
##Next, i will sort the dataframe using Group and Species_ID values. 
```{r}
genotypes_sorted_by_ID <- genotypes_melted[order(genotypes_melted$Sample_ID),]
genotypes_sorted_by_Group <- genotypes_melted[order(genotypes_melted$Group),]
```
#And at last, make a graph that shows the proportion of homozygous and heterozygous sites as well as missing data in each species (you won't be able to see species names). For doing that, we first built a new dataframe with all the counting values (for homozygous, heterozygous and NA) per Sample_ID, and then we melt the results
```{r}
counting_ID <- ddply(genotypes_sorted_by_ID, c("Sample_ID"), summarise, counting_homozygous=sum(isHomozygous, na.rm=TRUE), counting_heterozygous=sum(!isHomozygous, na.rm=TRUE), isNA=sum(is.na(isHomozygous)))

counting_ID_melt <- melt(counting_ID, measure.vars = c("counting_homozygous", "counting_heterozygous", "isNA"))
ggplot(counting_ID_melt,aes(x = Sample_ID, y= value, fill=variable)) + geom_bar(stat = "identity", position = "stack")

"The plot shows how the homozygous counting is bigger through all SNPs, having low counting for missing data. There are also a group of SNPs that seem not to have heterozygous alleles.
```

#The same process was made when the data is sorted according to the groups:
```{r}
counting_Group <- ddply(genotypes_sorted_by_Group, c("Group"), summarise, counting_homozygous=sum(isHomozygous, na.rm=TRUE), counting_heterozygous=sum(!isHomozygous, na.rm=TRUE), isNA=sum(is.na(isHomozygous)))


counting_Group_melt <- melt(counting_Group, measure.vars = c("counting_homozygous", "counting_heterozygous", "isNA"))

ggplot(counting_Group_melt,aes(x = Group, y= value, fill=variable)) + geom_bar(stat = "identity", position ="stack")

"This graph shows how the groups that contribute the most to the SNP number also contribute to the number of heterozygous and homozygous""
```
##Finally, a visualization of the data based on the calculations of the observed heterozygocity per locus. 

#I will calculate the observed heterozygocity per SNP as the rate between the number of heterozygous divided by the total of genotyped individuals. After that,I will construch  scatter plot. First, i willmelt the `fang_joint` file, using as measure variables all the genotyped individuals
```{r}
headers_names_joint<- colnames(fang_joint)[-c(1:3)]

genotypes_melted_joint <- melt(fang_joint, measure.vars = headers_names_joint)
```
##Then i will calculate the number of heterozygous, according to the SNP_ID:
```{r}
genotypes_melted_joint[ genotypes_melted_joint == "?/?" ] = NA
genotypes_melted_joint$isHomozygous <- (genotypes_melted_joint$value=="A/A" | genotypes_melted_joint$value=="C/C" | genotypes_melted_joint$value=="G/G" | genotypes_melted_joint$value=="T/T")
genotypes_sorted_by_SNP <- genotypes_melted_joint[order(genotypes_melted_joint$SNP_ID),]
```
#After this, we can calculate the observed heterozygosity per SNP and make our plot**
```{r}
Observed_Het_per_locus <- ddply(genotypes_sorted_by_SNP, c("SNP_ID"), summarise, heterozygocity_count=sum(!isHomozygous, na.rm=TRUE), total_count=sum(!is.na(isHomozygous)))

Observed_Het_per_locus$Obs_heterozygocity <- (Observed_Het_per_locus$heterozygocity_count/Observed_Het_per_locus$total_count)
Observed_Het_per_locus_melt <- melt(Observed_Het_per_locus, measure.vars = "Obs_heterozygocity")
ggplot(Observed_Het_per_locus_melt,aes(x = SNP_ID, y= value, fill=variable)) + geom_point()

"The observed heterozygocity range for each SNP, is mainly between 0 and 0.27, and one of them has a velue of 1.0"
```


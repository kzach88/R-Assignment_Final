---
title: "RcodeANDOutput"
author: "Zach"
date: "October 22, 2019"
output: pdf_document
---
##Import the original data files
```{r}
fang <- read.table("https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep="\t", stringsAsFactors = FALSE)
snps <- read.table("https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/snp_position.txt",  header = TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE)
```
```{r}
library(tidyverse)
library(reshape2)
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

Size: 11.05MB"
This file has 986 colums and 2782 rows
This file has 16 Groups in the Group column"
Based on the head command, genotypes has missing values coded with '?

snp_position.txt
This file is 82.76KB"
It has 986 rows and 15 columns"
This file has 339 candidates and 644 random SNPS"
SNP Position has column names for SNP ID, marker ID, Chromosome, Position, alternative and multiple positions, amplicon, feature name, gene"

##PART 2

##Data Processing

##First, I will subset the fang dataframe into two objects, one for maize (ZMMIL, ZMMLR, and ZMMMR) and one for teosinte (ZMPBA, ZMPIL, and ZMPJA) ans then select Maize or Teosinte Groups based on Group column"

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
##join the transposed files with the snp_id_chr_pos dataframe
```{r}
maize_joint <- merge(x = snp_id_chr_pos, y = transposedmaize, by.x = "SNP_ID", by.y ="row.names", all.y = TRUE)
teosinte_joint <- merge(x = snp_id_chr_pos, y = transposedteosinte, by.x = "SNP_ID", by.y ="row.names", all.y = TRUE)
```
##remove the rows that lack of information for chromosome and position
```{r}

teosinte_joint$Position <- as.numeric(as.character(teosinte_joint$Position))

teosinte_joint$Chromosome <- as.numeric(as.character(teosinte_joint$Chromosome))



teosinte_joint <- teosinte_joint[!(teosinte_joint$Position== "unknown") | !(teosinte_joint$Chromosome== "unknown") | !(teosinte_joint$Chromosome== "multiple"),]

teosinte_joint <- teosinte_joint[!(is.na(teosinte_joint$Position)),]



maize_joint$Position <- as.numeric(as.character(maize_joint$Position))

maize_joint$Chromosome <- as.numeric(as.character(maize_joint$Chromosome))



maize_joint <- tripsacum_joint[!(maize_joint$Position== "unknown") | !(maize_joint$Chromosome== "unknown") | !(maize_joint$Chromosome== "multiple"),]

maize_joint <- maize_joint[!(is.na(maize_joint$Position)),]

```
##sort the files
```{r}
teosinte_joint <- teosinte_joint[order(teosinte_joint$Position),]

maize_joint <- maize_joint[order(maize_joint$Position),]
```
##split the dataframes according to the chromosome number
```{r}
teosinte_split <- split(teosinte_joint, teosinte_joint$Chromosome)

allNames <- names(teosinte_split)

 for(thisName in allNames){

     saveName = paste0('Chr_', thisName, '_teosinteIncreasing.txt')

     write.table(teosinte_split[[thisName]], file = saveName, quote = FALSE, sep="\t", row.names = FALSE)

 }



maize_split <- split(maize_joint, maize_joint$Chromosome)

allNames <- names(maize_split)

 for(thisName in allNames){

     saveName = paste0('Chr_', thisName, '_maizeIncreasing.txt')

     write.table(maize_split[[thisName]], file = saveName, quote = FALSE, sep="\t", row.names = FALSE)

 }
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

```{r}
teosinte_joint_Decrease <- teosinte_joint[order(teosinte_joint$Position, decreasing=T),]

maize_joint_Decrease <- maize_joint[order(maize_joint$Position, decreasing=T),]
```
##now replace the values and create final files
```{r}
teosinte_joint_Decrease[] <- lapply(teosinte_joint_Decrease, as.character)

teosinte_joint_Decrease[teosinte_joint_Decrease == '?/?'] <- '-/-'



maize_joint_Decrease[] <- lapply(maize_joint_Decrease, as.character)

maize_joint_Decrease[maize_joint_Decrease == '?/?'] <- '-/-'

teosinte_split_Decrease <- split(teosinte_joint_Decrease, teosinte_joint_Decrease$Chromosome)

allNames <- names(teosinte_split_Decrease)

 for(thisName in allNames){

     saveName = paste0('Chr_', thisName, '_teosinte_Decreasing.txt')

     write.table(teosinte_split_Decrease[[thisName]], file = saveName, quote = FALSE, sep="\t", row.names = FALSE)

 }

maize_split_Decrease <- split(maize_joint_Decrease, maize_joint_Decrease$Chromosome)

allNames <- names(maize_split_Decrease)

 for(thisName in allNames){

     saveName = paste0('Chr_', thisName, '_maize_Decreasing.txt')

     write.table(maize_split_rev[[thisName]], file = saveName, quote = FALSE, sep="\t", row.names = FALSE)

 }
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
ggplot(fang_joint, aes((Chromosome))) + geom_bar()
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

"The plot shows how the homozygous counting is bigger through all SNPs, having low counting for missing data. There are also a group of SNPs that seem not to have heterozygous alleles"
```

#The same process was made when the data is sorted according to the groups:
```{r}
counting_Group <- ddply(genotypes_sorted_by_Group, c("Group"), summarise, counting_homozygous=sum(isHomozygous, na.rm=TRUE), counting_heterozygous=sum(!isHomozygous, na.rm=TRUE), isNA=sum(is.na(isHomozygous)))
counting_Group_melt <- melt(counting_Group, measure.vars = c("counting_homozygous", "counting_heterozygous", "isNA"))
ggplot(counting_Group_melt,aes(x = Group, y= value, fill=variable)) + geom_bar(stat = "identity", position ="stack")

"This graph shows how the groups that contribute the most to the SNP number also contribute to the number of heterozygous and homozygous"
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


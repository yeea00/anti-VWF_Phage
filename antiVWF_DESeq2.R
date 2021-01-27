# Analyze counts with DESeq2
library(DESeq2)
library(readxl)
library(data.table)
library(openxlsx)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import counts and sample information
All.counts <- lapply(excel_sheets("All_Counts.xlsx"), read_excel, path = "All_Counts.xlsx")
names(All.counts) <- excel_sheets("All_Counts.xlsx")  # rename df in list with sheet names
All.counts <- All.counts[-2] # remove replicate unselected phage counts (replicate A is paired with selected)
Info <- read.delim(file = "antiVWF_info.txt", header = T, sep = "\t")

# Separate counts based on strand and convert NA to 0
SS <- lapply(All.counts, function(x) x[,grepl("Location", colnames(x)) | grepl("SS", colnames(x))])  # extract location and counts for sense
SS <- plyr::join_all(SS, by = "Location")  # combine lists into dataframe
SS[is.na(SS)] <- 0  # NA to 0
colnames(SS)[-1] <- gsub("_SS", "", colnames(SS)[-1])  # remove "_SS" from colnames
AS <- lapply(All.counts, function(x) x[,grepl("Location", colnames(x)) | grepl("AS", colnames(x))])  # extract location and counts for antisense
AS <- plyr::join_all(AS, by = "Location")  # combine lists into dataframe
AS[is.na(AS)] <- 0  # NA to 0
colnames(AS)[-1] <- gsub("_AS", "", colnames(AS)[-1])  # remove "_AS" from colnames

# Run counts through DESeq2
dds.SS.mx <- DESeqDataSetFromMatrix(as.matrix(data.frame(SS[,-1], row.names = SS[,1])), 
                                      colData = Info, design = ~Treatment)
dds.SS <- DESeq(dds.SS.mx)
dds.AS.mx <- DESeqDataSetFromMatrix(as.matrix(data.frame(AS[,-1], row.names = AS[,1])), 
                                    colData = Info, design = ~Treatment)
dds.AS <- DESeq(dds.AS.mx)

#############################
# Negative vs. Input
#############################
# Sense strand
res.Neg.SS <- as.data.frame(results(dds.SS, pAdjustMethod = "fdr", independentFiltering=F, 
                                         contrast = c("Treatment", "Neg", "Input")))
res.Neg.SS <- setDT(res.Neg.SS, keep.rownames = T)
res.Neg.SS$log.padj <- -log10(res.Neg.SS$padj)
res.Neg.SS$Strand = "Sense"
res.Neg.SS <- rename(res.Neg.SS, c("rn" = "Location"))
res.Neg.SS <- plyr::join(All.counts[[1]][,1:6], res.Neg.SS, by = "Location")

# Antisense strand
res.Neg.AS <- as.data.frame(results(dds.AS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "Neg", "Input")))
res.Neg.AS <- setDT(res.Neg.AS, keep.rownames = T)
res.Neg.AS$log.padj <- -log10(res.Neg.AS$padj)
res.Neg.AS$Strand = "Antisense"
res.Neg.AS <- rename(res.Neg.AS, c("rn" = "Location"))
res.Neg.AS <- plyr::join(All.counts[[1]][,1:6], res.Neg.AS, by = "Location")

# Combine into single dataframe
res.Neg <- rbind(res.Neg.SS, res.Neg.AS)
res.Neg$Strand <- factor(res.Neg$Strand, levels = unique(res.Neg$Strand))
res.Neg$Frame <- substring(res.Neg$VWF_aa, regexpr("[.]", res.Neg$VWF_aa) + 1)

#############################
# VWFab vs. Input
#############################
# Sense strand
res.VWFab.SS <- as.data.frame(results(dds.SS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFab", "Input")))
res.VWFab.SS <- setDT(res.VWFab.SS, keep.rownames = T)
res.VWFab.SS$log.padj <- -log10(res.VWFab.SS$padj)
res.VWFab.SS$Strand = "Sense"
res.VWFab.SS <- rename(res.VWFab.SS, c("rn" = "Location"))
res.VWFab.SS <- plyr::join(All.counts[[1]][,1:6], res.VWFab.SS, by = "Location")

# Antisense strand
res.VWFab.AS <- as.data.frame(results(dds.AS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFab", "Input")))
res.VWFab.AS <- setDT(res.VWFab.AS, keep.rownames = T)
res.VWFab.AS$log.padj <- -log10(res.VWFab.AS$padj)
res.VWFab.AS$Strand = "Antisense"
res.VWFab.AS <- rename(res.VWFab.AS, c("rn" = "Location"))
res.VWFab.AS <- plyr::join(All.counts[[1]][,1:6], res.VWFab.AS, by = "Location")

# Combine into single dataframe
res.VWFab <- rbind(res.VWFab.SS, res.VWFab.AS)
res.VWFab$Strand <- factor(res.VWFab$Strand, levels = unique(res.VWFab$Strand))
res.VWFab$Frame <- substring(res.VWFab$VWF_aa, regexpr("[.]", res.VWFab$VWF_aa) + 1)

#############################
# VWFi_I_1 vs. Input
#############################
# Sense strand
res.VWFi_I_1.SS <- as.data.frame(results(dds.SS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFi_I_1", "Input")))
res.VWFi_I_1.SS <- setDT(res.VWFi_I_1.SS, keep.rownames = T)
res.VWFi_I_1.SS$log.padj <- -log10(res.VWFi_I_1.SS$padj)
res.VWFi_I_1.SS$Strand = "Sense"
res.VWFi_I_1.SS <- rename(res.VWFi_I_1.SS, c("rn" = "Location"))
res.VWFi_I_1.SS <- plyr::join(All.counts[[1]][,1:6], res.VWFi_I_1.SS, by = "Location")

# Antisense strand
res.VWFi_I_1.AS <- as.data.frame(results(dds.AS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFi_I_1", "Input")))
res.VWFi_I_1.AS <- setDT(res.VWFi_I_1.AS, keep.rownames = T)
res.VWFi_I_1.AS$log.padj <- -log10(res.VWFi_I_1.AS$padj)
res.VWFi_I_1.AS$Strand = "Antisense"
res.VWFi_I_1.AS <- rename(res.VWFi_I_1.AS, c("rn" = "Location"))
res.VWFi_I_1.AS <- plyr::join(All.counts[[1]][,1:6], res.VWFi_I_1.AS, by = "Location")

# Combine into single dataframe
res.VWFi_I_1 <- rbind(res.VWFi_I_1.SS, res.VWFi_I_1.AS)
res.VWFi_I_1$Strand <- factor(res.VWFi_I_1$Strand, levels = unique(res.VWFi_I_1$Strand))
res.VWFi_I_1$Frame <- substring(res.VWFi_I_1$VWF_aa, regexpr("[.]", res.VWFi_I_1$VWF_aa) + 1)

#############################
# VWFi_II_1 vs. Input
#############################
# Sense strand
res.VWFi_II_1.SS <- as.data.frame(results(dds.SS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFi_II_1", "Input")))
res.VWFi_II_1.SS <- setDT(res.VWFi_II_1.SS, keep.rownames = T)
res.VWFi_II_1.SS$log.padj <- -log10(res.VWFi_II_1.SS$padj)
res.VWFi_II_1.SS$Strand = "Sense"
res.VWFi_II_1.SS <- rename(res.VWFi_II_1.SS, c("rn" = "Location"))
res.VWFi_II_1.SS <- plyr::join(All.counts[[1]][,1:6], res.VWFi_II_1.SS, by = "Location")

# Antisense strand
res.VWFi_II_1.AS <- as.data.frame(results(dds.AS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFi_II_1", "Input")))
res.VWFi_II_1.AS <- setDT(res.VWFi_II_1.AS, keep.rownames = T)
res.VWFi_II_1.AS$log.padj <- -log10(res.VWFi_II_1.AS$padj)
res.VWFi_II_1.AS$Strand = "Antisense"
res.VWFi_II_1.AS <- rename(res.VWFi_II_1.AS, c("rn" = "Location"))
res.VWFi_II_1.AS <- plyr::join(All.counts[[1]][,1:6], res.VWFi_II_1.AS, by = "Location")

# Combine into single dataframe
res.VWFi_II_1 <- rbind(res.VWFi_II_1.SS, res.VWFi_II_1.AS)
res.VWFi_II_1$Strand <- factor(res.VWFi_II_1$Strand, levels = unique(res.VWFi_II_1$Strand))
res.VWFi_II_1$Frame <- substring(res.VWFi_II_1$VWF_aa, regexpr("[.]", res.VWFi_II_1$VWF_aa) + 1)

#############################
# VWFi_II_2 vs. Input
#############################
# Sense strand
res.VWFi_II_2.SS <- as.data.frame(results(dds.SS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFi_II_2", "Input")))
res.VWFi_II_2.SS <- setDT(res.VWFi_II_2.SS, keep.rownames = T)
res.VWFi_II_2.SS$log.padj <- -log10(res.VWFi_II_2.SS$padj)
res.VWFi_II_2.SS$Strand = "Sense"
res.VWFi_II_2.SS <- rename(res.VWFi_II_2.SS, c("rn" = "Location"))
res.VWFi_II_2.SS <- plyr::join(All.counts[[1]][,1:6], res.VWFi_II_2.SS, by = "Location")

# Antisense strand
res.VWFi_II_2.AS <- as.data.frame(results(dds.AS, pAdjustMethod = "fdr", independentFiltering=F, 
                                    contrast = c("Treatment", "VWFi_II_2", "Input")))
res.VWFi_II_2.AS <- setDT(res.VWFi_II_2.AS, keep.rownames = T)
res.VWFi_II_2.AS$log.padj <- -log10(res.VWFi_II_2.AS$padj)
res.VWFi_II_2.AS$Strand = "Antisense"
res.VWFi_II_2.AS <- rename(res.VWFi_II_2.AS, c("rn" = "Location"))
res.VWFi_II_2.AS <- plyr::join(All.counts[[1]][,1:6], res.VWFi_II_2.AS, by = "Location")

# Combine into single dataframe
res.VWFi_II_2 <- rbind(res.VWFi_II_2.SS, res.VWFi_II_2.AS)
res.VWFi_II_2$Strand <- factor(res.VWFi_II_2$Strand, levels = unique(res.VWFi_II_2$Strand))
res.VWFi_II_2$Frame <- substring(res.VWFi_II_2$VWF_aa, regexpr("[.]", res.VWFi_II_2$VWF_aa) + 1)

#############################
# Export all results to excel sheet
#############################
res.all <- list("Neg" = res.Neg, "VWFab" = res.VWFab,
                   "VWFi_I_1" = res.VWFi_I_1, "VWFi_II_1" = res.VWFi_II_1, "VWFi_II_2" = res.VWFi_II_2)  # define Excel sheets
write.xlsx(res.all, file = "DESeq2 Results.xlsx")  # write Excel file

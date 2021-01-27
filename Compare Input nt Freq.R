# Quantify nucleotide frequency at each position in unselected phage
library(data.table)
library(readxl)
library(ggplot2)
library(ggpmisc)
library(scales)
library(gridExtra)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import Input_A counts and calcualte frequencing per column
In.A.counts <- read_excel("All_Counts.xlsx", sheet = "Input_A", col_names = T)
In.A.freq <- In.A.counts
In.A.freq[,-c(1:6)][is.na(In.A.freq[,-c(1:6)])] <- 0
In.A.freq[,-c(1:6)] <- lapply(In.A.freq[,-c(1:6)], function(x) x/sum(x))
In.A.freq$Sense.Mean <- rowMeans(In.A.freq[,c(7,9,11)], na.rm = T)
In.A.freq$Sense.SD <- apply(In.A.freq[,c(7,9,11)], 1, sd)
In.A.freq$Antisense.Mean <- rowMeans(In.A.freq[,c(8,10,12)], na.rm = T)
In.A.freq$Antisense.SD <- apply(In.A.freq[,c(8,10,12)], 1, sd)
In.A.freq.SS <- In.A.freq[,-c(7:12,15,16)] 
In.A.freq.SS <- plyr::rename(In.A.freq.SS, c("Sense.Mean" = "Mean", "Sense.SD" = "SD")) 
In.A.freq.SS$Strand <- "Sense"
In.A.freq.AS <- In.A.freq[,-c(7:14)] 
In.A.freq.AS <- plyr::rename(In.A.freq.AS, c("Antisense.Mean" = "Mean", "Antisense.SD" = "SD")) 
In.A.freq.AS$Strand <- "Antisense"
In.A.freq.strand <- rbind(In.A.freq.SS, In.A.freq.AS)
In.A.freq.strand$Strand <- factor(In.A.freq.strand$Strand, levels = unique(In.A.freq.strand$Strand))
In.A.freq.strand$Batch <- "A"

# import Input_B counts and calcualte frequencing per column
In.B.counts <- read_excel("All_Counts.xlsx", sheet = "Input_B", col_names = T)
In.B.freq <- In.B.counts
In.B.freq[,-c(1:6)][is.na(In.B.freq[,-c(1:6)])] <- 0
In.B.freq[,-c(1:6)] <- lapply(In.B.freq[,-c(1:6)], function(x) x/sum(x))
In.B.freq$Sense.Mean <- rowMeans(In.B.freq[,c(7,9,11,13)], na.rm = T)
In.B.freq$Sense.SD <- apply(In.B.freq[,c(7,9,11,13)], 1, sd)
In.B.freq$Antisense.Mean <- rowMeans(In.B.freq[,c(8,10,12,14)], na.rm = T)
In.B.freq$Antisense.SD <- apply(In.B.freq[,c(8,10,12,14)], 1, sd)
In.B.freq.SS <- In.B.freq[,-c(7:14,17,18)] 
In.B.freq.SS <- plyr::rename(In.B.freq.SS, c("Sense.Mean" = "Mean", "Sense.SD" = "SD")) 
In.B.freq.SS$Strand <- "Sense"
In.B.freq.AS <- In.B.freq[,-c(7:16)] 
In.B.freq.AS <- plyr::rename(In.B.freq.AS, c("Antisense.Mean" = "Mean", "Antisense.SD" = "SD")) 
In.B.freq.AS$Strand <- "Antisense"
In.B.freq.strand <- rbind(In.B.freq.SS, In.B.freq.AS)
In.B.freq.strand$Strand <- factor(In.B.freq.strand$Strand, levels = unique(In.B.freq.strand$Strand))
In.B.freq.strand$Batch <- "B"

# combine frequencies into single df
In.freq.strand <- rbind(In.A.freq.strand, In.B.freq.strand)
In.freq.strand$Strand <- factor(In.freq.strand$Strand, levels = unique(In.freq.strand$Strand))

# calculate Pearson correlation and p-value
cor.sense <- cor.test(In.freq.strand$Mean[In.freq.strand$Strand == "Sense" & In.freq.strand$Batch == "A"], 
                      In.freq.strand$Mean[In.freq.strand$Strand == "Sense" & In.freq.strand$Batch == "B"], 
                      method = "pearson")
cor.antisense <- cor.test(In.freq.strand$Mean[In.freq.strand$Strand == "Antisense" & In.freq.strand$Batch == "A"], 
                          In.freq.strand$Mean[In.freq.strand$Strand == "Antisense" & In.freq.strand$Batch == "B"], 
                          method = "pearson")

####################################
# Correlations graph
####################################
# arrange new df for correlation plots
In.freq.cor <- as.data.frame(In.A.freq.strand[,1:9])
In.freq.cor <- plyr::rename(In.freq.cor, c("Mean" = "Mean_A", "SD" = "SD_A"))
In.freq.cor <- setDT(In.freq.cor)[setDT(In.B.freq.strand[,c(1,7:9)]), on = c("Strand", "Location")]
In.freq.cor <- plyr::rename(In.freq.cor, c("Mean" = "Mean_B", "SD" = "SD_B"))
In.freq.cor$Strand <- factor(In.freq.cor$Strand, levels = unique(In.freq.cor$Strand))

# plot Sense separately from Antisense 
p.cor <- ggplot(data = In.freq.cor, aes(x = Mean_A, y = Mean_B, shape = Strand)) +
  geom_jitter(size = 3, stroke = 0.33, color = "black") +
  scale_shape_manual(values = c(21,23)) +
  facet_grid(Strand~., scales = "free") +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, aes(label = ..rr.label..), parse = T) +
  scale_x_continuous(name = "Mean Frequency, Replicate A", limits = c(NA, NA), labels = scientific) +
  scale_y_continuous(name = "Mean Frequency, Replicate B", limits = c(NA, NA), labels = scientific) +
  theme_bw() + labs(tag = "C") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 9, face = "bold", angle = 90),
        strip.background = element_rect(fill = "gray80"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black", angle = 90, hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 9, color = "black"),
        axis.title.y = element_text(face = "bold", size = 9, color = "black"))

####################################
# Box plot of nt frequencies according to Replicate, Frame, and Strand
####################################
In.freq.box <- In.freq.strand[In.freq.strand$Source == "CDS",]
In.freq.box$Frame <- substring(In.freq.box$VWF_aa, regexpr("\\.", In.freq.box$VWF_aa) + 1)
p.box <- ggplot()+
  geom_boxplot(data = In.freq.box[In.freq.box$Strand == "Sense",], aes(x=Batch, y=Mean, color=Frame, fill = Frame), outlier.shape = 21) +
  geom_boxplot(data = In.freq.box[In.freq.box$Strand == "Antisense",], aes(x=Batch, y=Mean, color=Frame, fill = Frame), outlier.shape = 23) +
  facet_grid(.~Strand, scales = "free") +
  scale_color_manual(values = c("steelblue4", "sienna2", "black")) +
  scale_fill_manual(values = alpha(c("steelblue4", "sienna2", "yellow1"), 0.4)) +
  theme_bw() + labs(tag = "B", x = "Unselected Phage Replicates", y = "Mean Frequency") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(10, "mm"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "gray80"),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(face = "bold", size = 9, color = "black"),
        axis.title.y = element_text(face = "bold", size = 9, color = "black")) +
  scale_y_continuous(trans = "log2", breaks = c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)) 

####################################
# Scatter plot of nt frequencies (mean +/- SD)
####################################
# for graphing log of frequencies (not correlations), turn 0 into 1e-9
In.freq.strand$Mean[In.freq.strand$Mean == 0] <- 1e-9

# Coordinantes for VWF CDS (in top portion of facet)
ann.rect.VWF <- data.frame(xmin = 732, xmax = 9170, ymin = 0.004, ymax = 0.006)
ann.text.VWF <- data.frame(x = 4950.5, y = 0.012, label = c("VWF"))
ann.rect.Ori <- data.frame(xmin = 9707, xmax = 10326, ymin = 0.004, ymax = 0.006)
ann.text.Ori <- data.frame(x = 9950, y = 0.012, label = c("ori"))
ann.rect.amp <- data.frame(xmin = 10481, xmax = 11341, ymin = 0.004, ymax = 0.006)
ann.text.amp <- data.frame(x = 11100, y = 0.012, label = c("AmpR"))

# graph average frequencies
Batch.labels <- c(A = "Unselected Phage, Replicate A", B = "Unselected Phage, Replicate B")
p.In <- ggplot() + 
  geom_errorbar(data = In.freq.strand[In.freq.strand$Batch == "A",], aes(x=Location, ymin=Mean-SD, ymax=Mean+SD), width = 100, size = 0.1, color = "black") +
  geom_errorbar(data = In.freq.strand[In.freq.strand$Batch == "B",], aes(x=Location, ymin=Mean-SD, ymax=Mean+SD), width = 100, size = 0.1, color = "black") +
  geom_jitter(data = In.freq.strand[In.freq.strand$Batch == "A",], aes(x=Location, y=Mean, shape=Strand, fill = Source), size = 1, stroke = 0.3, color = "black") + 
  geom_jitter(data = In.freq.strand[In.freq.strand$Batch == "B",], aes(x=Location, y=Mean, shape=Strand, fill = Source), size = 1, stroke = 0.3, color = "black") + 
  scale_shape_manual(values = c(21,23)) +
  scale_fill_manual(values = c("olivedrab1", "white"), guide = "none") +
  facet_grid(.~Batch, scales = "free_y", labeller = labeller(Batch = Batch.labels)) +
  labs(shape="Strand") +
  theme_bw() + labs(tag = "A") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0, "mm"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "gray80"),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(face = "bold", size = 9, color = "black"),
        axis.title.y = element_text(face = "bold", size = 9, color = "black")) +
  scale_x_continuous(name = "Nucleotide Position in VWF Plasmid", limits = c(NA,12000), 
                     breaks = c(1, 3000, 6000, 9000, 12000)) +
  scale_y_log10(name = "Mean Frequency", limits = c(NA, NA)) +
  geom_rect(data = ann.rect.VWF, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "olivedrab1", color = "black", size = 0.3) +
  geom_text(data = ann.text.VWF, aes(x = x, y = y, label = label), fontface = "bold", size = 2) +
  geom_rect(data = ann.rect.Ori, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black", size = 0.3) +
  geom_text(data = ann.text.Ori, aes(x = x, y = y, label = label), fontface = "bold", size = 2) +
  geom_rect(data = ann.rect.amp, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black", size = 0.3) +
  geom_text(data = ann.text.amp, aes(x = x, y = y, label = label), fontface = "bold", size = 2) 

####################################
# Combine plots
####################################
p.in.all <- grid.arrange(arrangeGrob(p.In, p.box, ncol = 1), arrangeGrob(p.cor), ncol=2, widths = 2:1)
ggsave("Inputs_Figure.tiff", plot = p.in.all, width = 190, height = 150, units = "mm", device = "tiff", compression = "lzw")

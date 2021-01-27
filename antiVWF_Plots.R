# Plot enrichment results
library(ggplot2)
library(readxl)
library(ggpubr)
library(data.table)
library(scales)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import counts and define samples
all.data <- lapply(excel_sheets("DESeq2 Results.xlsx"), read_excel, path = "DESeq2 Results.xlsx")
names(all.data) <- excel_sheets("DESeq2 Results.xlsx")  # rename df in list with sheet names
samples <- list("Neg", "anti-VWF", "I-1", "II-1", "II-2")
all.data <- Map(cbind, all.data, Sample = samples)

#######################################################
# Volcano as parts (due to large difference in y-scale)
#######################################################
# Volcano negative control
volcano.1 <- rbindlist(all.data[c("Neg")])
volcano.1 <- volcano.1[volcano.1$Source == "CDS"]
volcano.1$Strand <- factor(volcano.1$Strand, levels = unique(volcano.1$Strand))
volcano.1$Sample <- factor(volcano.1$Sample, levels = unique(volcano.1$Sample))
p.volcano.1 <- ggplot(volcano.1, aes(log2FoldChange, log.padj)) +
  geom_vline(xintercept = 0, color = "black") + geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_jitter(aes(shape = Strand, fill = Frame, color = Frame), size = 2, stroke = 1) +
  facet_grid(Strand~Sample, scales = "free", labeller = as_labeller(c("Neg" = "Negative"))) +
  scale_shape_manual(values = c(21,23)) +
  scale_color_manual(values = c("steelblue4", "sienna2", "black")) +
  scale_fill_manual(values = alpha(c("steelblue4", "sienna2", "yellow1"), 0.6)) +
  scale_y_continuous(name = expression(bold(paste("-log"["10"], "[q-value]")))) +
  theme_classic() + guides(shape = F, fill = guide_legend(override.aes = list(shape=21,size=3))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9, face = "bold", color = "black"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 9, face = "bold"),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "black", size = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5))

# Volcano polyclonal rabbit anti-VWF
volcano.2 <- rbindlist(all.data[c("VWFab")])
volcano.2 <- volcano.2[volcano.2$Source == "CDS"]
volcano.2$Strand <- factor(volcano.2$Strand, levels = unique(volcano.2$Strand))
volcano.2$Sample <- factor(volcano.2$Sample, levels = unique(volcano.2$Sample))
p.volcano.2 <- ggplot(volcano.2, aes(log2FoldChange, log.padj)) +
  geom_vline(xintercept = 0, color = "black") + geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_jitter(aes(shape = Strand, fill = Frame, color = Frame), size = 2, stroke = 1) +
  facet_grid(Strand~Sample, scales = "free") +
  scale_shape_manual(values = c(21,23)) +
  scale_color_manual(values = c("steelblue4", "sienna2", "black")) +
  scale_fill_manual(values = alpha(c("steelblue4", "sienna2", "yellow1"), 0.6)) +
  scale_y_continuous(name = expression(bold(paste("-log"["10"], "[q-value]")))) +
  theme_classic() + guides(shape = F, fill = guide_legend(override.aes = list(shape=21,size=3))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9, face = "bold", color = "black"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 9, face = "bold"),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "black", size = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5))

# Volcano anti-VWF alloantibodies
volcano.3 <- rbindlist(all.data[c("VWFi_I_1", "VWFi_II_1", "VWFi_II_2")])
volcano.3 <- volcano.3[volcano.3$Source == "CDS"]
volcano.3$Strand <- factor(volcano.3$Strand, levels = unique(volcano.3$Strand))
volcano.3$Sample <- factor(volcano.3$Sample, levels = unique(volcano.3$Sample))
p.volcano.3 <- ggplot(volcano.3, aes(log2FoldChange, log.padj)) +
  geom_vline(xintercept = 0, color = "black") + geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_jitter(aes(shape = Strand, fill = Frame, color = Frame), size = 2, stroke = 1) +
  facet_grid(Strand~Sample, scales = "free") +
  scale_shape_manual(values = c(21,23)) +
  scale_color_manual(values = c("steelblue4", "sienna2", "black")) +
  scale_fill_manual(values = alpha(c("steelblue4", "sienna2", "yellow1"), 0.6)) +
  scale_x_continuous(name = expression(bold(paste("log"["2"], "[FC] (Selected vs. Unselected)"))),
                     breaks = c(-4, -2, 0, 2, 4, 6, 8)) +
  scale_y_continuous(name = expression(bold(paste("-log"["10"], "[q-value]")))) +
  theme_classic() + guides(shape = F, fill = guide_legend(override.aes = list(shape=21,size=3))) +
  theme(legend.position = "none",
        legend.box.spacing = unit(0.01, "cm"),
        axis.text.x = element_text(size = 9, face = "bold", color = "black"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 9, face = "bold"),
        strip.text.y = element_text(size = 9, face = "bold", angle = 90),
        strip.background = element_rect(fill = "gray92", color = "black", size = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5))

# combine and save arranged plots
p.volcano <- ggarrange(p.volcano.1, p.volcano.2, p.volcano.3, nrow = 1, widths = c(1,1,3), common.legend = T, legend = "top")
p.volcano <- annotate_figure(p.volcano,
                             left = text_grob(expression(bold(paste("-log"["10"], "[q-value]"))),
                                              size = 9, face = "bold", rot = 90),
                             bottom = text_grob(expression(bold(paste("log"["2"], "[FC] (Selected vs. Unselected)"))),
                                                size = 9, face = "bold"))
ggsave("Volcano.tiff", plot = p.volcano, width = 173, height = 75, units = "mm", device = "tiff", compression = "lzw")

################################################################################################
################################################################################################
# Plot enriched N-term and C-term for all anti-VWF (rabbit polyclonal and VWF inhibitors)
################################################################################################
# Comibine all anti-VWF into single dataframe
VWF.plot <- rbindlist(all.data[c("VWFab", "VWFi_I_1", "VWFi_II_1", "VWFi_II_2")])
VWF.plot <- subset(VWF.plot, (Strand == "Sense" & Frame == 1) | (Strand == "Antisense" & Frame == 2))
VWF.plot$Strand <- factor(VWF.plot$Strand, levels = unique(VWF.plot$Strand))
VWF.plot$Location.aa <- as.numeric(substr(VWF.plot$VWF_aa, 2, nchar(VWF.plot$VWF_aa)-2))
VWF.plot$End[VWF.plot$Frame == 1] <- "N-term"
VWF.plot$End[VWF.plot$Frame == 2] <- "C-term"
VWF.plot$End <- factor(VWF.plot$End, levels = unique(VWF.plot$End))

###############################
# Plot 1 - annotations and maps
# Coordinantes for VWF domains (in top portion of facet)
ann.rect.VWF <- data.frame(xmin = c(23, 388, 764, 866, 1238, 1474, 1681, 1875, 2256, 2721), 
                           xmax = c(387, 763, 865, 1237, 1473, 1680, 1874, 2255, 2720, 2813), 
                           ymin = rep(8.9, length.out = 10), ymax = rep(10.6, length.out = 10))
ann.text.VWF <- data.frame(x = c(193.5, 575, 814, 1051, 1360, 1581, 1782, 2066, 2488, 2767), 
                           y = rep(9.75, length.out = 10), 
                           label = c("D1", "D2", "D'", "D3", "A1", "A2", "A3", "D4", "C1-C6", "CK"))

# Coordinates for maximum and minimum length of displayed VWF fragment
VWF.size.line <- data.frame(x = c(1, 100, 1, 333), y = c(12, 12, 13, 13), 
                            End = factor(c("N-term", "C-term"), levels = unique(VWF.plot$End)),
                            Frag.size = c("Min", "Min", "Max", "Max"))
VWF.size.text <- data.frame(x = c(200, 450), y = c(12, 13), 
                            label = c("Min", "Max"))

# Coordinates for proteolytic VWF fragments
Frag.size.rect <- data.frame(xmin = c(764, 1674, 2129, 764, 1243, 1212, 1036, 1437), 
                             xmax = c(2128, 2128, 2813, 1035, 1481, 1491, 1274, 1491),
                             ymin = c(11, 12.5, 11, 12.5, 12.5, 14, 15.5, 15.5),
                             ymax = c(12.4, 13.9, 12.4, 13.9, 13.9, 15.4, 16.9, 16.9))
Frag.size.text <- data.frame(x = c(1445, 1900, 2471, 899, 1362, 1351, 1550), 
                             y = c(11.75, 13.25, 11.75, 13.25, 13.25, 14.75, 16.25), 
                             label = c("f3", "f1", "f2", "f4", "f5", "f6", "f7"))
Frag.size.line <- data.frame(x = c(1274, 1437), y = c(16.25, 16.25))

p1 <- ggplot()+ 
  geom_rect(data = ann.rect.VWF, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "olivedrab1", color = "black", alpha = 0.4) +
  geom_text(data = ann.text.VWF, aes(x = x, y = y, label = label), fontface = "bold", size = 2.5) +
  geom_line(data = VWF.size.line[VWF.size.line$Frag.size == "Min",], aes(x = x, y = y), linetype = "solid", color = "black") +
  geom_line(data = VWF.size.line[VWF.size.line$Frag.size == "Max",], aes(x = x, y = y), linetype = "solid", color = "black") +
  geom_point(data = VWF.size.line, aes(x = x, y = y, shape = End), fill = "white", size = 1.5, stroke = 0.7, color = "black", show.legend = F) +
  scale_shape_manual(values = c(21,23)) + 
  geom_text(data = VWF.size.text, aes(x = x, y = y, label = label), fontface = "bold", size = 2.5) +
  geom_rect(data = Frag.size.rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black") +
  geom_text(data = Frag.size.text, aes(x = x, y = y, label = label), fontface = "bold", size = 2.5) +
  geom_line(data = Frag.size.line, aes(x = x, y = y), linetype = "solid", color = "black") +
  scale_x_continuous(limits = c(NA,2850)) +
  theme_void() + theme(legend.position = "none")

###############################
# Plot 2 - polyclonal rabbit anti-VWF
# Coordinates for VWF epitopes of Tan et al. (2007)
Tan.rect.VWF <- data.frame(xmin = c(809, 1262, 1788, 2285, 2403, 2462, 2677, 2701), 
                           xmax = c(826, 1313, 1865, 2319, 2447, 2498, 2686, 2729), 
                           ymin = rep(9, length.out = 8), ymax = rep(9.5, length.out = 8))
Tan.text.VWF <- data.frame(x = c(817.5, 1287.5, 1826.5, 2275, 2405, 2510, 2650, 2755), 
                           y = rep(10.2, length.out = 8), 
                           label = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8"))

p2 <- ggplot()+ 
  geom_rect(data = Tan.rect.VWF, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", color = "black") +
  geom_text(data = Tan.text.VWF, aes(x = x, y = y, label = label), fontface = "bold", size = 2.5) + 
  geom_jitter(data = VWF.plot[VWF.plot$Sample == "anti-VWF"], aes(x=Location.aa, y=log2FoldChange, fill=log.padj, shape=End), 
              size = 1.5, stroke = 0.33, color = "black") + 
  scale_shape_manual(values = c(21,23)) + 
  scale_fill_gradientn(name=expression(atop(paste("-log"["10"],"      "),"(q-value)")), 
                       colors = c("steelblue4", "gold"), na.value="white",
                       values = rescale(c(min(VWF.plot$log.padj, na.rm = T), 30, 
                                          max(VWF.plot$log.padj, na.rm = T))),
                       breaks = c(1,2,5,10,20,40,80), trans = "log", limits=c(1,80), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  geom_hline(yintercept = log(1.5, base = 2), color = "black", linetype = "dashed") +
  scale_x_continuous(name = "VWF Position (aa)", limits = c(NA,2850), breaks = c(1, 500, 1000, 1500, 2000, 2500, 2813)) +
  labs(y = expression(bold(paste("log"["2"],"[Fold Change (Selected/Unselected)]"))), 
       shape="VWF End") +
  geom_text(aes(x = 1, y = 9, label = "anti-VWF", hjust = "left", fontface = "bold")) +
  theme_bw() + guides(shape = guide_legend(order = 1, keyheight = unit(0, "cm"))) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

###############################
# Plot 3 - I-1
# Coordinantes for patient I-1 deletions
ann.tri.I <- data.frame(x = 812, y = 6)
ann.line.I <- data.frame(x = c(730, 814), y = c(7, 7))
ann.text.I <- data.frame(x = c(840, 730), y = c(6.2, 8), label = c("c.2435del", "g.Ex17_Ex18del"))

p3 <- ggplot()+ 
  geom_point(data = ann.tri.I, aes(x = x, y = y), shape = 17, color = "firebrick") +
  geom_line(data = ann.line.I, aes(x = x, y = y), linetype = "solid", color = "firebrick") +
  geom_text(data = ann.text.I, aes(x = x, y = y, label = label), hjust = 0, fontface = "bold", size = 2.5) +
  geom_jitter(data = VWF.plot[VWF.plot$Sample == "I-1"], aes(x=Location.aa, y=log2FoldChange, fill=log.padj, shape=End), 
              size = 1.5, stroke = 0.33, color = "black") + 
  scale_shape_manual(values = c(21,23)) + 
  scale_fill_gradientn(name=expression(atop(paste("-log"["10"],"      "),"(q-value)")), 
                       colors = c("steelblue4", "gold"), na.value="white",
                       values = rescale(c(min(VWF.plot$log.padj, na.rm = T), 30, 
                                          max(VWF.plot$log.padj, na.rm = T))),
                       breaks = c(1,2,5,10,20,40,80), trans = "log", limits=c(1,80), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  geom_hline(yintercept = log(1.5, base = 2), color = "black", linetype = "dashed") +
  scale_x_continuous(name = "VWF Position (aa)", limits = c(NA,2850), breaks = c(1, 500, 1000, 1500, 2000, 2500, 2813)) +
  labs(y = expression(bold(paste("log"["2"],"[Fold Change (Selected/Unselected)]"))), 
       shape="VWF End") +
  geom_text(aes(x = 1, y = 6, label = "I-1", hjust = "left", fontface = "bold")) +
  theme_bw() + guides(shape = guide_legend(order = 1, keyheight = unit(0, "cm"))) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.spacing.y = unit(0.1,"cm"),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

###############################
# Plot 4 - II-1
# Coordinantes for patients II-1 and II-2 genetic deletions
ann.line.II <- data.frame(x = rep(c(220, 2629), length.out =2), y = rep(9.5, length.out = 4), 
                          Sample = factor(c("II-1", "II-1", "II-2", "II-2"), levels = unique(VWF.plot$Sample)))
ann.text.II <- data.frame(x = 1424.5, y = 10.5, label = "c.658_7887del",
                          Sample = factor(c("II-1", "II-1", "II-2", "II-2"), levels = unique(VWF.plot$Sample)))

p4 <- ggplot()+ 
  geom_line(data = ann.line.II, aes(x = x, y = y), linetype = "solid", color = "firebrick") +
  geom_text(data = ann.text.II, aes(x = x, y = y, label = label), fontface = "bold", size = 2.5) +
  geom_jitter(data = VWF.plot[VWF.plot$Sample == "II-1"], aes(x=Location.aa, y=log2FoldChange, fill=log.padj, shape=End), 
              size = 1.5, stroke = 0.33, color = "black") + 
  scale_shape_manual(values = c(21,23)) + 
  scale_fill_gradientn(name=expression(atop(paste("-log"["10"],"      "),"(q-value)")), 
                       colors = c("steelblue4", "gold"), na.value="white",
                       values = rescale(c(min(VWF.plot$log.padj, na.rm = T), 30, 
                                          max(VWF.plot$log.padj, na.rm = T))),
                       breaks = c(1,2,5,10,20,40,80), trans = "log", limits=c(1,80), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  geom_hline(yintercept = log(1.5, base = 2), color = "black", linetype = "dashed") +
  scale_x_continuous(name = "VWF Position (aa)", limits = c(NA,2850), breaks = c(1, 500, 1000, 1500, 2000, 2500, 2813)) +
  labs(y = expression(bold(paste("log"["2"],"[Fold Change (Selected/Unselected)]"))), 
       shape="VWF End") +
  geom_text(aes(x = 1, y = 8, label = "II-1", hjust = "left", fontface = "bold")) +
  theme_bw() + guides(shape = guide_legend(order = 1, keyheight = unit(0, "cm"))) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

###############################
# Plot 5 - II-2
p5 <- ggplot()+ 
  geom_line(data = ann.line.II, aes(x = x, y = y), linetype = "solid", color = "firebrick") +
  geom_text(data = ann.text.II, aes(x = x, y = y, label = label), fontface = "bold", size = 2.5) +
  geom_jitter(data = VWF.plot[VWF.plot$Sample == "II-2"], aes(x=Location.aa, y=log2FoldChange, fill=log.padj, shape=End), 
              size = 1.5, stroke = 0.33, color = "black") + 
  scale_shape_manual(values = c(21,23)) + 
  scale_fill_gradientn(name=expression(atop(paste("-log"["10"],"      "),"(q-value)")), 
                       colors = c("steelblue4", "gold"), na.value="white",
                       values = rescale(c(min(VWF.plot$log.padj, na.rm = T), 30, 
                                          max(VWF.plot$log.padj, na.rm = T))),
                       breaks = c(1,2,5,10,20,40,80), trans = "log", limits=c(1,80), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  geom_hline(yintercept = log(1.5, base = 2), color = "black", linetype = "dashed") +
  scale_x_continuous(name = "VWF Position (aa)", limits = c(NA,2850), breaks = c(1, 500, 1000, 1500, 2000, 2500, 2813)) +
  labs(y = expression(bold(paste("log"["2"],"[Fold Change (Selected/Unselected)]"))), 
       shape="VWF End") +
  geom_text(aes(x = 1, y = 8, label = "II-2", hjust = "left", fontface = "bold")) +
  theme_bw() + guides(shape = guide_legend(order = 1, keyheight = unit(0, "cm"))) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Combine and save enrichment plots
p.VWF <- egg::ggarrange(p1, p2, p3, p4, p5, ncol = 1, heights = c(0.5, 1, 1, 1, 1))
p.VWF <- annotate_figure(p.VWF, left = text_grob(expression(bold(paste("log"["2"],"[Fold Change (Selected/Unselected)]"))), rot = 90),
                         bottom = text_grob("VWF Position (aa)", face = "bold", hjust = 0.7))
ggsave("Enrichments.tiff", plot = p.VWF, width = 173, height = 160, units = "mm", device = "tiff", compression = "lzw")


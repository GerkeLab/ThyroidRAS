library(here)
library(tidyverse)
library(ggsci)
library(gridExtra)
library(knitr)
library(ggdendro)
library(reshape2)
library(grid)
library(cowplot)
library(limma)

# import data
load(here("data","Full_data.RData"))

load(here("data","Results_data.RData"))

# subset clinical into two data frame with patients in the same order as they appear in 
# the gene expression data 
clinicalio <- clinical[clinical$PATIENT_ID %in% io360$PATIENT_ID,]
clinicalio <- clinicalio[match(io360$PATIENT_ID,clinicalio$PATIENT_ID),]

clinicalp <- clinical[clinical$PATIENT_ID %in% pathways$PATIENT_ID,]
clinicalp <- clinicalp[match(pathways$PATIENT_ID,clinicalp$PATIENT_ID),]

#####-----------------------------------------------------------------------------------------#####
#                                             FIGURE 1
#####-----------------------------------------------------------------------------------------#####

# make the data into long format for plotting 
mutations_long <- melt(clinical[c("PATIENT_ID","KRAS","HRAS","NRAS","RET",
                                  "TERT","EIF1AX")],
                       id="PATIENT_ID")

# fix annotation names according to comments from Dr.Chung/Prera
mutations_long <- mutations_long %>%
  mutate(value = ifelse(value=="PTC1","RET-PTC1",value)) %>%
  mutate(value = ifelse(value=="A113 splice","A113fs/splice", value)) %>%
  mutate(value = ifelse(is.na(value),"",value)) %>%
  mutate(value = case_when(
    PATIENT_ID == "ID-65" & variable %in% c("RET","TERT","EIF1AX") ~ "Unmeasured",
    TRUE ~ value 
  ))

# need to refactor and add additional details inorder to get faux ~split~ legend 
mutations_long$value = factor(mutations_long$value,
                              levels=c("RAS mutations","G12C","G12D","G12V",
                                       "G13D","Q61K","Q61R"," ","Other alterations",
                                       "C250T","338TT","A113fs/splice","RET-PTC1",
                                        "", "Unmeasured"))

# force ordering of x axis to Juan's preference
mutations_long$variable = factor(mutations_long$variable,
                                 levels=c("KRAS","HRAS","NRAS",
                                          "TERT","EIF1AX","RET"))

# plot mutations 
ggplot(mutations_long, aes(x=variable,y=PATIENT_ID)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_manual(values = c("white","#DF8F44FF","#00A1D5FF","#B24745FF","#79AF97FF","#6A6599FF",
                               "#925E9FFF","white","white","#0073C2FF","#EFC000FF",
                               "#8F7700FF","#CD534CFF","white","grey50"),
                    labels = c("RAS mutations","G12C","G12D","G12V",
                               "G13D","Q61K","Q61R"," ","Other alterations",
                               "C250T","338TT","A113fs/splice","RET-PTC1",
                               "", "Unmeasured","  "),
                    drop=FALSE) + 
  xlab("Gene") + ylab("Patient ID") + 
  labs(fill="") + 
  theme_classic()

rm(mutations_long)


#####-----------------------------------------------------------------------------------------#####
#                                             FIGURE 2A
#####-----------------------------------------------------------------------------------------#####
# unsupervised clustering of pathways panel w/ annotations for :
# ras mutation, bethesda, histology, nuclear atypia and age

path_sig_genes <- pathResults[round(pathResults$pval,2)<0.05,]$gene
path_not_sig_genes <- pathResults[round(pathResults$pval,2)>=0.05,]$gene

# center and scale 
pathwaysc <- scale(pathways[,2:ncol(pathways)])
row.names(pathwaysc) <- pathways$PATIENT_ID

# dendogram / unsupervised clustering - samples 
# based on all samples 
path.dendro <- as.dendrogram(hclust(d=dist(x=pathwaysc))) # get order
path.dendro.plot <- ggdendrogram(data = path.dendro, rotate = TRUE) # plot for clustering

# dendogram / unsupervised clustering - sig genes 
path.dendro.genes1 <- as.dendrogram(hclust(d=dist(x=t(pathwaysc[,path_sig_genes])))) # get order
path.dendro.genes.plot1 <- ggdendrogram(data = path.dendro.genes1, rotate = TRUE) # plot for clustering

# dendogram / unsupervised clustering - not sig genes 
path.dendro.genes2 <- as.dendrogram(hclust(d=dist(x=t(pathwaysc[,path_not_sig_genes])))) # get order
path.dendro.genes.plot2 <- ggdendrogram(data = path.dendro.genes2, rotate = TRUE) # plot for clustering

## make gene expression data long for plotting in heatmap 
path_long <- melt(pathwaysc,id="Row.names")
## id as factor to order matching dendrogram 
path_long$Var1 <- as.character(path_long$Var1)
path_long$Var1 <- factor(x = path_long$Var1,
                         levels = c(labels(path.dendro)),
                         ordered = TRUE)
## genes as factor to order matching dendrogram 
path_long$Var2 <- as.character(path_long$Var2)
path_long$Var2 <- factor(x = path_long$Var2,
                         levels = c(labels(path.dendro.genes1), labels(path.dendro.genes2)),
                         ordered = TRUE)

# clinical annotation
clinicalp$PATIENT_ID <- factor(x=clinicalp$PATIENT_ID, # id as factor to order matching dendrogram 
                             levels = c(labels(path.dendro)),
                             ordered = TRUE)

# main heatmap showing gene expression
gene_heatmap <- ggplot(data = path_long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="#0095c6",mid="white",high="#dd8a3b") +
  labs(y="Sample") + 
  guides(fill=FALSE) + 
  geom_vline(xintercept = length(labels(path.dendro.genes1))+0.5) + 
  theme(#legend.position = "top",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) 

# subset to data needed for additional heatmap annotations 
clinical_long <- melt(clinicalp[,c("PATIENT_ID","histology","ras_mutation",
                                   "Beth_cat","nuclear_atypia","tumor_size")],
                      id="PATIENT_ID")

ras_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="ras_mutation",],
                      aes(x = variable, y = `PATIENT_ID`)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#374E55FF","#B24745FF","#80796BFF")) + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

hist_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="histology",],
                       aes(x = variable, y = `PATIENT_ID`)) +
  geom_tile(aes(fill = value)) +
  # scale_fill_ucscgb() + # not enough colors in jama
  # scale_fill_hue(l=40, c=45) + 
  scale_fill_jama() + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

beth_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="Beth_cat",],
                       aes(x = variable, y = `PATIENT_ID`)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#79AF97FF","#B24745FF","#6A6599FF","#374E55FF","#80796BFF")) +
  # scale_fill_jama() +
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

nuclear_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="nuclear_atypia",],
                          aes(x = variable, y = `PATIENT_ID`)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#374E55FF", "#DF8F44FF","#00A1D5FF")) + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

size_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="tumor_size",],
                       aes(x = variable, y = `PATIENT_ID`)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#79AF97FF","#6A6599FF","#80796BFF")) +
  # scale_fill_jama() + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

# pdf(file="/Volumes/Lab_Gerke/thyroidInnovation/figures_for_manuscript/Figure2A_pathways_heatmap.pdf",
#     width=6.88, height=7.88)
plot_grid(gene_heatmap, ras_heatmap, beth_heatmap, hist_heatmap,
          nuclear_heatmap, size_heatmap, 
          align = "h", nrow = 1, 
          rel_widths = c(7.5/10,0.5/10,0.5/10,0.5/10,0.5/10,0.5/10),
          labels = c("","Ras","Beth.","Hist.","Atypia","Size"),
          label_size = 6,
          label_x = c(0.1,-0.01,-0.3,-0.25,-0.35,-0.2))
# dev.off()

rm(path.dendro, path.dendro.plot, gene_heatmap, beth_heatmap, hist_heatmap,
   nuclear_heatmap, size_heatmap, ras_heatmap, path.dendro.genes1, path.dendro.genes.plot1,
   clinical_long, path.dendro.genes2, path.dendro.genes.plot2)

#####-----------------------------------------------------------------------------------------#####
#                                             FIGURE 2B
#####-----------------------------------------------------------------------------------------#####
# unsupervised clustering of io360 panel w/ annotations for :
# ras mutation, bathesda, histology, nuclear atypia and tumor size

io_sig_genes <- ioResults[round(ioResults$pval,2)<0.05,]$gene
io_not_sig_genes <- ioResults[round(ioResults$pval,2)>=0.05,]$gene

# center and scale 
io360c <- scale(io360[,2:ncol(io360)])
row.names(io360c) <- io360$PATIENT_ID

# dendogram / unsupervised clustering
io360.dendro <- as.dendrogram(hclust(d=dist(x=io360c))) # get order
io360.dendro.plot <- ggdendrogram(data = io360.dendro, rotate = TRUE) # plot for clustering

# dendogram / unsupervised clustering - sig genes 
io360.dendro.genes1 <- as.dendrogram(hclust(d=dist(x=t(io360c[,io_sig_genes])))) # get order
io360.dendro.genes.plot1 <- ggdendrogram(data = io360.dendro.genes1, rotate = TRUE) # plot for clustering
# dendogram / unsupervised clustering - not sig genes 
io360.dendro.genes2 <- as.dendrogram(hclust(d=dist(x=t(io360c[,io_not_sig_genes])))) # get order
io360.dendro.genes.plot2 <- ggdendrogram(data = io360.dendro.genes2, rotate = TRUE) # plot for clustering


# heatmap 
## make data long
io360_long <- melt(io360c,id="Row.names")
## id as factor to order matching dendrogram 
io360_long$Var1 <- as.character(io360_long$Var1)
io360_long$Var1 <- factor(x = io360_long$Var1,
                          levels = c(labels(io360.dendro)),
                          ordered = TRUE)
## gene as factor to order matching dendrogram 
io360_long$Var2 <- as.character(io360_long$Var2)
io360_long$Var2 <- factor(x = io360_long$Var2,
                          levels = c(labels(io360.dendro.genes1),labels(io360.dendro.genes2)),
                          ordered = TRUE)

# clinical annotation
clinicalio$PATIENT_ID <- factor(x=clinicalio$PATIENT_ID, # id as factor to order matching dendrogram 
                              levels = c(labels(io360.dendro)),
                              ordered = TRUE)

# main heatmap showing gene expression
gene_heatmap <- ggplot(data = io360_long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="#0095c6",mid="white",high="#dd8a3b") +
  labs(y="Sample") + 
  guides(fill=FALSE) +
  geom_vline(xintercept = length(labels(io360.dendro.genes1))+0.5) + 
  theme(#legend.position = "top",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) 


clinical_long <- melt(clinicalio[,c("PATIENT_ID","histology","ras_mutation",
                                    "Beth_cat","nuclear_atypia","tumor_size")],
                      id="PATIENT_ID")

ras_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="ras_mutation",],
                      aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#374E55FF","#B24745FF","#80796BFF")) + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

hist_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="histology",],
                       aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_jama() +
  # scale_fill_hue(l=40, c=45) + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

beth_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="Beth_cat",],
                       aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#79AF97FF","#B24745FF","#6A6599FF","#374E55FF","#80796BFF")) +
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

nuclear_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="nuclear_atypia",],
                          aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#374E55FF", "#DF8F44FF","#00A1D5FF")) + 
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

size_heatmap <- ggplot(data = clinical_long[clinical_long$variable=="tumor_size",],
                       aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#79AF97FF","#6A6599FF","#80796BFF")) +
  guides(fill=FALSE) + 
  theme(axis.line = element_blank(),
        plot.margin=margin(l=-0.35,unit="cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

# pdf(file="/Volumes/Lab_Gerke/thyroidInnovation/figures_for_manuscript/Figure2B_io360_heatmap.pdf",
#     width=6.88, height=7.88)
plot_grid(gene_heatmap, ras_heatmap, beth_heatmap, hist_heatmap,
          nuclear_heatmap, size_heatmap, 
          align = "h", nrow = 1, 
          rel_widths = c(7.5/10,0.5/10,0.5/10,0.5/10,0.5/10,0.5/10),
          labels = c("","Ras","Beth.","Hist.","Atypia","Size"),
          label_size = 6,
          label_x = c(0.1,-0.01,-0.3,-0.25,-0.35,-0.2))
# dev.off()

rm(io360.dendro, io360.dendro.plot, gene_heatmap, beth_heatmap, hist_heatmap,
   nuclear_heatmap, size_heatmap, ras_heatmap, io360.dendro.genes1, io360.dendro.genes.plot1,
   io360.dendro.genes2, io360.dendro.genes.plot2)

#####-----------------------------------------------------------------------------------------#####
#                                           FIGURE 2 LEGEND
#####-----------------------------------------------------------------------------------------#####

# legend plot 
# taken from: https://stackoverflow.com/questions/13143894/how-do-i-position-two-legends-independently-in-ggplot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend) }

l1 <- ggplot(data = path_long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="#0095c6",mid="white",high="#dd8a3b") + 
  theme(legend.position="bottom") + 
  labs(fill="Expression")

l2 <- ggplot(data = clinical_long[clinical_long$variable=="ras_mutation",],
             aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#374E55FF","#B24745FF","#80796BFF")) + 
  theme(legend.position="bottom") + 
  labs(fill="Ras")

l3 <-  ggplot(data = clinical_long[clinical_long$variable=="Beth_cat",],
              aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#79AF97FF","#B24745FF","#6A6599FF","#374E55FF","#80796BFF")) +
  theme(legend.position="bottom") + 
  labs(fill="Bethesda")

l4 <- ggplot(data = clinical_long[clinical_long$variable=="histology",],
             aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = str_wrap(value,width = 40))) +
  scale_fill_jama(labels=c("No Neoplastic Appearance","Benign","Indeterminant",
                           "Low-risk Malignant","High-risk Malignant")) + 
  labs(fill="Histology") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=2))
# theme(legend.key.height=unit(0.5, "cm"), legend.key.width=unit(0.5, "cm"))

l5 <- ggplot(data = clinical_long[clinical_long$variable=="nuclear_atypia",],
             aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#374E55FF", "#DF8F44FF","#00A1D5FF")) + 
  theme(legend.position="bottom") + 
  labs(fill="Nuclear Atypia")

l6 <- ggplot(data = clinical_long[clinical_long$variable=="tumor_size",],
             aes(x = variable, y = PATIENT_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values=c("#79AF97FF","#6A6599FF","#80796BFF"),
                    labels=c("0-2","2-4",">4")) +
  theme(legend.position="bottom") + 
  labs(fill="Tumor Size")

l1 <- g_legend(l1)
l2 <- g_legend(l2)
l3 <- g_legend(l3)
l4 <- g_legend(l4)
l5 <- g_legend(l5)
l6 <- g_legend(l6)

legend1_grid <- cowplot::plot_grid(l1, align = "v", nrow = 2)

legends <- legend1_grid + 
  ggplot2::annotation_custom(grob = l2, 
                             xmin = 0.225, xmax = 0.225, ymin = 0.67, ymax = 0.67) + 
  ggplot2::annotation_custom(grob = l3, 
                             xmin = 0.25, xmax = 0.25, ymin = 0.6, ymax = 0.6) + 
  ggplot2::annotation_custom(grob = l4, 
                             xmin = 0.377, xmax = 0.377, ymin = 0.5, ymax = 0.5) + 
  ggplot2::annotation_custom(grob = l5, 
                             xmin = 0.285, xmax = 0.285, ymin = 0.4, ymax = 0.4) + 
  ggplot2::annotation_custom(grob = l6, 
                             xmin = 0.225, xmax = 0.225, ymin = 0.33, ymax = 0.33)

legendsV2 <- legend1_grid + 
  ggplot2::annotation_custom(grob = l2, 
                             xmin = 0.6, xmax = 0.6, ymin = 0.75, ymax = 0.75) + 
  ggplot2::annotation_custom(grob = l3, 
                             xmin = 0.25, xmax = 0.25, ymin = 0.68, ymax = 0.68) + 
  ggplot2::annotation_custom(grob = l4, 
                             xmin = 0.38, xmax = 0.38, ymin = 0.6, ymax = 0.6) + 
  ggplot2::annotation_custom(grob = l5, 
                             xmin = 0.285, xmax = 0.285, ymin = 0.52, ymax = 0.52) + 
  ggplot2::annotation_custom(grob = l6, 
                             xmin = 0.78, xmax = 0.78, ymin = 0.52, ymax = 0.52)

# pdf(file="/Volumes/Lab_Gerke/thyroidInnovation/figures_for_manuscript/heatmap_legend.pdf",
#     width=6.88, height=7.88)
legends
# dev.off()

rm(l1,l2,l3,l4,l5,l6,legend1_grid)

#####-----------------------------------------------------------------------------------------#####
#                                             FIGURE 3A
#####-----------------------------------------------------------------------------------------#####

path_diff_genes <- pathResults[round(pathResults$adj_pval,2)<0.20,]$gene

path_tmp <- as.data.frame(cbind(status=clinicalp$histology,pathways))
melted_path <- melt(path_tmp[,c("PATIENT_ID","status",path_diff_genes)],
                    id.vars = c("PATIENT_ID","status"))

ggplot(melted_path[melted_path$variable %in% path_diff_genes,], aes(x=variable, y=value, fill=as.factor(status)))+
  scale_fill_manual(values=rep("white",length(path_diff_genes)+1))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(status), col=as.factor(status)),alpha=1)+
  geom_boxplot(outlier.shape=NA,alpha=0)+
  theme(axis.text.x = element_text(color = "black", size=10),#, angle=90, hjust = 1),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="Gene",y="Expression") +
  scale_color_jama()

rm(path_diff_genes,path_tmp,melted_path)

#####-----------------------------------------------------------------------------------------#####
#                                             FIGURE 3B
#####-----------------------------------------------------------------------------------------#####

io_diff_genes <- ioResults[round(ioResults$adj_pval,2)<0.20,]$gene

io_tmp <- as.data.frame(cbind(status=clinicalio$histology,io360))
melted_io <- melt(io_tmp[,c("PATIENT_ID","status",io_diff_genes)],
                  id.vars = c("PATIENT_ID","status"))

ggplot(melted_io[melted_io$variable %in% io_diff_genes,], aes(x=variable, y=value, fill=as.factor(status)))+
  scale_fill_manual(values=rep("white",length(io_diff_genes)+1))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(status), col=as.factor(status)),alpha=1)+
  geom_boxplot(outlier.shape=NA,alpha=0)+
  theme(axis.text.x = element_text(color = "black", size=10),#, angle=90, hjust = 1),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="Gene",y="Expression") +
  scale_color_jama()

rm(io_diff_genes,io_tmp,melted_io)

##### --------------------------------------------------------------------------------------- #####
#                                             FIGURE 4
##### --------------------------------------------------------------------------------------- #####

clinical %>%
  mutate(hscore = as.numeric(hscore)) %>%
  ggplot(aes(x = as.factor(histology), y = hscore, fill = as.factor(histology))) + 
  scale_fill_manual(values = rep("white", 5)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.75,
                                              jitter.height =  0,
                                              dodge.width = 0.75),
              aes(fill = as.factor(histology), col = as.factor(histology)),
              size = 3) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + 
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title = element_text(face = "bold", size = 15)) + 
  guides(fill = FALSE, colour = FALSE) + 
  labs(x = "", y = "H-score") +
  ggtitle("Angiopoeitin-2 IHC staining in \n RAS-mutant thyroid tumors") + 
  scale_color_jama() + 
  scale_x_discrete(labels = c("No Neoplastic \n Appearance", "Benign", "Indeterminant",
                              "Low-risk \n Malignant", "High-risk \n Malignant"))
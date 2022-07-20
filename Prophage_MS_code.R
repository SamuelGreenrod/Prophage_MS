# MS R script

library("ape")
library("Biostrings")
library("pgirmess")
library("ggplot2")
library("ggtree")
library("ggnewscale")
library("magrittr")
library("ggpubr")
library("pheatmap")
library("broom")
library("datarium")
library("rstatix")
library("tidyverse")
library("treeio")
source("functions/summarySE.R")
library("ggbeeswarm")
library("reshape2")
library(gridExtra)
library(nlme)
library("biobase")
library(piecewiseSEM)
library("scatterpie")
library("maps")
packageVersion("maps")
library("paco")
library("phytools")
library("RColorBrewer")
library("FSA")
library(vegan)
library("piecewiseSEM")
library(phangorn)
library(MASS)
library(dplyr)
library(tidyr)
library(patchwork)


### Revisions ###

# Shared prophage hit
library("ggVennDiagram")

Shared <- read.csv("PHASTER_vs_PHISPY_vs_virsorter2CheckV_venndiagram.csv",fileEncoding="UTF-8-BOM")
PHASTER <- Shared$PHASTER_list
PhiSpy <- Shared$Phispy_list
Virsorter <- Shared$Virsorter_checkV_list

x <- list(
  "PHASTER" = PHASTER, 
  "PhiSpy" = PhiSpy, 
  "Virsorter" = Virsorter
)


ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 1, set_name_size = 7,text_size = 5
)


VennPlotData(PHASTER)
setClass("PHASTER",PHASTER)

A <- ggVennDiagram(list(PHASTER, PhiSpy, Virsorter),category.names = c("PHASTER","PhiSpy","Virsorter2"),set_size = 5)+
  ggplot2::scale_fill_gradient(low="grey95",high = "firebrick3")+
  labs(fill="Prophage count")+ theme(legend.title=element_text(size=14))
geom_sf(size = 1.2, color = "black", data = venn_setedge(data))

ggVennDiagram(
  list(PHASTER, PhiSpy, Virsorter), label_alpha = 0.7,
  category.names = c("PHASTER","PhiSpy","Virsorter2"),set_size = 5)+
  ggplot2::scale_fill_gradient(low="grey95",high = "firebrick3")+
  geom_sf(size = 1.2, color = "black", data = venn_setedge(data))+
  labs(fill="Prophage count")+ theme(legend.title=element_text(size=14))

Shared <- read.csv("PHASTER_vs_PHISPY_vs_virsorter2CheckV_intactvenndiagram.csv",fileEncoding="UTF-8-BOM")
PHASTER <- Shared$PHASTER_list2
PhiSpy <- Shared$Phispy_list2
Virsorter <- Shared$Virsorter_checkV_list2

data2 <- list(PHASTER, PhiSpy, Virsorter)
ggVennDiagram(
  data2, label_alpha = 0.7,
  category.names = c("PHASTER","PhiSpy","Virsorter2"),set_size = 5)+
  ggplot2::scale_fill_gradient(low="grey95",high = "navyblue")+
  geom_sf(size = 1.2, color = "black", data = venn_setedge(data))+
  labs(fill="Prophage count")+ theme(legend.title=element_text(size=14))

ggVennDiagram(list(PHASTER, PhiSpy, Virsorter),category.names = c("PHASTER","PhiSpy","Virsorter2"),set_size = 5)+
  ggplot2::scale_fill_gradient(low="grey95",high = "grey40")+
  labs(fill="Intact prophage count")+ theme(legend.title=element_text(size=14))
geom_sf(size = 1.2, color = "black", data = venn_setedge(data))

Shared <- read.csv("Intact_prophages_annoted_vs_toolcomparison_venndiagram.csv",fileEncoding="UTF-8-BOM")
PHASTER_only <- Shared$Only_PHASTER
More_than_one_tool <- Shared$More_than_2_tools
Annoted <- Shared$Annotated

x <- list(
  "Only found with PHASTER" = PHASTER_only, 
  "Found with >1 tool" = More_than_one_tool, 
  "Taxonomic annotation" = Annoted
)


ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 1, set_name_size = 7,text_size = 5
)


### Cornerstone phage gene content

genes <- read.csv("Intact_prophage_cornerstone_phage_genecontent.csv",fileEncoding="UTF-8-BOM")
genes_onlycategory <- select(genes,-Gene_annotation)

genes_onlycategory_transpose <- t(genes_onlycategory)
genes_long <- pivot_longer(genes_onlycategory, c(PHASTERYO002_1:PHASTERYO384_1),names_to = "Prophage label", values_to = "")


### Figures and stats in order ###

# Figure 1
stats <- read.csv("Intact_v_incomplete.csv", header=TRUE,fileEncoding="UTF-8-BOM")
stats$Completeness <- factor(stats$Completeness, levels = c("Incomplete","Intact"))

stats$GC_range <- factor(stats$GC_range, levels = c("52_to_57","57_to_59","59_to_61","61_to_63","63_to_65","65_to_67","67_to_69","More_than_69"))
ggplot(data=stats, aes(x = GC_range, y= frequency(GC_range))) + 
  geom_bar(aes(fill=Completeness),stat='identity')+
  coord_cartesian(ylim = c(0, 270))+ 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("GC content of prophages (%)")+
  ylab("Number of prophages")+
  labs(fill="Completeness")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"))+
  scale_fill_manual(values=c("#88CCEE","#6699CC","#332288"))+ 
  theme(legend.position = c(0.17, 0.87))+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))

stats$Length_range <- factor(stats$Length_range, levels = c("3-5","5-15","15-25","25-35","35-45","45-55","More-55"))

ggplot(data=stats, aes(x = Length_range, y= frequency(Length_range))) +
  geom_bar(aes(fill=Completeness),stat='identity')+
  coord_cartesian(ylim = c(0, 270))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Length of prophages (kb)")+
  labs(fill="Completeness")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"))+
  scale_fill_manual(values=c("#88CCEE","#6699CC","#332288"))+
  theme(legend.position="none")+
  theme(axis.title.y = element_blank())

A + B + plot_layout(guides="collect")
b
ggarrange(A,B,common.legend = TRUE,legend="none")

incomplete_length <- subset(stats, Completeness=="Incomplete")

incomplete_length$Length <- as.numeric(incomplete_length$Length)

intact <- subset(stats, Completeness=="Intact")

intact_length$Length <- as.numeric(intact_length$Length)
mean(intact_length$Length)

mean(incomplete_length$Length)
hist(incomplete_length$Length)
hist(intact_length$Length)
var.test(incomplete_length$Length,intact_length$Length)
wilcox.test(incomplete_length$Length,intact_length$Length,paired=FALSE)
ggplot(aes(x=))

mod <- lm()


bc <- boxcox(Prophage_length_.bp. ~ Prophage_ID_Mashdist_0.1+Completeness,data=Que_incl_int_comp)
(lambda <- bc$x[which.max(bc$y)])

Length_aov <- aov(((Prophage_length_.bp.^lambda-1)/lambda) ~ Prophage_ID_Mashdist_0.1 + Completeness, data = Que_incl_int_comp)
summary(Length_aov)

TukeyHSD(Length_aov, which = "Completeness")

plot(res.aov2, 1)
plot(res.aov2, 2)



### Mash matrix

library(pheatmap)

matrix <- read.csv("Intact_prophage_validated_mashmatrix.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)
m <- pheatmap(mat, show_rownames = TRUE, show_colnames = FALSE)
m
k <- data.frame(rownames(mat[m$tree_row[["order"]],]))
write.csv(k,'Intact_mash_Heatmap_order.csv')


matrix <- read.csv("Intact_prophage_validated_mashmatrix.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix[, -1];
row.names(matrix2) <- matrix$X
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix <- melt(matrix3)
head(melted_matrix)

ggplot(data = melted_matrix, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev) +
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank())+
  scale_fill_gradient(low="grey20",high="grey100")+
  labs(fill="Genetic dissimilarity \n(Mash distance)")+
  theme(legend.title = element_text(colour="black", size=13, face="bold"),panel.border = element_rect(colour = "black", fill=NA))



## Mash tree
library(patchwork)
### Figure 3

tree <- read.tree("Intact_prophages_validated_mashtree.dnd")
rooted.tree <- root(tree, which(tree$tip.label == "Burkholderia_KS10"))
mashtree <- ggtree(rooted.tree,layout = "rectangular")+ theme(plot.margin=margin(0,0,0,10))+ ylab("Prophage Mash distance tree")+
  theme(axis.title=element_text(size=12,face="bold"))
mashtree
mashtree <- ggtree(rooted.tree,layout = "rectangular")+ theme(plot.margin=margin(0,0,0,10))+ 
  geom_balance(node=478, fill="#0072B2", color=NA,alpha=.2,extendto=1.1) +
  geom_balance(node=361, fill="#009E73", color=NA,alpha=.2,extendto=1.1)+
  geom_balance(node=422, fill="#D55E00", color=NA,alpha=.2,extendto=1.1)+
  geom_balance(node=467, fill="#009E73", color=NA,alpha=.2,extendto=1.1)+
  theme(panel.border = element_rect(colour = "black", fill=NA))

mashtree

tip_labs <- get_taxa_name(mashtree)
write.csv(tip_labs,"Mashtree_tiplabel.csv")

mashtree_rev <- mashtree + coord_flip()
mashtree_rev
new_tree <- mashtree_rev + scale_x_reverse()
new_tree

# Phage cornerstone genes
matrix <- read.csv("Phage_cornerstone_genes.csv",fileEncoding="UTF-8-BOM")
matrix_subset <- subset(matrix,select=c("Category","Phage_cell_lysis","Phage_DNA_replication_and_packaging","Phage_structural_protein"))
matrix2 <- matrix_subset[, -1];
row.names(matrix2) <- matrix_subset$Category
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix <- melt(matrix3)
head(melted_matrix)

cornerstone <- ggplot(data = melted_matrix, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev)+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="top",axis.ticks.y=element_blank(),legend.spacing.y=unit(0.5,"cm"))+  theme(plot.margin=margin(0,0,0,20))+
  labs(fill="Copy number")+
  scale_fill_gradient(low="grey90",high="grey20")+
  xlab("Prophage cornerstone gene categories")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.key=element_rect(colour="black",size=1),panel.border = element_rect(colour = "black", fill=NA))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(breaks=c("Phage_cell_lysis","Phage_DNA_replication_and_packaging","Phage_structural_protein"),
                   labels=c("Cell lysis", "DNA replication + packaging","Structural genes"))+
  theme(axis.title=element_text(size=12,face="bold"))
mashtree + cornerstone + plot_layout(widths = c(0.5, 0.6))


# Gene presence/absence and GC length plot

matrix <- read.csv("Intact_prophage_genepresenceabsence.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix[, -1];
row.names(matrix2) <- matrix$Prophage
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix <- melt(matrix3)
head(melted_matrix)

tmp = group_by(melted_matrix,Var1, Var2) %>% 
  summarise(., PA=as.factor(ifelse(sum(value)>0,1,0)))

gene_matrix<- ggplot(data = tmp, aes(x=Var2, y=Var1, fill=PA)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev)+ 
  scale_x_discrete(limits = rev)+ 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank())+
  scale_fill_manual(labels = c("Gene present", "Gene absent"),values=c("grey98","navyblue"))+
  theme(legend.key=element_rect(colour="black", size=1),panel.border = element_rect(colour = "black", fill=NA))+
  theme(plot.margin=margin(20,5,0,50))+
  theme(legend.title=element_text(size=10,face="bold"))+
  labs(fill="Gene presence")+
  coord_flip()


lengthGC <- read.csv("Intact_prophage_lengthGC.csv",fileEncoding="UTF-8-BOM")
lengthGC$Prophage <- factor(lengthGC$Prophage, as.character(lengthGC$Prophage))

length_data <- subset(lengthGC, select=c("Prophage","Length"))
GC_data <- subset(lengthGC, select=c("Prophage","GC"))

length_heatmap <- ggplot(length_data, aes(x = 1,y = Prophage, fill=Length)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev)+ 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),legend.position="right")+
  theme(plot.margin=margin(0,0,0,5))+
  scale_fill_gradient(low="blue4", high="steelblue1",limits=c(5000,65000))+
  theme(legend.title=element_text(size=10,face="bold"),panel.border = element_rect(colour = "black", fill=NA))+
  labs(fill="Length (bp)")+
  coord_flip()

GC_heatmap <- ggplot(GC_data, aes(x = 1,y = Prophage, fill=GC)) + 
  geom_tile() + 
  scale_y_discrete(limits = rev)+ 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank())+ 
  theme(plot.margin=margin(0,0,0,5),legend.position="right")+
  scale_fill_gradient(low="yellow", high="red",limits=c(55,70))+
  theme(legend.title=element_text(size=10,face="bold"),panel.border = element_rect(colour = "black", fill=NA))+
  labs(fill="GC content (%)")+
  coord_flip()

new_tree + gene_matrix + length_heatmap + GC_heatmap + plot_layout(heights = c(0.5,1,0.15,0.15,0.15),guides = 'collect', ncol=1)


genecontent <- read.csv("Cluster_gene_content.csv",fileEncoding="UTF-8-BOM")
gene_content_long

## Intact vs incomplete GC length - FigureS2

# GC content and length difference between intact, questionable, and incomplete prophages

Que_incl_int_comp <- read.csv("Intact_vs_incomplete_GClength_onlyphaster.csv",fileEncoding="UTF-8-BOM")
Que_incl_int_comp$Completeness <- factor(Que_incl_int_comp$Completeness, levels = c("Intact","Incomplete"))
Que_incl_int_comp$Taxonomic_identity <- factor(Que_incl_int_comp$Taxonomic_identity, levels = c("Ralstonia phage RSM3","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified C","Unclassified D","Unclassified F"))


# GC content

b <- ggplot(data=Que_incl_int_comp, aes(x = Taxonomic_identity, y= Length,fill=Completeness))  + 
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Prophage length (bp)") +
  xlab("Prophage group") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_discrete(name = "Completeness", labels = c("Intact", "Incomplete"),drop=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) + 
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) + 
  labs(fill="Completeness")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(breaks=c("Ralstonia phage RSM3","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified C","Unclassified D","Unclassified F"),
                   labels=c(expression(paste(phi," RSM3")),expression(paste(phi," RSS1")),expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),"Unclassified D","Unclassified C","Unclassified F"))



group_by(Que_incl_int_comp, Taxonomic_identity, Completeness) %>%
  summarise(
    count = n(),
    mean = mean(Length, na.rm = TRUE),
    sd = sd(Length, na.rm = TRUE)
  )

bc <- boxcox(Length ~ Taxonomic_identity+Completeness,data=Que_incl_int_comp)
(lambda <- bc$x[which.max(bc$y)])

Length_aov <- aov(((Length^lambda-1)/lambda) ~ Taxonomic_identity + Completeness, data = Que_incl_int_comp)
summary(Length_aov)

TukeyHSD(Length_aov, which = "Completeness")

plot(Length_aov, 1)
plot(Length_aov, 2)


# GC

a <- ggplot(data=Que_incl_int_comp, aes(x = Taxonomic_identity, y= GC,fill=Completeness))  + 
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Prophage GC content (%)") +
  xlab("Prophage group") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_discrete(name = "Completeness", labels = c("Intact", "Incomplete"),drop=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) + 
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
  labs(fill="Completeness")


group_by(Que_incl_int_comp, Taxonomic_identity, Completeness) %>%
  summarise(
    count = n(),
    mean = mean(GC, na.rm = TRUE),
    sd = sd(GC, na.rm = TRUE)
  )


bc <- boxcox(GC ~ Taxonomic_identity+Completeness,data=Que_incl_int_comp)
(lambda <- bc$x[which.max(bc$y)])

res.aov2 <- aov(((GC^lambda-1)/lambda) ~ Taxonomic_identity + Completeness, data = Que_incl_int_comp)
summary(res.aov2)

TukeyHSD(res.aov2, which = "Completeness")

plot(res.aov2, 1)
plot(res.aov2, 2)

ggarrange(b,a,common.legend = TRUE,legend="right")

a /b + plot_layout(guides="collect")


#Figure S3

# Phylotype prophage abundance figures
Phylotype_prophage_number <- read.csv("Phylotype_prophage_number.csv",fileEncoding="UTF-8-BOM")
Phylotype_prophage_number_long <- gather(Phylotype_prophage_number, Proph_completeness, Proph_freq, Intact:Incomplete, factor_key=TRUE)

ggplot(data=Phylotype_prophage_number_long, aes(x = Phylotype, y= Proph_freq,fill=Proph_completeness))  + 
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Number of prophages per genome") +
  xlab("Phylotype") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_discrete(name = "Completeness", labels = c("Intact","Incomplete")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) + 
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) + 
  labs(fill="Completeness")                                                                                                                                                                                                                                                                                                                                                                                


### CHECK WITH VILLE - DIFFERENT DISTRIBUTIONS SO KW MIGHT NOT BE BEST
intact_prophages <- Phylotype_prophage_number_long[which(Phylotype_prophage_number_long$Proph_completeness=='Intact'),]

hist(intact_prophages[which(intact_prophages$Phylotype=='I'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='IIA'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='IIB'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='III'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='IV'),]$Proph_freq)

kruskal.test(Proph_freq ~ Phylotype, data = intact_prophages)
dunnTest(Proph_freq ~ Phylotype, data = intact_prophages)

incomplete_prophages <- Phylotype_prophage_number_long[which(Phylotype_prophage_number_long$Proph_completeness=='Incomplete'),]
kruskal.test(Proph_freq ~ Phylotype, data = incomplete_prophages)
dunnTest(Proph_freq ~ Phylotype, data = incomplete_prophages)

questionable_prophages <- Phylotype_prophage_number_long[which(Phylotype_prophage_number_long$Proph_completeness=='Questionable..n...73.'),]
kruskal.test(Proph_freq ~ Phylotype, data = questionable_prophages)
dunnTest(Proph_freq ~ Phylotype, data = questionable_prophages)


# Figure S4

tree <- read.tree("Intact_prophages_validated_ref_mashtree.dnd")
rooted.tree <- root(tree, which(tree$tip.label == "Burkholderia_KS10"))

tip_colours <- read.csv("Figure_S3_tiplabels.csv",fileEncoding="UTF-8-BOM")
colours <- data.frame(tip_colours$Colour)

mashtree <- ggtree(rooted.tree,layout = "rectangular")+ ylab("Prophage Mash distance tree")+theme(axis.title=element_text(size=12,face="bold"))+ geom_tiplab(aes(
  subset=(grepl('Ralstonia',label,fixed=TRUE)==TRUE)),colour=colours$tip_colours.Colour, size = 4)+ ggplot2::xlim(0, 2)#geom_tiplab(colour=colours$tip_colours.Colour, size = 5) + ggplot2::xlim(0, 2)
mashtree

tip_labs <- get_taxa_name(mashtree)
write.csv(tip_labs,"Mashtree_refs__tiplabel.csv")



# World map plot

#Figure 4

install.packages("Polychrome")
library(Polychrome)
help(Polychrome)
data(alphabet)
names(alphabet) <- NULL


c19 <- c("dodgerblue2", "#E31A1C", "green4", "maroon", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "steelblue4", "orchid1", "darkturquoise", "khaki2", "brown","yellow4","#6A3D9A","green1")
world <- map_data('world')

p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="antiquewhite", color="black") +
  coord_quickmap()+theme(panel.background = element_rect(fill = "aliceblue"),panel.border = element_rect(colour = "black", fill=NA, size=1))+ labs(x = "Longitude") + labs(y = "Latitude") +
  theme(axis.text.y = element_text(size = 10, color = "black")) + 
  theme(axis.text.x = element_text(size = 10, color = "black")) + 
  theme(axis.title = element_text( face="bold", size=14))+theme(legend.title = element_text(color="black",size=12,face="bold"))+theme(legend.text = element_text(size=12))+labs(fill="Prophage")

p
pies <- read.csv("scatterpie_continent_moved2.csv")

pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

library("colorBlindness")
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

colours <- paletteMartin
newcolours <- colours[-c(1,12)]
newcolours
names(newcolours) <- NULL   

newcolours

p + geom_scatterpie(aes(x=?..long,y=lat,group=region,r=15),data=pies,cols=c("RS551","PE226","RSM3","RSS1","RSS30","RSA1","RsoM1USA","RSY1","RS138","Dina","Unclassified.C","Unclassified.A","Unclassified.F"),alpha=.8)+ 
  scale_fill_manual(values = newcolours,labels=c(expression(paste(phi," RS551")),expression(paste(phi," PE226")),expression(paste(phi," RSM3")),expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),expression(paste(phi," RS138")),expression(paste(phi," Dina")),"Unclassified C","Unclassified A","Unclassified F"))+ 
  annotate(geom="text", x=-10, y=-60, label="(87)",color="black",fontface="bold")+ 
  annotate(geom="text", x=180, y=5, label="(63)",color="black",fontface="bold")+ 
  annotate(geom="text", x=-30, y=25, label="(89)",color="black",fontface="bold")+ 
  annotate(geom="text", x=80, y=-60, label="(15)",color="black",fontface="bold")+ 
  annotate(geom="text", x=-165, y=-20, label="(28)",color="black",fontface="bold")+ 
  annotate(geom="text", x=-120, y=-60, label="(26)",color="black",fontface="bold")+
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  theme(legend.justification = "top") +
  theme(legend.text.align = 0)
#scale_colour_manual(labels=c(expression(paste(phi," RS551")),expression(paste(phi," PE226")),expression(paste(phi," RSM3")),expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSS-TH1")),expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),expression(paste(phi," RS138")),expression(paste(phi," Dina")),"Unclassified C","Unclassified A"))
#theme(legend.position="right")+
#guides(fill=guide_legend(ncol=10))

## Same again but as a bar chart

library(tidyr)
Prophage_continent_proportion <- read.csv("Continent_distributions_prophages.csv",fileEncoding="UTF-8-BOM")
Gathered <- pivot_longer(Prophage_continent_proportion, cols=c(Africa:Oceania),names_to = "Continents",values_to="Proportion")

Gathered$Prophage <- factor(Gathered$Prophage,levels = c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))


ggplot(data=Gathered, aes(x=Prophage,y=Proportion)) + 
  scale_y_continuous(labels = scales::percent, limits=c(0,1))+
  geom_bar(stat="identity", color="black",fill="lightblue")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  ylab("Percent of total prophages in continent (%)") +
  xlab("Prophage")+
  theme(axis.title=element_text(size=13,face="bold"),plot.title = element_text(size=15,face="bold"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(legend.position="none")+
  facet_grid(rows=vars(Continents))+
  scale_x_discrete(breaks=c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"),
                   labels=c(expression(paste(phi," RS551")),expression(paste(phi," PE226")),expression(paste(phi," RSM3")),expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSS-TH1")), "Unclassified B",expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),"Unclassified D",expression(paste(phi," RS138")),"Unclassified G",expression(paste(phi," Dina")),"Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))




# Figure 5A - re-done

install.packages("rlang")
require(ggtree)
install.packages("phytools")
require(phytools)
library(patchwork)
tree <- read.tree("core_tree.treefile")
rooted.tree <- midpoint.root(tree)
tip_labs <- rooted.tree$tip.label

to_drop <- c("UY031","K60","GMI1000","CMR15","PSI07")
rooted.tree_reduced <- drop.tip(rooted.tree,to_drop)

ggtree(rooted.tree_reduced,layout = "rectangular",ladderize=TRUE)+ geom_text(aes(label=node), hjust=-.3)

tree_fig <- ggtree(rooted.tree_reduced,layout = "rectangular",ladderize=TRUE)+ 
  theme(plot.margin=margin(0,0,0,0)) + ylab("RSSC phylogeny")+
  theme(axis.title=element_text(size=12,face="bold"))

tree_fig_collapse <- tree_fig %>%
  collapse(node=201,"max")%>%
  collapse(node=265,"max")

tree_fig_collapse_colour <- tree_fig_collapse %>%
  + geom_hilight(node=268, fill="steelblue", alpha=.3, extend=1) +
  geom_hilight(node=363, fill="darkgreen", alpha=.3, extend=1)+
  geom_hilight(node=203, fill="yellow", alpha=.3, extend=1)+
  geom_hilight(node=196, fill="red", alpha=.3, extend=1)+
  geom_hilight(node=258, fill="orange", alpha=.3, extend=1)
tree_fig_collapse_colour

matrix <- read.csv("R_pickettii_tree_presence_absence_matrix_validated_intact_prophages.csv")
matrix2 <- matrix[, -1];
row.names(matrix2) <- matrix$?..Isolate
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix <- melt(matrix3)
head(melted_matrix)
melted_matrix$value <- as.factor(melted_matrix$value)


All_intact_phylogeny <- ggplot(data = melted_matrix, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev)+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),legend.spacing.y=unit(0.5,"cm"))+
  scale_fill_manual(values=c("grey98","grey40","firebrick"),labels=c("Absent","Present (one copy)","Present (two copies)")) + 
  theme(plot.margin=margin(0,0,0,30))+
  labs(fill="Prophage")+
  #xlab("Prophage type")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.key=element_rect(colour="black",size=1),panel.border = element_rect(colour = "black", fill=NA))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(breaks=c("Ralstonia.phage.Rs551","Ralstonia.phage.PE226","Ralstonia.phage.RSM3","Ralstonia.phage.RSS1","Ralstonia.phage.RSS30","Ralstonia.phage.RSS.TH1","Unclassified.B","Ralstonia.phage.phiRSA1","Ralstonia.phage.RsoM1USA","Ralstonia.phage.RSY1","Unclassified.D","Ralstonia.phage.RS138","Unclassified.G","Ralstonia.phage.Dina","Unclassified.C","Unclassified.H","Unclassified.A","Unclassified.F","Unclassified.I","Unclassified.J"),
                   labels=c(expression(paste(phi," RS551")),expression(paste(phi," PE226")),expression(paste(phi," RSM3")),expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSS-TH1")), "Unclassified B",expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),"Unclassified D",expression(paste(phi," RS138")),"Unclassified G",expression(paste(phi," Dina")),"Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))+
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.justification = "top")


Fig_5a <- tree_fig + All_intact_phylogeny + plot_layout(widths = c(0.4, 1))+plot_layout(guides="collect")


#Figure 5 supplementary incomplete

tree <- read.tree("core_tree.treefile")
rooted.tree <- midpoint.root(tree)
tip_labs <- rooted.tree$tip.label

to_drop <- c("UY031","K60","GMI1000","CMR15","PSI07")
rooted.tree_reduced <- drop.tip(rooted.tree,to_drop)

tree_fig <- ggtree(rooted.tree_reduced,layout = "rectangular",ladderize=TRUE)+ 
  theme(plot.margin=margin(0,0,0,0)) + ylab("RSSC phylogeny")+
  theme(axis.title=element_text(size=12,face="bold"))


matrix <- read.csv("R_pickettii_tree_presence_absence_matrix_incomplete_prophages.csv")
matrix2 <- matrix[, -1];
row.names(matrix2) <- matrix$?..Isolate
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix <- melt(matrix3)
head(melted_matrix)
melted_matrix$value <- as.factor(melted_matrix$value)


All_intact_phylogeny <- ggplot(data = melted_matrix, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev)+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),legend.spacing.y=unit(0.5,"cm"))+
  scale_fill_manual(values=c("grey98","grey40","salmon","firebrick"),labels=c("Absent","One copy","Two copies","Three copies")) + 
  theme(plot.margin=margin(0,0,0,30))+
  labs(fill="Prophage \npresence",x="Prophage type")+
  #xlab("Prophage type")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.key=element_rect(colour="black",size=1),panel.border = element_rect(colour = "black", fill=NA))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(breaks=c("Ralstonia.phage.Rs551","Ralstonia.phage.PE226","Ralstonia.phage.RSM3","Ralstonia.phage.RSS1","Ralstonia.phage.RSS30","Ralstonia.phage.RSS.TH1","Unclassified.B","Ralstonia.phage.phiRSA1","Ralstonia.phage.RsoM1USA","Ralstonia.phage.RSY1","Unclassified.D","Ralstonia.phage.RS138","Unclassified.G","Ralstonia.phage.Dina","Unclassified.C","Unclassified.H","Unclassified.A","Unclassified.F","Unclassified.I","Unclassified.J"),
                   labels=c(expression(paste(phi," RS551")),expression(paste(phi," PE226")),expression(paste(phi," RSM3")),expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSS-TH1")), "Unclassified B",expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),"Unclassified D",expression(paste(phi," RS138")),"Unclassified G",expression(paste(phi," Dina")),"Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))+
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.justification = "top")


tree_fig + All_intact_phylogeny + plot_layout(widths = c(0.4, 1),guides="collect")


## Figure_5 supplementary figure - Figure S7

Prophage_phylotype_proportion <- read.csv("Phylotype_distributions_prophages_dividephylotypeabundance.csv",fileEncoding="UTF-8-BOM")
Gathered <- pivot_longer(Prophage_phylotype_proportion, cols=c(I:IV),names_to = "Phylotypes",values_to="Proportion")

Gathered$Prophage <- factor(Gathered$Prophage,levels = c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

cols <- c("#66CC33","#CCFFFF","#FF9900","#CC0000","#003399","#663300","#FFFF00","#660066","#006600","#999999")

ggplot(data=Gathered, aes(x=Prophage,y=Proportion,fill=Cluster)) + 
  scale_y_continuous(labels = scales::percent, limits=c(0,1))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=cols)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  ylab("Percent of isolates containing prophage (%)") +
  xlab("Prophage")+
  theme(axis.title=element_text(size=13,face="bold"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.justification = "top")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"))+
  facet_grid(rows=vars(Phylotypes))+
  scale_x_discrete(breaks=c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"),
                   labels=c(expression(paste(phi," RS551")),expression(paste(phi," RSM3")),"Unclassified B",expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSS-TH1")),expression(paste(phi," PE226")),expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),"Unclassified D",expression(paste(phi," RS138")),"Unclassified G",expression(paste(phi," Dina")),"Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

#Figure_S8

BC_and_mash <- read.csv("Phylotype_BC_mash_regression.csv")

ggplot(data=BC_and_mash, aes(x = Phylotype, y= BC))  + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  # geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+ylab("Average prophage diversity") + xlab("Phylotype") + theme_bw()+ geom_smooth(method = "lm")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title=element_text(size=12,face="bold"))+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill="Completeness")+ theme(legend.title = element_text(color = "black", size = 10,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                            legend.text = element_text(color = "black",size=10,face="bold"))                                                                                                                                                                                                                                                                                                                                                                                

ggplot(data=BC_and_mash, aes(x = Phylotype, y= Mash))  + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  # geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+ylab("Average host diversity") + xlab("Phylotype") + theme_bw()+ geom_smooth(method = "lm")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title=element_text(size=12,face="bold"))+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill="Completeness")  + theme(legend.title = element_text(color = "black", size = 10,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                          legend.text = element_text(color = "black",size=10,face="bold"))                                                                                                                                                                                                                                                                                                                                                                              

ggarrange(a,b)

kruskal.test(BC_and_mash$Mash~BC_and_mash$Phylotype)
dunnTest(BC_and_mash$Mash~BC_and_mash$Phylotype)
kruskalmc(BC_and_mash$Mash~BC_and_mash$Phylotype,probs=0.05)

kruskal.test(BC_and_mash$BC~BC_and_mash$Phylotype)
dunnTest(BC_and_mash$BC~BC_and_mash$Phylotype)
kruskalmc(BC_and_mash$BC~BC_and_mash$Phylotype,probs=0.05)

# Figure 5b

# Make Bray-Curtis distance matrices
matrix <- read.csv("PhyloI_presenceabsence.csv",fileEncoding="UTF-8-BOM")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)
x <- vegdist(mat, method="bray")
x[is.na(x)] <- 0
PhyloI_BC <- as.matrix(x)
write.csv(PhyloI_BC, "PhyloI_BCmatrix.csv")


matrix <- read.csv("PhyloIIA_presenceabsence.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)
x <- vegdist(mat, method="bray")
x[is.na(x)] <- 0
PhyloIIA_BC <- as.matrix(x)
write.csv(PhyloIIA_BC, "PhyloIIA_BCmatrix.csv")


matrix <- read.csv("PhyloIIB_presenceabsence.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)
x <- vegdist(mat, method="bray")
x[is.na(x)] <- 0
PhyloIIB_BC <- as.matrix(x)
write.csv(PhyloIIB_BC, "PhyloIIB_BCmatrix.csv")

# Take averages of each row and add to spreadsheet with host mash distances
reg <- read.csv("Phylotype_BC_mash_regression.csv")
lm2 <- lme(BC ~ log(Mash),data=reg, random=~1|Phylotype)

coef(lm2)
par(mfrow=c(2,2))
hist(lm2$residuals)
qqnorm(lm2$residuals)
qqline(lm2$residuals)
plot(lm2$fitted,lm2$residuals)
summary(lm2)
anova(lm2)

library("piecewiseSEM")
install.packages("glue")
detach("glue", unload=TRUE,character.only=TRUE)
rsquared(lm2)

my_y_title <- expression(paste("Average ", italic("R. solanacearum"), " dissimilarity (log",[10]")"))

library("ggtext")
ggplot(reg,aes(x=log(Mash),y=BC))+ 
  geom_point(aes(fill=Phylotype),colour="black",pch=21, size=3)+
  xlim(c(-6.5,-3.5))+ylim(c(0.125,1))+
  labs(y="Average prophage dissimilarity",x="Average <i>R. solanacearum</i> dissimilarity (log<sub>10</sub>)")+
  theme_bw()+
  theme(axis.title.x = element_markdown())+
  geom_segment(aes(x=-6.06,xend=-3.888,y=0.125,yend=1),size=1)+
  theme(legend.position = c(0.15,0.83))+
  theme(axis.title=element_text(size=13,face="bold"))+
  theme(legend.title = element_text(colour="black", size=16,face="bold"),legend.text = element_text(size=10))
# + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




## Congruence figure 5C
library(vegan)
# Figure 5c
## PACo analysis of phylotype I clade

# All isolate phylogeny
tree <- read.tree("core_tree.treefile")
rooted.tree <- midpoint.root(tree)

# Extract phylotype I clade
rooted.tree.phyloI <- extract.clade(rooted.tree,209)
rooted.tree.phyloI.drop <- drop.tip(rooted.tree.phyloI,"GMI1000")

matrix <- read.csv("PhyloI_presenceabsence.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")


PhyloI_BC <- as.matrix(x)
PhyloI_BC[is.na(PhyloI_BC)] <- 0
PhyloI_UPGMA <- upgma(PhyloI_BC)

# Start running PACo
host.D <- cophenetic(rooted.tree.phyloI.drop)
BC.D <- cophenetic(PhyloI_UPGMA)

HP <- read.csv("PACo_binarymatrix_phyloI.csv")
row.names(HP) <- HP$?..
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

## PACo analysis of phylotype IIA clade
rooted.tree.phyloIIA <- extract.clade(rooted.tree,372)
rooted.tree.phyloIIA.drop <- drop.tip(rooted.tree.phyloIIA,"K60")
circ1 <- ggtree(rooted.tree.phyloIIA,layout = "rectangular",ladderize=TRUE)+ theme(plot.margin=margin(0,0,0,0))+ geom_text(aes(label=node))
circ1

matrix <- read.csv("PhyloIIA_presenceabsence.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloIIA_BC <- as.matrix(x)
PhyloIIA_BC[is.na(PhyloIIA_BC)] <- 0
PhyloIIA_UPGMA <- upgma(PhyloIIA_BC)

# Start running PACo
host.D <- cophenetic(rooted.tree.phyloIIA.drop)
BC.D <- cophenetic(PhyloIIA_UPGMA)

HP <- read.csv("PACo_binarymatrix_phyloIIA.csv")
row.names(HP) <- HP$?..
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

## PACo analysis of phylotype IIB clade

rooted.tree.phyloIIB <- extract.clade(rooted.tree,276)
rooted.tree.phyloIIB.drop <- drop.tip(rooted.tree.phyloIIB,"UY031")

matrix <- read.csv("PhyloIIB_presenceabsence.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloIIB_BC <- as.matrix(x)
PhyloIIB_BC[is.na(PhyloIIB_BC)] <- 0
PhyloIIB_UPGMA <- upgma(PhyloIIB_BC)

# Start running PACo
host.D <- cophenetic(rooted.tree.phyloIIB.drop)
BC.D <- cophenetic(PhyloIIB_UPGMA)

HP <- read.csv("PACo_binarymatrix_phyloIIB.csv")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof




#PACo analysis of all isolates together

tree <- read.tree("core_tree.treefile")
rooted.tree <- midpoint.root(tree)
tip <- c("K60","GMI1000","UY031","CMR15","PSI07")
rooted.tree.drop <- drop.tip(rooted.tree,tip)

# Use presence/absence matrix with reference isolates removed
matrix <- read.csv("R_pickettii_tree_presence_absence_matrix_forcongruence.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
All_BC <- as.matrix(x)
All_BC[is.na(All_BC)] <- 0
AllBC_UPGMA <- upgma(All_BC)

# Start running PACo
host.D <- cophenetic(rooted.tree.drop)
BC.D <- cophenetic(AllBC_UPGMA)

HP <- read.csv("PACo_binarymatrix_All.csv")
row.names(HP) <- HP$?..
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
weight <- (res^-2)/50

# Make figure 5c tanglegram
rooted.tree.drop2 <- rotateNodes(rooted.tree.drop,"all")
rooted.tree.drop3 <- untangle(rooted.tree.drop2)

AllBC_UPGMA2 <- rotateNodes(AllBC_UPGMA,"all")
AllBC_UPGMA3 <- untangle(AllBC_UPGMA2)

Figure_5C <- cophyloplot(rooted.tree.drop3, AllBC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
                         lwd=1, col='steelblue', length.line=0, gap=0, space=70,rotate=TRUE)+coord_flip()

Figure_5C

Fig_5a
Fig_5a + Fig_5b + plot_layout(ncol=2,widths=c(1,0.5))



# Figure 6

#VIGA results
VIGA <- read.csv("VIGA_proportion_table.csv",fileEncoding="UTF-8-BOM")
VIGA_long <- gather(VIGA, Annotation_hit, Annotation_freq, DNA.replication:Other, factor_key=TRUE)

VIGA_long$Annotation_hit <- factor(VIGA_long$Annotation_hit,levels=c("Other","Hypothetical.protein.Unknown.function","Virulence","Bacterial.stress.tolerance","Membrane.bound.protein","Protein.secretion","Cellular.metabolism","Post.translational.modification","Protein.degradation","DNA.methyltransferase","RNA.binding","DNA.binding","DNA.replication","Phage.cell.lysis","Phage.structural.protein"))
VIGA_long$Prophage <- factor(VIGA_long$Prophage, levels=c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

VIGA_long$Annotation_hit <- gsub('.', ' ', VIGA_long$Annotation_hit, fixed=TRUE)
VIGA_long$Annotation_hit <- gsub('Hypothetical protein ', '', VIGA_long$Annotation_hit, fixed=TRUE)

VIGA_long$Annotation_hit <- factor(VIGA_long$Annotation_hit,levels=c("Other","Unknown function","Virulence","Bacterial stress tolerance","Membrane bound protein","Protein secretion","Cellular metabolism","Post translational modification","Protein degradation","DNA methyltransferase","RNA binding","DNA binding","DNA replication","Phage cell lysis","Phage structural protein"))


ggplot(data=VIGA_long, aes(x=Prophage,y=Annotation_freq,fill=Annotation_hit)) +
  scale_y_continuous(labels = scales::percent, position="right")+
  geom_bar(aes(fill=Annotation_hit),position='fill',stat='identity', na.rm = FALSE,colour="black")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = mycolors)+
  ylab("") +
  xlab("Prophage")+ 
  labs(fill = "Gene function")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"), axis.ticks.x = element_blank())+
  theme(legend.title = element_text(colour="black", size=16,face="bold"), legend.position="left",legend.text = element_text(size=10))+
  guides(fill = guide_legend(label.position = "left", title.position="top", title.hjust=0.99))+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
  theme(plot.margin=margin(20,0,0,0))


# CAZy results
matrix_cazy <- read.csv("CAZy_heatmap.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix_cazy[, -1];
row.names(matrix2) <- matrix_cazy$Genes
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix_cazy <- melt(matrix3)
head(melted_matrix)

CAZy <- ggplot(data = melted_matrix_cazy, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(colour="black")+ 
  labs(fill="")+
  scale_y_discrete(limits = rev,breaks=c("CBM50","CBM50|GH23","GH1|GH95","GH100","GH104","GH13_30","GH23","GH24","GH58","GH73","GT1","GT2","GT4"),
                   labels=c("LysM domain (CBM50)","Lysozyme & LysM domain (CBM50|GH23)","Various (GH1|GH95)","Alkaline and neutral invertase (GH100)","Peptidoglycan lytic transglycosylase (GH104)","??-glucosidase (GH13_30)","Lysozyme (GH23)","Lysozyme (GH24)","Endo-N-acetylneuraminidase (GH58)","Lysozyme (GH73)","Various (GT1)","Various (GT2)","Glycogen synthase (GT4)"))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_gradient(low="grey100",high="darkgreen",limits=c(0, 1.5), breaks=seq(0,1.5,by=0.5))+
  theme(legend.justification = "top")+
  theme(axis.text.y=element_blank())+
  coord_fixed()


# RalstoT3E results
matrix_ralsto <- read.csv("Ralsto_T3E_heatmap.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix_ralsto[, -1];
row.names(matrix2) <- matrix_ralsto$Genes
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix_ralsto <- melt(matrix3)
head(melted_matrix)

RalstoT3E <- ggplot(data = melted_matrix_ralsto, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(colour="black")+ 
  scale_y_discrete(limits = rev)+ 
  theme(axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_gradient(low="grey100",high="firebrick",limits=c(0, 2), breaks=seq(0,2,by=0.5))+
  labs(fill="")+
  theme(axis.text.y=element_blank())+
  theme(legend.justification = "top")+
  coord_fixed()

# PHI-base results
matrix_phi <- read.csv("PHI_base_heatmap.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix_phi[, -1];
row.names(matrix2) <- matrix_phi$Genes
matrix3 <- as.matrix(matrix2)
head(matrix3)
melted_matrix_phi <- melt(matrix3)
head(melted_matrix_phi)

Phi_base <- ggplot(data = melted_matrix_phi, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(colour="black")+ 
  theme(axis.title.y=element_blank(),legend.position="right",axis.title.x=element_blank(),axis.text.x=element_blank())+
  scale_fill_gradient(low="grey100",high="darkblue",limits=c(0, 1.25), breaks=seq(0,1.25,by=0.5))+
  theme(legend.title=element_blank(), axis.ticks.x=element_blank())+
  coord_fixed()+
  labs(fill="")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  theme(legend.justification = "top")+
  theme(axis.text.y=element_blank())+
  scale_y_discrete(limits = rev,breaks=c("AldB_(PSPTO_2673)","EadM","EntA_(VK055_1922)","hupB_(XAC1081)","lolA","lpxO","msh2","PD0928","rarA","T6SS1","Trr1","vgrG-1"),
                   labels=c("Aldehyde dehydrogenase (PHI:7788; AldB","DNA methyltransferase (PHI:8921; EadM)","2,3-dihydroxybenzoate-2,3-dehydrogenase (PHI:7485; EntA)","Histone-like DNA binding protein (PHI:8676; HupB)","Outer-membrane lipoprotein carrier protein (PHI:8855; LolA)","Acyl hydroxylase (PHI:7948; LpxO)","DNA mismatch repair protein (PHI:7207; Msh2)","Zot-like toxin (PHI:4984; PD0928)","Replication-associated recombination protein (PHI:6447; RarA)"," Type VI secretion system protein (PHI:4558; T6SS1)","Thioredoxin reductase (PHI:6470; Trr1)","Type VI secretion system protein (PHI:7092; VgrG)"))

Prophage_phylo <- read.csv("Prophage_phylotype_heatmap.csv",fileEncoding="UTF-8-BOM")
Prophage_phylo$Phage_hits <- factor(Prophage_phylo$Phage_hits, levels=c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

colour_phylo <- c("khaki","lightsalmon","darkorchid3","forestgreen","gray80")
options(repr.plot.height =1)

options(repr.plot.width=6, repr.plot.height=4)

windows.options(height=3)
options(repr.plot.height = 1.5)
ggplot(data=Prophage_phylo, aes(x=Phage_hits,y=frequency(Phylotype),fill=Phylotype)) + 
  scale_y_continuous(labels = scales::percent)+
  geom_bar(position='fill',stat='identity', na.rm = FALSE)+
  scale_fill_manual(values=colour_phylo)+
  theme_bw() +
  xlab("Prophage type")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.y=element_blank())+
  theme(axis.title=element_text(size=12,face="bold"))+
  scale_x_discrete(breaks=c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"),
                   labels=c(expression(paste(phi," RS551")),expression(paste(phi," RSM3")),"Unclassified B",expression(paste(phi," RSS1")),expression(paste(phi," RSS30")),expression(paste(phi," RSS-TH1")),expression(paste(phi," PE226")),expression(paste(phi," RSA1")),expression(paste(phi," RsoM1USA")),expression(paste(phi," RSY1")),"Unclassified D",expression(paste(phi," RS138")),"Unclassified G",expression(paste(phi," Dina")),"Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12))


CAZy + RalstoT3E + Phi_base +plot_layout(ncol=1)
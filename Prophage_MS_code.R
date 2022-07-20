#### MS R script ####

### Dependencies ###

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
library("gridExtra")
library("nlme")
library("biobase")
library("piecewiseSEM")
library("scatterpie")
library("maps")
library("paco")
library("phytools")
library("RColorBrewer")
library("FSA")
library("vegan")
library("piecewiseSEM")
library("phangorn")
library("MASS")
library("dplyr")
library("tidyr")
library("patchwork")
library("ggvenn")
library("pheatmap")
library("Polychrome")
library("tidyr")
library("ggtext")

### Code for analysis ###

# Figure S2

prophage_mash_tree <- read.tree("Intact_prophages_validated_mashtree.dnd")
rooted_prophage_mash_tree <- root(prophage_mash_tree, which(tree$tip.label == "Burkholderia_KS10"))

rooted_prophage_mash_tree_figure <- ggtree(rooted.tree,layout = "rectangular")+ theme(plot.margin=margin(0,0,0,10))+ ylab("Prophage Mash distance tree")+
  theme(axis.title=element_text(size=12,face="bold"))

tip_labs <- get_taxa_name(rooted_prophage_mash_tree_figure)
write.csv(tip_labs,"Mashtree_tiplabel.csv")

phage_cornerstone_gene_matrix <- read.csv("Phage_cornerstone_genes.csv",fileEncoding="UTF-8-BOM")

phage_cornerstone_gene_matrix_subset <- subset(phage_cornerstone_gene_matrix,select=c("Category","Phage_cell_lysis","Phage_DNA_replication_and_packaging","Phage_structural_protein"))

phage_cornerstone_gene_matrix_subset2 <- phage_cornerstone_gene_matrix_subset[, -1];

row.names(phage_cornerstone_gene_matrix_subset2) <- phage_cornerstone_gene_matrix_subset$Category

phage_cornerstone_gene_matrix_subset3 <- as.matrix(phage_cornerstone_gene_matrix_subset2)
phage_cornerstone_gene_matrix_subset_melted <- melt(phage_cornerstone_gene_matrix_subset3)

cornerstone_prophage_genes_fig <- ggplot(data = phage_cornerstone_gene_matrix_subset_melted, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev)+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="top",axis.ticks.y=element_blank(),legend.spacing.y=unit(0.5,"cm"))+  theme(plot.margin=margin(0,0,0,20))+
  labs(fill="Copy number")+
  scale_fill_gradient(low="grey90",high="grey20")+
  xlab("Prophage cornerstone gene categories")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.key=element_rect(colour="black",size=1),panel.border = element_rect(colour = "black", fill=NA))+
  scale_x_discrete(breaks=c("Phage_cell_lysis","Phage_DNA_replication_and_packaging","Phage_structural_protein"),
                   labels=c("Cell lysis", "DNA replication + packaging","Structural genes"))+
  theme(axis.title=element_text(size=12,face="bold"))

rooted_prophage_mash_tree_figure + cornerstone_prophage_genes_fig + plot_layout(widths = c(0.5, 0.6))

# Figure 1

intvincl_general <- read.csv("Intact_v_incomplete.csv", header=TRUE,fileEncoding="UTF-8-BOM")
intvincl_general$Completeness <- factor(intvincl_general$Completeness, levels = c("Incomplete","Intact"))

intvincl_general$GC_range <- factor(intvincl_general$GC_range, levels = c("52_to_57","57_to_59","59_to_61","61_to_63","63_to_65","65_to_67","67_to_69","More_than_69"))

Fig1A <- ggplot(data=intvincl_general, aes(x = GC_range, y= frequency(GC_range))) + 
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

intvincl_general$Length_range <- factor(intvincl_general$Length_range, levels = c("3-5","5-15","15-25","25-35","35-45","45-55","More-55"))

Fig1B <- ggplot(data=intvincl_general, aes(x = Length_range, y= frequency(Length_range))) +
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

Fig1A + Fig1B + plot_layout(guides="collect")


# Figure S3

# GC content and length difference between related intact and incomplete prophages

intvincl_related <- read.csv("Intact_vs_incomplete_GClength_onlyphaster.csv",fileEncoding="UTF-8-BOM")
intvincl_related$Completeness <- factor(intvincl_related$Completeness, levels = c("Intact","Incomplete"))
intvincl_related$Taxonomic_identity <- factor(intvincl_related$Taxonomic_identity, levels = c("Ralstonia phage RSM3","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified C","Unclassified D","Unclassified F"))

# GC content

FigS3A <- ggplot(data=intvincl_related, aes(x = Taxonomic_identity, y= GC,fill=Completeness))  + 
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

group_by(intvincl_related, Taxonomic_identity, Completeness) %>%
  summarise(
    count = n(),
    mean = mean(GC, na.rm = TRUE),
    sd = sd(GC, na.rm = TRUE)
  )

bc <- boxcox(GC ~ Taxonomic_identity+Completeness,data=intvincl_related)
(lambda <- bc$x[which.max(bc$y)])

res.aov2 <- aov(((GC^lambda-1)/lambda) ~ Taxonomic_identity + Completeness, data = intvincl_related)
summary(res.aov2)

TukeyHSD(res.aov2, which = "Completeness")

plot(res.aov2, 1)
plot(res.aov2, 2)

# Length

FigS3B <- ggplot(data=intvincl_related, aes(x = Taxonomic_identity, y= Length,fill=Completeness))  + 
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

group_by(intvincl_related, Taxonomic_identity, Completeness) %>%
  summarise(
    count = n(),
    mean = mean(Length, na.rm = TRUE),
    sd = sd(Length, na.rm = TRUE)
  )

bc <- boxcox(Length ~ Taxonomic_identity+Completeness,data=intvincl_related)
(lambda <- bc$x[which.max(bc$y)])

Length_aov <- aov(((Length^lambda-1)/lambda) ~ Taxonomic_identity + Completeness, data = intvincl_related)
summary(Length_aov)

TukeyHSD(Length_aov, which = "Completeness")

plot(Length_aov, 1)
plot(Length_aov, 2)

FigS3A/FigS3B + plot_layout(guides="collect")


# Figure 2

mash_matrix <- read.csv("Intact_prophage_validated_mashmatrix.csv")
rownames(mash_matrix) <- mash_matrix[, 1];
mash_matrix2 <- mash_matrix[, -1];
mash_matrix3 <- as.matrix(mash_matrix2)

mash_heatmap <- pheatmap(mash_matrix3, show_rownames = TRUE, show_colnames = FALSE)
mash_heatmap

heatmap_labels <- data.frame(rownames(mash_matrix3[mash_heatmap$tree_row[["order"]],]))
write.csv(heatmap_labels,'Intact_mash_Heatmap_order.csv')

mash_matrix_reordered <- read.csv("Intact_prophage_validated_mashmatrix_reorder.csv",fileEncoding="UTF-8-BOM")
mash_matrix_reordered2 <- mash_matrix_reordered[, -1];
row.names(mash_matrix_reordered2) <- mash_matrix_reordered$X
mash_matrix_reordered3 <- as.matrix(mash_matrix_reordered2)
mash_matrix_reordered3_melted <- melt(mash_matrix_reordered3)

head(mash_matrix_reordered3_melted)

ggplot(data = mash_matrix_reordered3_melted, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ 
  scale_y_discrete(limits = rev) +
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank())+
  scale_fill_gradient(low="grey20",high="grey100")+
  labs(fill="Genetic dissimilarity \n(Mash distance)")+
  theme(legend.title = element_text(colour="black", size=13, face="bold"),panel.border = element_rect(colour = "black", fill=NA))


### Figure 3

# Mash tree
prophage_mash_tree <- read.tree("Intact_prophages_validated_mashtree.dnd")
rooted_prophage_mash_tree <- root(prophage_mash_tree, which(tree$tip.label == "Burkholderia_KS10"))

rooted_prophage_mash_tree_fig_familylabel <- ggtree(rooted_prophage_mash_tree,layout = "rectangular")+ theme(plot.margin=margin(0,0,0,10))+ 
  geom_balance(node=478, fill="#0072B2", color=NA,alpha=.2,extendto=1.1) +
  geom_balance(node=361, fill="#009E73", color=NA,alpha=.2,extendto=1.1)+
  geom_balance(node=422, fill="#D55E00", color=NA,alpha=.2,extendto=1.1)+
  geom_balance(node=467, fill="#009E73", color=NA,alpha=.2,extendto=1.1)+
  theme(panel.border = element_rect(colour = "black", fill=NA))

rooted_prophage_mash_tree_fig_familylabel

tip_labs <- get_taxa_name(rooted_prophage_mash_tree_fig_familylabel)
write.csv(tip_labs,"Mashtree_tiplabel.csv")

rooted_prophage_mash_tree_fig_familylabel_rev <- rooted_prophage_mash_tree_fig_familylabel + coord_flip()

rooted_prophage_mash_tree_fig_familylabel_rev_inverted <- rooted_prophage_mash_tree_fig_familylabel_rev + scale_x_reverse()

# Gene presence/absence and GC/length plot

gene_presence_matrix <- read.csv("Intact_prophage_genepresenceabsence.csv",fileEncoding="UTF-8-BOM")
gene_presence_matrix2 <- gene_presence_matrix[, -1];
row.names(gene_presence_matrix2) <- gene_presence_matrix$Prophage
gene_presence_matrix3 <- as.matrix(gene_presence_matrix2)
gene_presence_matrix3_melted <- melt(gene_presence_matrix3)

tmp = group_by(gene_presence_matrix3_melted,Var1, Var2) %>% 
  summarise(., PA=as.factor(ifelse(sum(value)>0,1,0)))

gene_presence_matrix_fig <- ggplot(data = tmp, aes(x=Var2, y=Var1, fill=PA)) + 
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

rooted_prophage_mash_tree_fig_familylabel_rev_inverted + gene_presence_matrix_fig + length_heatmap + GC_heatmap + plot_layout(heights = c(0.5,1,0.15,0.15,0.15),guides = 'collect', ncol=1)


#Figure S4

phylotype_prophage_number <- read.csv("Phylotype_prophage_number.csv",fileEncoding="UTF-8-BOM")
phylotype_prophage_number_long <- gather(phylotype_prophage_number, Proph_completeness, Proph_freq, Intact:Incomplete, factor_key=TRUE)

ggplot(data=phylotype_prophage_number_long, aes(x = Phylotype, y= Proph_freq,fill=Proph_completeness))  + 
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

intact_prophages <- phylotype_prophage_number_long[which(phylotype_prophage_number_long$Proph_completeness=='Intact'),]

hist(intact_prophages[which(intact_prophages$Phylotype=='I'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='IIA'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='IIB'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='III'),]$Proph_freq)
hist(intact_prophages[which(intact_prophages$Phylotype=='IV'),]$Proph_freq)

kruskal.test(Proph_freq ~ Phylotype, data = intact_prophages)
dunnTest(Proph_freq ~ Phylotype, data = intact_prophages)

incomplete_prophages <- phylotype_prophage_number_long[which(phylotype_prophage_number_long$Proph_completeness=='Incomplete'),]

hist(incomplete_prophages[which(intact_prophages$Phylotype=='I'),]$Proph_freq)
hist(incomplete_prophages[which(intact_prophages$Phylotype=='IIA'),]$Proph_freq)
hist(incomplete_prophages[which(intact_prophages$Phylotype=='IIB'),]$Proph_freq)
hist(incomplete_prophages[which(intact_prophages$Phylotype=='III'),]$Proph_freq)
hist(incomplete_prophages[which(intact_prophages$Phylotype=='IV'),]$Proph_freq)

kruskal.test(Proph_freq ~ Phylotype, data = incomplete_prophages)
dunnTest(Proph_freq ~ Phylotype, data = incomplete_prophages)


# Figure S5

prophage_mash_tree_refs <- read.tree("Intact_prophages_validated_ref_mashtree.dnd")
prophage_mash_tree_refs_rooted <- root(prophage_mash_tree_refs, which(prophage_mash_tree_refs$tip.label == "Burkholderia_KS10"))

tip_colours <- read.csv("Figure_S3_tiplabels.csv",fileEncoding="UTF-8-BOM")
colours <- data.frame(tip_colours$Colour)

prophage_mash_tree_refs_labelled <- ggtree(prophage_mash_tree_refs_rooted,layout = "rectangular")+ ylab("Prophage Mash distance tree")+theme(axis.title=element_text(size=12,face="bold"))+ geom_tiplab(aes(
  subset=(grepl('Ralstonia',label,fixed=TRUE)==TRUE)),colour=colours$tip_colours.Colour, size = 4)+ ggplot2::xlim(0, 2)#geom_tiplab(colour=colours$tip_colours.Colour, size = 5) + ggplot2::xlim(0, 2)


tip_labs <- get_taxa_name(prophage_mash_tree_refs_labelled)
write.csv(tip_labs,"Mashtree_refs__tiplabel.csv")


# Figure 4

world <- map_data('world')

world_map <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="antiquewhite", color="black") +
  coord_quickmap()+theme(panel.background = element_rect(fill = "aliceblue"),panel.border = element_rect(colour = "black", fill=NA, size=1))+ labs(x = "Longitude") + labs(y = "Latitude") +
  theme(axis.text.y = element_text(size = 10, color = "black")) + 
  theme(axis.text.x = element_text(size = 10, color = "black")) + 
  theme(axis.title = element_text( face="bold", size=14))+theme(legend.title = element_text(color="black",size=12,face="bold"))+theme(legend.text = element_text(size=12))+labs(fill="Prophage")

pie_charts <- read.csv("scatterpie_continent_moved2.csv")

colours <- paletteMartin
newcolours <- colours[-c(1,12)]
names(newcolours) <- NULL   

world_map + geom_scatterpie(aes(x=?..long,y=lat,group=region,r=15),data=pie_charts,cols=c("RS551","PE226","RSM3","RSS1","RSS30","RSA1","RsoM1USA","RSY1","RS138","Dina","Unclassified.C","Unclassified.A","Unclassified.F"),alpha=.8)+ 
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

# Fig. S6

prophage_continent_prophage_proportion <- read.csv("Continent_distributions_prophages.csv",fileEncoding="UTF-8-BOM")
prophage_continent_prophage_proportion_gathered <- pivot_longer(prophage_continent_prophage_proportion, cols=c(Africa:Oceania),names_to = "Continents",values_to="Proportion")

prophage_continent_prophage_proportion_gathered$Prophage <- factor(prophage_continent_prophage_proportion_gathered$Prophage,levels = c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

ggplot(data=prophage_continent_prophage_proportion_gathered, aes(x=Prophage,y=Proportion)) + 
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


# Figure 5A

rssc_tree <- read.tree("core_tree.treefile")
rssc_tree_rooted <- midpoint.root(rssc_tree)
tip_labs <- rssc_tree_rooted$tip.label

to_drop <- c("UY031","K60","GMI1000","CMR15","PSI07")
rssc_tree_rooted_dropped <- drop.tip(rssc_tree_rooted,to_drop)

rssc_tree_rooted_dropped_fig <- ggtree(rooted.tree_reduced,layout = "rectangular",ladderize=TRUE)+ 
  theme(plot.margin=margin(0,0,0,0)) + ylab("RSSC phylogeny")+
  theme(axis.title=element_text(size=12,face="bold"))

rssc_tree_rooted_dropped_fig_collapsed <- rssc_tree_rooted_dropped_fig %>%
  collapse(node=201,"max")%>%
  collapse(node=265,"max")

rssc_tree_rooted_dropped_fig_collapsed_highlighted <- rssc_tree_rooted_dropped_fig_collapsed %>%
  + geom_hilight(node=268, fill="steelblue", alpha=.3, extend=1) +
  geom_hilight(node=363, fill="darkgreen", alpha=.3, extend=1)+
  geom_hilight(node=203, fill="yellow", alpha=.3, extend=1)+
  geom_hilight(node=196, fill="red", alpha=.3, extend=1)+
  geom_hilight(node=258, fill="orange", alpha=.3, extend=1)


prophage_presence_matrix <- read.csv("R_pickettii_tree_presence_absence_matrix_validated_intact_prophages.csv",fileEncoding="UTF-8-BOM")
prophage_presence_matrix2 <- prophage_presence_matrix[, -1];
row.names(prophage_presence_matrix2) <- prophage_presence_matrix$Isolate
prophage_presence_matrix3 <- as.matrix(prophage_presence_matrix2)
prophage_presence_matrix3_melted <- melt(prophage_presence_matrix3)

prophage_presence_matrix3_melted$value <- as.factor(prophage_presence_matrix3_melted$value)

prophage_presence_fig <- ggplot(data = melted_matrix, aes(x=Var2, y=Var1, fill=value)) + 
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


rssc_tree_rooted_dropped_fig_collapsed_highlighted + prophage_presence_fig + plot_layout(widths = c(0.4, 1))+plot_layout(guides="collect")


#Figure S7

rssc_tree <- read.tree("core_tree.treefile")
rssc_tree_rooted <- midpoint.root(rssc_tree)
tip_labs <- rssc_tree_rooted$tip.label

to_drop <- c("UY031","K60","GMI1000","CMR15","PSI07")
rssc_tree_rooted_dropped <- drop.tip(rssc_tree_rooted,to_drop)

rssc_tree_rooted_dropped_fig <- ggtree(rooted.tree_reduced,layout = "rectangular",ladderize=TRUE)+ 
  theme(plot.margin=margin(0,0,0,0)) + ylab("RSSC phylogeny")+
  theme(axis.title=element_text(size=12,face="bold"))

incomplete_prophage_presence_matrix <- read.csv("R_pickettii_tree_presence_absence_matrix_incomplete_prophages.csv",fileEncoding="UTF-8-BOM")
incomplete_prophage_presence_matrix2 <- incomplete_prophage_presence_matrix[, -1];
row.names(incomplete_prophage_presence_matrix2) <- incomplete_prophage_presence_matrix$Isolate
incomplete_prophage_presence_matrix3 <- as.matrix(incomplete_prophage_presence_matrix2)
incomplete_prophage_presence_matrix3 <- melt(incomplete_prophage_presence_matrix3)


incomplete_prophage_presence_fig <- ggplot(data = incomplete_prophage_presence_matrix3, aes(x=Var2, y=Var1, fill=value)) + 
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


rssc_tree_rooted_dropped_fig + incomplete_prophage_presence_fig + plot_layout(widths = c(0.4, 1),guides="collect")


## Figure S8

prophage_phylotype_proportion <- read.csv("Phylotype_distributions_prophages_dividephylotypeabundance.csv",fileEncoding="UTF-8-BOM")
prophage_phylotype_proportion_gathered <- pivot_longer(prophage_phylotype_proportion, cols=c(I:IV),names_to = "Phylotypes",values_to="Proportion")

prophage_phylotype_proportion_gathered$Prophage <- factor(prophage_phylotype_proportion_gathered$Prophage,levels = c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

colours <- c("#66CC33","#CCFFFF","#FF9900","#CC0000","#003399","#663300","#FFFF00","#660066","#006600","#999999")

ggplot(data=prophage_phylotype_proportion_gathered, aes(x=Prophage,y=Proportion,fill=Cluster)) + 
  scale_y_continuous(labels = scales::percent, limits=c(0,1))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=colours)+
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


#Figure_S9

BC_and_mash <- read.csv("Phylotype_BC_mash_regression.csv")

FigS9A <- ggplot(data=BC_and_mash, aes(x = Phylotype, y= BC))  +
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+ylab("Average prophage diversity") + 
  xlab("Phylotype") + 
  theme_bw()+ 
  geom_smooth(method = "lm")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title=element_text(size=12,face="bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(fill="Completeness")+ 
  theme(legend.title = element_text(color = "black", size = 10,face="bold"))

FigS9B <- ggplot(data=BC_and_mash, aes(x = Phylotype, y= Mash))  + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  # geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+ylab("Average host diversity") + xlab("Phylotype") + theme_bw()+ geom_smooth(method = "lm")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title=element_text(size=12,face="bold"))+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill="Completeness")  + theme(legend.title = element_text(color = "black", size = 10,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                          legend.text = element_text(color = "black",size=10,face="bold"))                                                                                                                                                                                                                                                                                                                                                                              

FigS9A + FigS9B

kruskal.test(BC_and_mash$Mash~BC_and_mash$Phylotype)
dunnTest(BC_and_mash$Mash~BC_and_mash$Phylotype)
kruskalmc(BC_and_mash$Mash~BC_and_mash$Phylotype,probs=0.05)

kruskal.test(BC_and_mash$BC~BC_and_mash$Phylotype)
dunnTest(BC_and_mash$BC~BC_and_mash$Phylotype)
kruskalmc(BC_and_mash$BC~BC_and_mash$Phylotype,probs=0.05)


# Figure 5B

# Make Bray-Curtis distance matrices
phyloI_presence_matrix <- read.csv("PhyloI_presenceabsence.csv",fileEncoding="UTF-8-BOM")
rownames(phyloI_presence_matrix) <- phyloI_presence_matrix[, 1];
phyloI_presence_matrix2 <- phyloI_presence_matrix[, -1];
phyloI_presence_matrix3 <- as.matrix(phyloI_presence_matrix2)

phyloI_bc <- vegdist(phyloI_presence_matrix3, method="bray")
phyloI_bc[is.na(phyloI_bc)] <- 0
phyloI_bc_matrix <- as.matrix(phyloI_bc)
write.csv(phyloI_bc_matrix, "PhyloI_BCmatrix.csv")


phyloIIA_presence_matrix <- read.csv("PhyloIIA_presenceabsence.csv")
rownames(phyloIIA_presence_matrix) <- phyloIIA_presence_matrix[, 1];
phyloIIA_presence_matrix2 <- phyloIIA_presence_matrix[, -1];
phyloIIA_presence_matrix3 <- as.matrix(phyloIIA_presence_matrix2)

phyloIIA_bc <- vegdist(phyloIIA_presence_matrix3, method="bray")
phyloIIA_bc[is.na(phyloIIA_bc)] <- 0
phyloIIA_bc_matrix <- as.matrix(phyloIIA_bc)
write.csv(PhyloIIA_BC, "PhyloIIA_BCmatrix.csv")


phyloIIB_presence_matrix <- read.csv("PhyloIIB_presenceabsence.csv")
rownames(phyloIIB_presence_matrix) <- phyloIIB_presence_matrix[, 1];
phyloIIB_presence_matrix2 <- phyloIIB_presence_matrix[, -1];
phyloIIB_presence_matrix3 <- as.matrix(phyloIIB_presence_matrix2)

phyloIIB_bc <- vegdist(phyloIIB_presence_matrix3, method="bray")
phyloIIB_bc[is.na(phyloIIB_bc)] <- 0
phyloIIB_bc_matrix <- as.matrix(phyloIIB_bc)
write.csv(PhyloIIB_BC, "PhyloIIB_BCmatrix.csv")


# Take averages of each row and add to spreadsheet with host mash distances
bc_mash_regression <- read.csv("Phylotype_BC_mash_regression.csv")
lm2 <- lme(BC ~ log(Mash),data=bc_mash_regression, random=~1|Phylotype)

coef(lm2)
par(mfrow=c(2,2))
hist(lm2$residuals)
qqnorm(lm2$residuals)
qqline(lm2$residuals)
plot(lm2$fitted,lm2$residuals)
summary(lm2)
anova(lm2)
rsquared(lm2)

ggplot(bc_mash_regression,aes(x=log(Mash),y=BC))+ 
  geom_point(aes(fill=Phylotype),colour="black",pch=21, size=3)+
  xlim(c(-6.5,-3.5))+ylim(c(0.125,1))+
  labs(y="Average prophage dissimilarity",x="Average <i>R. solanacearum</i> dissimilarity (log<sub>10</sub>)")+
  theme_bw()+
  theme(axis.title.x = element_markdown())+
  geom_segment(aes(x=-6.06,xend=-3.888,y=0.125,yend=1),size=1)+
  theme(legend.position = c(0.15,0.83))+
  theme(axis.title=element_text(size=13,face="bold"))+
  theme(legend.title = element_text(colour="black", size=16,face="bold"),legend.text = element_text(size=10))


## Figure 5C

# All isolate phylogeny
rssc_tree <- read.tree("core_tree.treefile")
rssc_tree_rooted <- midpoint.root(rssc_tree)

## PACo analysis of phylotype I clade

# Extract phylotype I clade
rssc_tree_rooted_phyloI <- extract.clade(rssc_tree_rooted,209)
rssc_tree_rooted_phyloI_drop <- drop.tip(rssc_tree_rooted_phyloI,"GMI1000")

phyloI_presence_matrix <- read.csv("PhyloI_presenceabsence.csv",fileEncoding="UTF-8-BOM")
rownames(phyloI_presence_matrix) <- phyloI_presence_matrix[, 1];
phyloI_presence_matrix2 <- phyloI_presence_matrix[, -1];
phyloI_presence_matrix3 <- as.matrix(phyloI_presence_matrix2)

phyloI_bc <- vegdist(phyloI_presence_matrix3, method="bray")
phyloI_bc[is.na(phyloI_bc)] <- 0
phyloI_bc_matrix <- as.matrix(phyloI_bc)

phyloI_UPGMA <- upgma(phyloI_bc_matrix)

# Start running PACo
host.D <- cophenetic(rssc_tree_rooted_phyloI_drop)
BC.D <- cophenetic(phyloI_UPGMA)

HP <- read.csv("PACo_binarymatrix_phyloI.csv",fileEncoding="UTF-8-BOM")
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


## PACo analysis of phylotype IIA clade
rssc_tree_rooted_phyloIIA <- extract.clade(rssc_tree_rooted,372)
rssc_tree_rooted_phyloIIA_drop <- drop.tip(rssc_tree_rooted_phyloIIA,"K60")

phyloIIA_presence_matrix <- read.csv("PhyloIIA_presenceabsence.csv")
rownames(phyloIIA_presence_matrix) <- phyloIIA_presence_matrix[, 1];
phyloIIA_presence_matrix2 <- phyloIIA_presence_matrix[, -1];
phyloIIA_presence_matrix3 <- as.matrix(phyloIIA_presence_matrix2)

phyloIIA_bc <- vegdist(phyloIIA_presence_matrix3, method="bray")
phyloIIA_bc[is.na(phyloIIA_bc)] <- 0
phyloIIA_bc_matrix <- as.matrix(phyloIIA_bc)

phyloIIA_UPGMA <- upgma(phyloIIA_bc_matrix)

# Start running PACo
host.D <- cophenetic(rssc_tree_rooted_phyloIIA_drop)
BC.D <- cophenetic(phyloIIA_UPGMA)

HP <- read.csv("PACo_binarymatrix_phyloIIA.csv",fileEncoding="UTF-8-BOM")
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

## PACo analysis of phylotype IIB clade

rssc_tree_rooted_phyloIIB <- extract.clade(rssc_tree_rooted,276)
rssc_tree_rooted_phyloIIB_drop <- drop.tip(rssc_tree_rooted_phyloIIB,"UY031")

phyloIIB_presence_matrix <- read.csv("PhyloIIB_presenceabsence.csv")
rownames(phyloIIB_presence_matrix) <- phyloIIB_presence_matrix[, 1];
phyloIIB_presence_matrix2 <- phyloIIB_presence_matrix[, -1];
phyloIIB_presence_matrix3 <- as.matrix(phyloIIB_presence_matrix2)

phyloIIB_bc <- vegdist(phyloIIB_presence_matrix3, method="bray")
phyloIIB_bc[is.na(phyloIIB_bc)] <- 0
phyloIIB_bc_matrix <- as.matrix(phyloIIB_bc)

phyloIIB_UPGMA <- upgma(phyloIIB_bc_matrix)

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

tip <- c("K60","GMI1000","UY031","CMR15","PSI07")
rssc_tree_rooted_dropped <- drop.tip(rssc_tree_rooted,tip)

# Use presence/absence matrix with reference isolates removed
all_presence_matrix <- read.csv("R_pickettii_tree_presence_absence_matrix_forcongruence.csv")
rownames(all_presence_matrix) <- all_presence_matrix[, 1];
all_presence_matrix2 <- all_presence_matrix[, -1];
all_presence_matrix3 <- as.matrix(all_presence_matrix2)

all_bc <- vegdist(all_presence_matrix3, method="bray")
all_bc[is.na(all_bc)] <- 0
all_bc_matrix <- as.matrix(all_bc)

all_bc_UPGMA <- upgma(all_bc_matrix)

# Start running PACo
host.D <- cophenetic(rssc_tree_rooted_dropped)
BC.D <- cophenetic(all_bc_UPGMA)

HP <- read.csv("PACo_binarymatrix_All.csv",fileEncoding="UTF-8-BOM")
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

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

# Tanglegram
rssc_tree_rooted_dropped2 <- rotateNodes(rssc_tree_rooted_dropped,"all")
rssc_tree_rooted_dropped3 <- untangle(rssc_tree_rooted_dropped2)

cophyloplot(rooted.tree.drop3, AllBC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
                         lwd=1, col='steelblue', length.line=0, gap=0, space=70,rotate=TRUE)+coord_flip()




# Figure 6

# VIGA results
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

CAZy + RalstoT3E + Phi_base +plot_layout(ncol=1)

Prophage_phylo <- read.csv("Prophage_phylotype_heatmap.csv",fileEncoding="UTF-8-BOM")
Prophage_phylo$Phage_hits <- factor(Prophage_phylo$Phage_hits, levels=c("Ralstonia phage Rs551","Ralstonia phage RSM3","Unclassified B","Ralstonia phage RSS1","Ralstonia phage RSS30","Ralstonia phage RSS-TH1","Ralstonia phage PE226","Ralstonia phage phiRSA1","Ralstonia phage RsoM1USA","Ralstonia phage RSY1","Unclassified D","Ralstonia phage RS138","Unclassified G","Ralstonia phage Dina","Unclassified C","Unclassified H","Unclassified A","Unclassified F","Unclassified I","Unclassified J"))

colour_phylo <- c("khaki","lightsalmon","darkorchid3","forestgreen","gray80")

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

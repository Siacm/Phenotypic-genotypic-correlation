library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library("readr")
library("data.table")
library("dplyr")
library(openxlsx)
library(tibble)
library(ggrepel)

# Classic palette BuPu, with 4 colors
coul = brewer.pal(4, "BuPu") 

# I can add more tones to this palette :
coul = colorRampPalette(coul)(25)

#########azithro pie chart##########
azith_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "azitrgenes_count", colNames = TRUE )
colnames(azith_rgenes)<-c("Gene", "Count")
azith_rgenes$Gene <- gsub(":", "", azith_rgenes$Gene)
azith_rgenes<-azith_rgenes[,c(1,2)]

azith_rgenes$fraction = azith_rgenes$Count / sum(azith_rgenes$Count)

azith_rgenes$ymax = cumsum(azith_rgenes$fraction)
azith_rgenes$ymin = c(0, head(azith_rgenes$ymax, n=-1))
colnames(azith_rgenes$Gene)

azith_rgenes$Gene <- factor(azith_rgenes$Gene, levels = rev(as.character(azith_rgenes$Gene)))

p1 = ggplot(azith_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

azith_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(6,"Blues"))(4), name = "Gene", labels = rev(paste(percent(azith_rgenes$fraction)," ", azith_rgenes$Gene," ","(n=",azith_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size=20, face = "bold")) +labs(x = NULL, y = NULL, title = "Azithromycin Resistance Genes      Total: 20 individual genes ")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))


#azith_pie +  geom_label_repel(aes(label=paste(percent(azith_rgenes$fraction),"\n", azith_rgenesc$Gene,"\n","n=",azith_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 2)
#azith_pie +  geom_label_repel(aes(label=paste(percent(azith_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=8)

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/azithro_pie.png', plot=azith_pie, width = 14, height = 14, dpi = 100) ###width, height changes the sizes of the labels too.



#########cef pie chart##########
cef_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "cefrgenes_count", colNames = TRUE )
colnames(cef_rgenes)<-c("Gene", "Count")
cef_rgenes$Gene <- gsub(" :", "", cef_rgenes$Gene)
cef_rgenes<-cef_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cef_rgenes$fraction = cef_rgenes$Count / sum(cef_rgenes$Count)

cef_rgenes$ymax = cumsum(cef_rgenes$fraction)
cef_rgenes$ymin = c(0, head(cef_rgenes$ymax, n=-1))
colnames(cef_rgenes$Gene)

cef_rgenes$Gene <- factor(cef_rgenes$Gene, levels = rev(as.character(cef_rgenes$Gene)))

p1 = ggplot(cef_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

cef_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Greens"))(10), name = "Gene", labels = rev(paste(percent(cef_rgenes$fraction)," ", cef_rgenes$Gene," ","(n=",cef_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size=20, face = "bold")) +labs(x = NULL, y = NULL, title = "Top 10 Cefotaxime Resistance Genes       Total: 395 individual genes")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))

#cef_pie +  geom_label_repel(aes(label=paste(percent(cef_rgenes$fraction),"\n", cef_rgenesc$Gene,"\n","n=",cef_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE,nudge_x = 2)
#cef_pie +  geom_label_repel(aes(label=paste(percent(cef_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=5, direction="x", fill="white")

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/individual_cef_pie.png', plot=cef_pie, width =14, height = 10, dpi = 100) ###changes the sizes of the labels too.


#########rpattern genes cip pie pheno resis chart##########
cip_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "ciprpattern_count", colNames = TRUE )
colnames(cip_rgenes)<-c("Gene", "Count")
cip_rgenes$Gene <- gsub(" :", "", cip_rgenes$Gene)
cip_rgenes<-cip_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cip_rgenes$fraction = cip_rgenes$Count / sum(cip_rgenes$Count)

cip_rgenes$ymax = cumsum(cip_rgenes$fraction)
cip_rgenes$ymin = c(0, head(cip_rgenes$ymax, n=-1))
colnames(cip_rgenes$Gene)

cip_rgenes$Gene <- factor(cip_rgenes$Gene, levels = rev(as.character(cip_rgenes$Gene)))

p1 = ggplot(cip_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

cip_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Oranges"))(10), name = "Gene", labels = rev(paste(percent(cip_rgenes$fraction)," ", cip_rgenes$Gene," ","(n=",cip_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size=20, face = "bold")) +labs(x = NULL, y = NULL, title = "Top 10 Ciprofloxacin Resistance Genes Total: X individual genes")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))

#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction),"\n", cip_rgenesc$Gene,"\n","n=",cip_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE,nudge_x = 2)
#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=8)

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/cipresisonlypattern_pie.png', plot=cip_pie, width = 20, height = 10, dpi = 100)


#########rpattern genes cip pie pheno resis chart##########

cip_rgenes$fraction = cip_rgenes$Count / sum(cip_rgenes$Count)

cip_rgenes$ymax = cumsum(cip_rgenes$fraction)
cip_rgenes$ymin = c(0, head(cip_rgenes$ymax, n=-1))
colnames(cip_rgenes$Gene)

cip_rgenes$Gene <- factor(cip_rgenes$Gene, levels = rev(as.character(cip_rgenes$Gene)))

p1 = ggplot(cip_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

cip_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Oranges"))(10), name = "Gene", labels = rev(paste(percent(cip_rgenes$fraction)," ", cip_rgenes$Gene," ","(n=",cip_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size=20, face = "bold")) +labs(x = NULL, y = NULL, title = "Top 10 Ciprofloxacin Resistance Genes Total:  X individual genes")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))


#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction),"\n", cip_rgenesc$Gene,"\n","n=",cip_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE,nudge_x = 2)
#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=8)

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/cipironlypattern_pie.png', plot=cip_pie, width = 20, height = 10, dpi = 100)



############ADDITIONAL PIE CHARTS##################

#########individual genes cip pie pheno resis chart##########
azith_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "ciprgenes_count", colNames = TRUE )
colnames(cip_rgenes)<-c("Gene", "Count")
cip_rgenes$Gene <- gsub(" :", "", cip_rgenes$Gene)
cip_rgenes<-cip_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cip_rgenes$fraction = cip_rgenes$Count / sum(cip_rgenes$Count)

cip_rgenes$ymax = cumsum(cip_rgenes$fraction)
cip_rgenes$ymin = c(0, head(cip_rgenes$ymax, n=-1))
colnames(cip_rgenes$Gene)

cip_rgenes$Gene <- factor(cip_rgenes$Gene, levels = rev(as.character(cip_rgenes$Gene)))

p1 = ggplot(cip_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

cip_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Oranges"))(10), name = "Gene", labels = rev(paste(percent(cip_rgenes$fraction)," ", cip_rgenes$Gene," ","(n=",cip_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size=20, face = "bold")) +labs(x = NULL, y = NULL, title = "Top 10 Ciprofloxacin Resistance Genes Total: X individual genes")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))

#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction),"\n", cip_rgenesc$Gene,"\n","n=",cip_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE,nudge_x = 2)
#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=8)

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/individualcipresisonly_pie.png', plot=cip_pie, width = 20, height = 10, dpi = 100)

#### individual genes cip pheno resis/intermediate pie chart##########
cip_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/perl_flago_cov/ciprofloxacin_tallied_results.xlsx", sheet = "count_phenoresisinter", colNames = FALSE )
colnames(cip_rgenes)<-c("Gene", "Count")
cip_rgenes$Gene <- gsub(" :", "", cip_rgenes$Gene)
cip_rgenes<-cip_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cip_rgenes$fraction = cip_rgenes$Count / sum(cip_rgenes$Count)

cip_rgenes$ymax = cumsum(cip_rgenes$fraction)
cip_rgenes$ymin = c(0, head(cip_rgenes$ymax, n=-1))
colnames(cip_rgenes$Gene)

cip_rgenes$Gene <- factor(cip_rgenes$Gene, levels = rev(as.character(cip_rgenes$Gene)))

p1 = ggplot(cip_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

cip_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Oranges"))(10), name = "Gene", labels = rev(paste(percent(cip_rgenes$fraction)," ", cip_rgenes$Gene," ","(n=",cip_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size=20, face = "bold")) +labs(x = NULL, y = NULL, title = "Top 10 Ciprofloxacin Resistance Genes Total:  X individual genes")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))

#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction),"\n", cip_rgenesc$Gene,"\n","n=",cip_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE,nudge_x = 2)
#cip_pie +  geom_label_repel(aes(label=paste(percent(cip_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=8)

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/cipresisinter_pie.png', plot=cip_pie, width = 20, height = 10, dpi = 100)


#########azithro_novel pie chart##########
azithnovel_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "novelazitrgenes_count", colNames = TRUE )
colnames(azithnovel_rgenes)<-c("Gene", "Count")
azithnovel_rgenes$Gene <- gsub(":", "", azithnovel_rgenes$Gene)
azithnovel_rgenes<-azithnovel_rgenes[,c(1,2)]

azithnovel_rgenes$fraction = azithnovel_rgenes$Count / sum(azithnovel_rgenes$Count)

azithnovel_rgenes$ymax = cumsum(azithnovel_rgenes$fraction)
azithnovel_rgenes$ymin = c(0, head(azithnovel_rgenes$ymax, n=-1))
colnames(azithnovel_rgenes$Gene)

azithnovel_rgenes$Gene <- factor(azithnovel_rgenes$Gene, levels = rev(as.character(azithnovel_rgenes$Gene)))

p1 = ggplot(azithnovel_rgenes, aes(fill=Gene, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='blue') +
  coord_polar(theta="y") +
  xlim(c(1, 4)) 

azithnovel_pie<-p1 +   guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(4,"Blues"))(5), name = "Gene", labels = rev(paste(percent(azithnovel_rgenes$fraction)," ", azithnovel_rgenes$Gene," ","(n=",azithnovel_rgenes$Count,")")))+ blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle("") +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(legend.title = element_text(size=20, face="bold")) +
  theme(legend.text = element_text(size = 20, face = "bold")) +labs(x = NULL, y = NULL, title = "Azithromycin Resistance Genes      Total: 30 individual genes ")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=22))


#azithnovel_pie +  geom_label_repel(aes(label=paste(percent(azithnovel_rgenes$fraction),"\n", azithnovel_rgenesc$Gene,"\n","n=",azithnovel_rgenesc$volume),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE,nudge_x = 2)
#azithnovel_pie +  geom_label_repel(aes(label=paste(percent(azithnovel_rgenes$fraction)),x=4,y=(ymin+ymax)/2),inherit.aes = TRUE, show.legend = FALSE, nudge_x = 1, size=8)

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/azithronovel_pie.png', plot=azithnovel_pie, width = 14, height = 14, dpi = 100)



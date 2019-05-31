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
getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx")
azith_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "azitrgenes_count", colNames = TRUE )
#colnames(azith_rgenes)<-c("Gene", "Count")
azith_rgenes$Gene <- gsub(":", "", azith_rgenes$Gene)
azith_rgenes<-azith_rgenes[,c(1,2)]

azith_rgenesc <- azith_rgenes %>% ###arranges the piechart by largest to smallest part
  group_by(Gene) %>%
  summarise(volume = sum(Count)) %>%
  mutate(share=volume/sum(volume)) %>%
  arrange(desc(volume))

#Pos <- c(cumsum(360*azith_rgenesc$share)-(360*azith_rgenesc$share/2))
#Pos <- c(ifelse(Pos<=180,Pos,Pos-180))
#azith_rgenesc$Pos<-Pos
azith_rgenesc$Gene <- factor(azith_rgenesc$Gene, levels = rev(as.character(azith_rgenesc$Gene)))

azithro_labels <- tibble(x.breaks = seq(1, 1.5, length.out = length(azith_rgenesc$volume)), ###needed for the box labels
                    y.breaks = cumsum(azith_rgenesc$volume) - azith_rgenesc$volume/2,
                    labels = paste(percent(azith_rgenesc$share),"\n", azith_rgenesc$Gene,"\n","n=",azith_rgenesc$volume),
                    Descripcion = azith_rgenesc$Gene)

azith_pie <- ggplot(azith_rgenesc, aes(x=1, y=volume, fill=Gene)) +
  geom_bar(stat="identity", color="black", width=1) +
  coord_polar(theta='y') +
  geom_label_repel(data = azithro_labels, aes(x = x.breaks, y = y.breaks, ###gives the box labels
                                         label = labels, fill = azith_rgenesc$Gene),
                   label.padding = unit(0.1, "lines"),
                   size = 5,
                   color="black",
                   show.legend = FALSE,
                   inherit.aes = FALSE)+
  #geom_text(aes(x=2.1,
  #              label=paste0(percent(azith_rgenesc$share)," ",azith_rgenesc$Gene," (n=",azith_rgenesc$volume,")"),
  #              angle=90-azith_rgenesc$Pos), position = position_stack(vjust = 0.5), color="red", size=4.5) +  #####to change the rotation, change X in angle=X-genes, change the position of labels at the beginning with x="length_you_want"
  labs(x = NULL, y = NULL, title = "Azithromycin Resistance Genes") +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(6,"Blues"))(4))+
  theme(axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20), legend.text=element_text(size=20)) +
  scale_x_discrete(limits=c(0, 1))
azith_pie
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/azithro_pie.png', plot=azith_pie, width = 14, height = 14, dpi = 100) ###width, height changes the sizes of the labels too.



#########cef pie chart##########
cef_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "cefrgenes_count", colNames = TRUE )
#colnames(cef_rgenes)<-c("Gene", "Count")
cef_rgenes$Gene <- gsub(" :", "", cef_rgenes$Gene)
cef_rgenes<-cef_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cef_rgenesc <- cef_rgenes %>% ###arranges the piechart by largest to smallest part
  group_by(Gene) %>%
  summarise(volume = sum(Count)) %>%
  mutate(share=volume/sum(volume)) %>%
  arrange(desc(volume))

#Pos <- c(cumsum(360*cef_rgenesc$share)-(360*cef_rgenesc$share/2))
#Pos <- c(ifelse(Pos<=180,Pos,Pos-180))
#cef_rgenesc$Pos<-Pos
cef_rgenesc$Gene <- factor(cef_rgenesc$Gene, levels = rev(as.character(cef_rgenesc$Gene)))

cef_labels <- tibble(x.breaks = seq(1, 1.5, length.out = length(cef_rgenesc$volume)), 
                    y.breaks = cumsum(cef_rgenesc$volume) - cef_rgenesc$volume/2,
                    labels = paste(percent(cef_rgenesc$share),"\n", cef_rgenesc$Gene,"\n","n=",cef_rgenesc$volume),
                    Descripcion = cef_rgenesc$Gene)

cef_pie <- ggplot(cef_rgenesc, aes(x=1, y=volume, fill=Gene)) +
  geom_bar(stat="identity", color="black", width=1) +
  coord_polar(theta='y') +
  geom_label_repel(data = cef_labels, aes(x = x.breaks, y = y.breaks, 
                                         label = labels, fill = cef_rgenesc$Gene),
                   label.padding = unit(0.1, "lines"),
                   size = 3.5,
                   color="black",
                   show.legend = FALSE,
                   inherit.aes = FALSE, nudge_x = 1)+
  #geom_text(aes(x=2.1,
  #              label=paste0(percent(cef_rgenesc$share)," ",cef_rgenesc$Gene," (n=",cef_rgenesc$volume,")"),
  #              angle=95-cef_rgenesc$Pos), position = position_stack(vjust = 0.5), color="red", size=4) +     #####to change the rotation, change X in angle=X-genes, change the position of labels at the beginning with x="length_you_want"
  guides(fill=guide_legend(override.aes=list(colour=NA)), reverse=T) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Greens"))(10))+
  labs(x = NULL, y = NULL, title = "Top 10 Cefotaxime Resistance Genes") +
  theme(axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20), legend.text=element_text(size=20)) +
  scale_x_discrete(limits=c(0, 1))
cef_pie
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/individual_cef_pie.png', plot=cef_pie, width =14, height = 10, dpi = 100) ###changes the sizes of the labels too.


#########rpattern genes cip pie pheno resis chart##########
cip_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "ciprpattern_count", colNames = TRUE )
#colnames(cip_rgenes)<-c("Gene", "Count")
cip_rgenes$Gene <- gsub(" :", "", cip_rgenes$Gene)
cip_rgenes<-cip_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cip_rgenesc <- cip_rgenes %>% ###arranges the piechart by largest to smallest part
  group_by(Gene) %>%
  summarise(volume = sum(Count)) %>%
  mutate(share=volume/sum(volume)) %>%
  arrange(desc(volume))

#Pos <- c(cumsum(360*cip_rgenesc$share)-(360*cip_rgenesc$share/2))
#Pos <- c(ifelse(Pos<=180,Pos,Pos-180))
#cip_rgenesc$Pos<-Pos
cip_rgenesc$Gene <- factor(cip_rgenesc$Gene, levels = rev(as.character(cip_rgenesc$Gene)))

cip_labels <- tibble(x.breaks = seq(1, 1.5, length.out = length(cip_rgenesc$volume)), 
                     y.breaks = cumsum(cip_rgenesc$volume) - cip_rgenesc$volume/2,
                     labels = paste(percent(cip_rgenesc$share),"\n", cip_rgenesc$Gene,"\n","n=",cip_rgenesc$volume),
                     Descripcion = cip_rgenesc$Gene)

cip_pie <- ggplot(cip_rgenesc, aes(x=1, y=volume, fill=Gene)) +
  geom_bar(stat="identity", color="black", width=1) +
  coord_polar(theta='y') +
  geom_label_repel(data = cip_labels, aes(x = x.breaks, y = y.breaks, 
                                          label = labels, fill = cip_rgenesc$Gene),
                   label.padding = unit(0.1, "lines"),
                   size = 3.5,
                   color="black",
                   show.legend = FALSE,
                   inherit.aes = FALSE)+
  #geom_text(aes(x=2.1,
  #              label=paste0(percent(cip_rgenesc$share)," ",cip_rgenesc$Gene," (n=",cip_rgenesc$volume,")"),
  #              angle=95-cip_rgenesc$Pos), position = position_stack(vjust = 0.5), color="red", size=4) +     #####to change the rotation, change X in angle=X-genes, change the position of labels at the beginning with x="length_you_want"
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Oranges"))(10))+
  labs(x = NULL, y = NULL, title = "Top 10 Ciprofloxacin Resistance Genes") +
  theme(axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20), legend.text=element_text(size=20)) +
  scale_x_discrete(limits=c(0, 1))
cip_pie
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/cipresisonlypattern_pie.png', plot=cip_pie, width = 20, height = 10, dpi = 100)


#########individual genes cip pie pheno resis chart##########
cip_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "ciprgenes_count", colNames = TRUE )
#colnames(cip_rgenes)<-c("Gene", "Count")
cip_rgenes$Gene <- gsub(" :", "", cip_rgenes$Gene)
cip_rgenes<-cip_rgenes[1:10,c(1,2)]  #####to get only top 10 genes

cip_rgenesc <- cip_rgenes %>% ###arranges the piechart by largest to smallest part
  group_by(Gene) %>%
  summarise(volume = sum(Count)) %>%
  mutate(share=volume/sum(volume)) %>%
  arrange(desc(volume))

#Pos <- c(cumsum(360*cip_rgenesc$share)-(360*cip_rgenesc$share/2))
#Pos <- c(ifelse(Pos<=180,Pos,Pos-180))
#cip_rgenesc$Pos<-Pos
cip_rgenesc$Gene <- factor(cip_rgenesc$Gene, levels = rev(as.character(cip_rgenesc$Gene)))

cip_labels <- tibble(x.breaks = seq(1, 1.5, length.out = length(cip_rgenesc$volume)), 
                     y.breaks = cumsum(cip_rgenesc$volume) - cip_rgenesc$volume/2,
                     labels = paste(percent(cip_rgenesc$share),"\n", cip_rgenesc$Gene,"\n","n=",cip_rgenesc$volume),
                     Descripcion = cip_rgenesc$Gene)

cip_pie <- ggplot(cip_rgenesc, aes(x=1, y=volume, fill=Gene)) +
  geom_bar(stat="identity", color="black", width=1) +
  coord_polar(theta='y') +
  geom_label_repel(data = cip_labels, aes(x = x.breaks, y = y.breaks, 
                                          label = labels, fill = cip_rgenesc$Gene),
                   label.padding = unit(0.1, "lines"),
                   size = 3.5,
                   color="black",
                   show.legend = FALSE,
                   inherit.aes = FALSE)+
  #geom_text(aes(x=2.1,
  #              label=paste0(percent(cip_rgenesc$share)," ",cip_rgenesc$Gene," (n=",cip_rgenesc$volume,")"),
  #              angle=95-cip_rgenesc$Pos), position = position_stack(vjust = 0.5), color="red", size=4) +     #####to change the rotation, change X in angle=X-genes, change the position of labels at the beginning with x="length_you_want"
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Oranges"))(10))+
  labs(x = NULL, y = NULL, title = "Top 10 Ciprofloxacin Resistance Genes") +
  theme(axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20), legend.text=element_text(size=20)) +
  scale_x_discrete(limits=c(0, 1))
cip_pie
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/cipironlypattern_pie.png', plot=cip_pie, width = 20, height = 10, dpi = 100)



############ADDITIONAL PIE CHARTS##################

#########azithro_novel pie chart##########
azithnovel_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "novelazitrgenes_count", colNames = TRUE )
#colnames(azithnovel_rgenes)<-c("Gene", "Count")
azithnovel_rgenes$Gene <- gsub(":", "", azithnovel_rgenes$Gene)
azithnovel_rgenes<-azithnovel_rgenes[,c(1,2)]

azithnovel_rgenesc <- azithnovel_rgenes %>% ###arranges the piechart by largest to smallest part
  group_by(Gene) %>%
  summarise(volume = sum(Count)) %>%
  mutate(share=volume/sum(volume)) %>%
  arrange(desc(volume))

#Pos <- c(cumsum(360*azithnovel_rgenesc$share)-(360*azithnovel_rgenesc$share/2))
#Pos <- c(ifelse(Pos<=180,Pos,Pos-180))
#azithnovel_rgenesc$Pos<-Pos
azithnovel_rgenesc$Gene <- factor(azithnovel_rgenesc$Gene, levels = rev(as.character(azithnovel_rgenesc$Gene)))

azithronovel_labels <- tibble(x.breaks = seq(1, 1.5, length.out = length(azithnovel_rgenesc$volume)),
                    y.breaks = cumsum(azithnovel_rgenesc$volume) - azithnovel_rgenesc$volume/2,
                    labels = paste(percent(azithnovel_rgenesc$share),"\n", azithnovel_rgenesc$Gene,"\n","n=",azithnovel_rgenesc$volume),
                    Descripcion = azithnovel_rgenesc$Gene)

azithnovel_pie <- ggplot(azithnovel_rgenesc, aes(x=1, y=volume, fill=Gene)) +
  geom_bar(stat="identity", color="black", width=1) +
  coord_polar(theta='y') +
  geom_label_repel(data = azithronovel_labels, aes(x = x.breaks, y = y.breaks, 
                                         label = labels, fill = azithnovel_rgenesc$Gene),
                   label.padding = unit(0.1, "lines"),
                   size = 5,
                   color="black",
                   show.legend = FALSE,
                   inherit.aes = FALSE)+
  #geom_text(aes(x=2,
  #              label=paste0(percent(azithnovel_rgenesc$share)," ",azithnovel_rgenesc$Gene," (n=",azithnovel_rgenesc$volume,")"),
  #              angle=90-azithnovel_rgenesc$Pos), position = position_stack(vjust = 0.5), color="red", size=4) +  #####to change the rotation, change X in angle=X-genes, change the position of labels at the beginning with x="length_you_want"
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(4,"Blues"))(5))+
  labs(x = NULL, y = NULL, title = "Azithromycin Resistance Genes") +
  theme(axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20), legend.text=element_text(size=20)) +
  scale_x_discrete(limits=c(0, 1))
azithnovel_pie
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/azithronovel_pie.png', plot=azithnovel_pie, width = 14, height = 1, dpi = 100)
library(ggplot2)
library(openxlsx)
library(RColorBrewer)


############top 10 genetically diverse serovars (genetically different AMR genes, could be from same abx class)
Rserovar_data<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "all_rsero_abx_rgene_count",colNames = T )
Rserovar_data$Gene <- gsub(" :", "", Rserovar_data$Gene)

# Stacked
ggplot(Rserovar_data, aes(fill=Antibiotic, y=Count, x=ST)) + 
  geom_bar( stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#############top 10 resistant serovars (highest number of resistance to an overall abx classes)
getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_antibiotic_counts070519.xlsx")
Rserovar_data<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_antibiotic_counts070519.xlsx", sheet = "top10_rsero_abx",colNames = T )
#Rserovar_data$Gene <- gsub(" :", "", Rserovar_data$Gene)

# Stacked
ggplot(Rserovar_data, aes(fill=Antibiotic, y=Count, x=ST)) + 
  geom_bar( stat="identity")+
  scale_fill_manual("legend", values = c("Azithromycin" = "#ffffb3","Ampicillin" = "#8dd3c7", "Cefotaxime" = "turquoise", "Meropenem" = "seagreen",
                                         "Chloramphenicol" = "#fdb462",  "Streptomycin" = "#bebada", "Tetracycline" = "#80b1d3", "Ciprofloxacin" = "#d9d9d9","Sulphathiozole" = "#fb8072"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################top 10 cef-resistant serovars (highest number of resistance to cef - with additional resistant genes)
getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Serovar_counts.xlsx")
Rcef_data<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Serovar_counts.xlsx", sheet = "cef_rgenes",colNames = T )
Rcef_data$Gene <- gsub(" :", "", Rcef_data$Gene)
Rcef_data<-Rcef_data[,c(1:4)]

# Stacked
ggplot(Rcef_data, aes(fill=Antibiotic, y=Count, x=reorder(ST, -Count, sum))) + 
  geom_bar( stat="identity", width = 0.8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Top 10 Serovars Resistant to Cefoxtaxime", x="Serovars", fill="")+
  scale_fill_manual("legend", values = c("Azithromycin" = "#ffffb3","Ampicillin" = "#8dd3c7", "Cefotaxime" = "turquoise", "Meropenem" = "seagreen",
                                         "Chloramphenicol" = "#fdb462",  "Streptomycin" = "#bebada", "Tetracycline" = "#80b1d3", "Ciprofloxacin" = "#d9d9d9","Sulphathiozole" = "#fb8072"))+
  theme(title = element_text(size = 18), legend.text=element_text(size=18), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15))+
  # zoom in to last 3 bars on the axis
  scale_x_discrete(limits = tail(levels(reorder(Rcef_data$ST, Rcef_data$Count, FUN = sum)), 10)) + #looks at top 10
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5)) #centre the title
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Rcef_serovar.png', width = 16, height = 9, dpi = 100)


################top 10 cefONLY-resistant serovars (highest number of resistance to cef - with additional resistant genes)
getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Serovar_counts.xlsx")
Rcef_data<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Serovar_counts.xlsx", sheet = "onlycef_rgenes",colNames = T )
Rcef_data$Gene <- gsub(" :", "", Rcef_data$Gene)
Rcef_data<-Rcef_data[,c(1:4)]

# Stacked
ggplot(Rcef_data, aes(fill=Gene, y=Count, x=reorder(ST, -Count, sum))) + 
  geom_bar( stat="identity", width = 0.8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Top 10 Serovars Resistant to Cefoxtaxime", x="Serovars", fill="")+
  theme(title = element_text(size = 18), legend.text=element_text(size=18), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15))+
  # zoom in to last 3 bars on the axis
  scale_x_discrete(limits = tail(levels(reorder(Rcef_data$ST, Rcef_data$Count, FUN = sum)), 10)) + #looks at top 10
  scale_fill_manual(values = rev(colorRampPalette(brewer.pal(12,"PRGn"))(12)), guide=guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5)) #centre the title
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/onlyRcef_serovar.png', width = 16, height = 9, dpi = 100)



##############top 10 overall-resistant serovars looking at all antiboitcs
getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_antibiotic_counts070519.xlsx")
Rallabx_data<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_antibiotic_counts070519.xlsx", sheet = "rsero_abxclass>2_count",colNames = T )
Rallabx_data<-Rallabx_data[,c(1:3)]

# Stacked
ggplot(Rallabx_data, aes(fill=factor(Antibiotic, levels=c("Ampicillin", "Cefotaxime","Meropenem", "Azithromycin", "Ciprofloxacin", "Chloramphenicol", "Streptomycin", "Sulphathiozole", "Tetracycline")), y=Count, x=reorder(ST, -Count, sum), label=Count)) + 
  geom_bar( stat="identity", width = 0.8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Top 10 Multidrug Resistant Serovars in this Study", x="Serovars", fill="")+
  scale_fill_manual("legend", values = c("Azithromycin" = "#ffffb3","Ampicillin" = "#8dd3c7", "Cefotaxime" = "turquoise", "Meropenem" = "seagreen",
                                         "Chloramphenicol" = "#fdb462",  "Streptomycin" = "#bebada", "Tetracycline" = "#80b1d3", "Ciprofloxacin" = "#d9d9d9","Sulphathiozole" = "#fb8072"))+
  theme(title = element_text(size = 18), legend.text=element_text(size=18), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15))+
  # zoom in to last 3 bars on the axis
  scale_x_discrete(limits = tail(levels(reorder(Rallabx_data$ST, Rallabx_data$Count, FUN = sum)), 10)) + #looks at top 10
  coord_flip()+
  #geom_text(size = 3, position = position_stack(vjust = 0.5))
  theme(plot.title = element_text(hjust = 0.5)) #centre the title
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Rallabx_serovar.png', width = 16, height = 9, dpi = 100)




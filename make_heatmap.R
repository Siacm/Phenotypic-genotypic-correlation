library("ggplot2")
library("plyr")
library("reshape2")
library("scales")
library(tidyr)
library(dplyr)
library(openxlsx)

getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_antibiotic_counts070519.xlsx")
all_rsero_abx_count<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_antibiotic_counts070519.xlsx", sheet = "all_rsero_abx_count",colNames = T)
all_rsero_abx_count<-all_rsero_abx_count[,c(1:3)]

ggplot(all_rsero_abx_count, aes(x = Antibiotic, y = reorder(ST,Count), fill = Antibiotic, alpha = Count)) +
  # Dummy tile geom for white background
  geom_tile(fill = "white", alpha = 1, color="transparent") +
  geom_tile() +
  coord_fixed(ratio = 0.4)+
  labs(x="",y="",title = "Resistance to Antibiotic by Serovar") +
  theme(panel.background = element_rect(fill = 'white'))+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank())+
  scale_alpha_continuous(breaks = seq(0, max(all_rsero_abx_count$Count), length.out = 10), 
                         limits = c(-200, NA))+ #,guide = "none")+
  scale_fill_manual(values=c("Ampicillin"="#8dd3c7","Cefotaxime"="turquoise","Meropenem"="seagreen","Ciprofloxacin"="#d9d9d9",
                                "Sulphathiozole"="#fb8072","Tetracycline"="#80b1d3","Chloramphenicol"="#fdb462","Azithromycin"="#ffffb3","Streptomycin"="#bebada"),
                    breaks=c("Ampicillin","Cefotaxime","Meropenem","Azithromycin","Streptomycin",
                             "Chloramphenicol","Tetracycline", "Ciprofloxacin","Sulphathiozole"))+
  scale_x_discrete(limits=c("Ampicillin","Cefotaxime","Meropenem","Azithromycin","Streptomycin",
                            "Chloramphenicol","Tetracycline", "Ciprofloxacin","Sulphathiozole"))

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/serovar_abx_heatmap.png', width = 20, height = 20, dpi = 200)

  
install.packages("BiocManager", dependencies=TRUE)
install.packages("main_gp_tablea.table", dependencies=TRUE)
install.packages("openxlsx", dependencies=TRUE)
install.packages("ggrepel", dependencies=TRUE)
install.packages("forcats", dependencies=TRUE)

library(openxlsx)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(forcats)
library(scales)
library(survival)
library(openxlsx)
library(data.table) 

specsens<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_origrepeat_070519.xlsx", sheet = "specsens" )

#azithsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[1]*100, digits = 4)),"% sensitivity","\n")
#ampsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[2]*100, digits = 4)),"% sensitivity","\n")
#cefsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[3]*100, digits = 4)),"% sensitivity","\n")
#mersens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[4]*100, digits = 4)),"% sensitivity","\n")
#strepsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[5]*100, digits = 4)),"% sensitivity","\n")
#chlorsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[6]*100, digits = 4)),"% sensitivity","\n")
#tetsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[7]*100, digits = 4)),"% sensitivity","\n")
#cipsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[8]*100, digits = 4)),"% sensitivity","\n")
#sulsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[9]*100, digits = 4)),"% sensitivity","\n")

#azithspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[1]*100, digits = 4)),"% specificity","\n")
#ampspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[2]*100, digits = 4)),"% specificity","\n")
#cefspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[3]*100, digits = 4)),"% specificity","\n")
#merspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[4]*100, digits = 4)),"% specificity","\n")
#strepspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[5]*100, digits = 4)),"% specificity","\n")
#chlorspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[6]*100, digits = 4)),"% specificity","\n")
#tetspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[7]*100, digits = 4)),"% specificity","\n")
#cipspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[8]*100, digits = 4)),"% specificity","\n")
#sulspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[9]*100, digits = 4)),"% specificity","\n")

azithsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[1]*100, digits = 4)),"% sensitivity","\n")
ampsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[2]*100, digits = 4)),"% sensitivity","\n")
cefsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[3]*100, digits = 4)),"% sensitivity","\n")
mersens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[4]*100, digits = 4)),"% sensitivity","\n")
strepsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[5]*100, digits = 4)),"% sensitivity","\n")
chlorsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[6]*100, digits = 4)),"% sensitivity","\n")
tetsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[7]*100, digits = 4)),"% sensitivity","\n")
cipsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[8]*100, digits = 4)),"% sensitivity","\n")
sulsens<-paste("\n",sprintf("%2.2f", round(specsens$Sensitivity[9]*100, digits = 4)),"% sensitivity","\n")

azithspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[1]*100, digits = 4)),"% specificity","\n")
ampspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[2]*100, digits = 4)),"% specificity","\n")
cefspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[3]*100, digits = 4)),"% specificity","\n")
merspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[4]*100, digits = 4)),"% specificity","\n")
strepspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[5]*100, digits = 4)),"% specificity","\n")
chlorspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[6]*100, digits = 4)),"% specificity","\n")
tetspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[7]*100, digits = 4)),"% specificity","\n")
cipspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[8]*100, digits = 4)),"% specificity","\n")
sulspec<-paste("\n",sprintf("%2.2f", round(specsens$Specificity[9]*100, digits = 4)),"% specificity","\n")

####to make barplot of resistance geno vs pheno
getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_origrepeat_070519.xlsx")
main_gpresis_table<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_origrepeat_070519.xlsx", sheet = "resistance_numbers" )

main_gpresis_table.g <- gather(main_gpresis_table, type, value, -Antibiotic)  ###to arrange barplots according to pheno r, geno r, pheno s, geno s...

main_gpresis_table.g$group<-c(rep(azithsens, 1), rep(ampsens, 1),rep(cefsens, 1), rep(mersens, 1), rep(strepsens, 1),rep(chlorsens, 1), rep(tetsens, 1), rep(cipsens, 1), rep(sulsens, 1),
                              rep(azithsens, 1), rep(ampsens, 1),rep(cefsens, 1), rep(mersens, 1), rep(strepsens, 1),rep(chlorsens, 1), rep(tetsens, 1), rep(cipsens, 1), rep(sulsens, 1))
main_gpresis_table.g$class<-c(rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1),
                              rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1))
main_gpresis_table.g$type <- factor(main_gpresis_table.g$type,levels = c("Phenotypically.resistant", "Genotypically.resistant"))

main_gpresis_table.g$class <- factor(main_gpresis_table.g$class, levels = c("Beta-lactam","Macrolide", "Aminoglycoside",
                                                                            "Chloramphenicol", "Tetracycline","Fluoroquinolone","Sulfonamide"))

p<-ggplot(main_gpresis_table.g, aes(fill=type, x=Antibiotic, y=value)) + geom_bar(aes(fill = type), stat = "identity", position = "dodge", colour="white") + geom_text(aes(label=value), position=position_dodge(width=1.2), vjust=-0.5)+ scale_fill_manual(values=c("#de2d26","#fc9272")) + theme(title = element_text(size = 18), legend.text=element_text(size=18), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15)) + labs(y = "Number of isolates", fill = "", title = "Phenotypic versus Genotypic Resistance (Total isolates: 4,726)") +
  theme(plot.title = element_text(hjust = 0.5))
p+facet_grid(.~class+group, scales= "free_x", space= "free_x")+
  theme(strip.text.x = element_text(size = 10.5))

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/gpresistance_image.png', width = 16, height = 9, dpi = 100)

####to make barplot of sensitive geno vs pheno
main_gpsens_table<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_origrepeat_070519.xlsx", sheet = "sensitive_numbers" )

main_gpsens_table.g <- gather(main_gpsens_table, type, value, -Antibiotic)  ###to arrange barplots according to pheno r, geno r, pheno s, geno s...
main_gpsens_table.g$group<-c(rep(azithspec, 1), rep(ampspec, 1),rep(cefspec, 1), rep(merspec, 1), rep(strepspec, 1),rep(chlorspec, 1), rep(tetspec, 1), rep(cipspec, 1), rep(sulspec, 1),
                             rep(azithspec, 1), rep(ampspec, 1),rep(cefspec, 1), rep(merspec, 1), rep(strepspec, 1),rep(chlorspec, 1), rep(tetspec, 1), rep(cipspec, 1), rep(sulspec, 1))
main_gpsens_table.g$class<-c(rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1),
                              rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1))
main_gpsens_table.g$type <- factor(main_gpsens_table.g$type,levels = c("Phenotypically.susceptible", "Genotypically.susceptible"))

main_gpsens_table.g$class <- factor(main_gpsens_table.g$class, levels = c("Beta-lactam","Macrolide", "Aminoglycoside",
                                                                            "Chloramphenicol", "Tetracycline","Fluoroquinolone","Sulfonamide"))

p<-ggplot(main_gpsens_table.g, aes(fill=type, x=Antibiotic, y=value)) + geom_bar(aes(fill = type), stat = "identity", position = "dodge", colour="white") + geom_text(aes(label=value), position=position_dodge(width=1.2), vjust=-0.5)+ scale_fill_manual(values=c("#31a354", "#a1d99b")) + theme(title = element_text(size = 20), legend.text=element_text(size=20), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15)) + labs(y = "Number of isolates", fill = "", title = "Phenotypic versus Genotypic Susceptibilities (Total isolates: 4,726)") +
  theme(plot.title = element_text(hjust = 0.5))
p+facet_grid(.~class+group, scales= "free_x", space= "free_x")+
  theme(strip.text.x = element_text(size = 10))

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/gpsensitive_image.png', width = 16, height = 9, dpi = 100)

####to make barplot of all geno vs pheno
all_gp_table<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_origrepeat_070519.xlsx", sheet = "all_numbers" )

all_gp_table.g <- gather(all_gp_table, type, value, -Antibiotic)  ###to arrange barplots according to pheno r, geno r, pheno s, geno s...
#all_gp_table.g$group<-c(rep("Part1", 1), rep("Part2", 1),rep("Part3", 1), rep("Part4", 1), rep("Part5", 1),rep("Part6", 1), rep("Part7", 1), rep("Part8", 1), rep("Part9", 1),
#                             rep("Part1", 1), rep("Part2", 1),rep("Part3", 1), rep("Part4", 1), rep("Part5", 1),rep("Part6", 1), rep("Part7", 1), rep("Part8", 1), rep("Part9", 1),
#                             rep("Part1", 1), rep("Part2", 1),rep("Part3", 1), rep("Part4", 1), rep("Part5", 1),rep("Part6", 1), rep("Part7", 1), rep("Part8", 1), rep("Part9", 1),
#                             rep("Part1", 1), rep("Part2", 1),rep("Part3", 1), rep("Part4", 1), rep("Part5", 1),rep("Part6", 1), rep("Part7", 1), rep("Part8", 1), rep("Part9", 1))
all_gp_table.g$class<-c(rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1),
                             rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1),
                             rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1),
                             rep("Macrolide", 1), rep("Beta-lactam", 1),rep("Beta-lactam", 1), rep("Beta-lactam", 1), rep("Aminoglycoside", 1),rep("Chloramphenicol", 1), rep("Tetracycline", 1), rep("Fluoroquinolone", 1), rep("Sulfonamide", 1))
all_gp_table.g$type <- factor(all_gp_table.g$type,levels = c("Phenotypically.resistant", "Genotypically.resistant","Phenotypically.susceptible", "Genotypically.susceptible"))

all_gp_table.g$class <- factor(all_gp_table.g$class, levels = c("Beta-lactam","Macrolide", "Aminoglycoside",
                                                                            "Chloramphenicol", "Tetracycline","Fluoroquinolone","Sulfonamide"))


p<-ggplot(all_gp_table.g, aes(fill=type, x=Antibiotic, y=value)) + geom_bar(aes(fill = type), stat = "identity", position = "dodge", colour="white") + geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.5, size=2.5) + scale_fill_manual(values=c("#de2d26","#fc9272","#31a354", "#a1d99b")) + theme(title = element_text(size = 16), legend.text=element_text(size=16), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15)) + labs(y = "Number of isolates", fill = "", title = "All Phenotypic versus Genotypic Sensitivities (Total isolates: 4,726)") +
  theme(plot.title = element_text(hjust = 0.5))
p+facet_grid(.~class, scales= "free_x", space= "free_x")+
  theme(strip.text.x = element_text(size = 10.5))

ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/allgp_image.png', width = 16, height = 9, dpi = 100)

####to make barplot of overall geno vs pheno
main_gp_table<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_origrepeat_070519.xlsx", sheet = "overall_numbers" )

main_gp_table.g <- gather(main_gp_table, type, value, -Antibiotics)  ###to arrange barplots according to pheno r, geno r, pheno s, geno s...

main_gp_table.g$type <- factor(main_gp_table.g$type,levels = c("Phenotypically.resistant", "Genotypically.resistant","Phenotypically.susceptible", "Genotypically.susceptible"))
ggplot(main_gp_table.g, aes(fill=type, x=Antibiotics, y=value)) + geom_bar(aes(fill = type), stat = "identity", position = "dodge", colour="white") + geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.5, size=2.5) + scale_fill_manual(values=c("#de2d26","#fc9272","#31a354", "#a1d99b")) + theme(title = element_text(size = 20), legend.text=element_text(size=20), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15)) + labs(y = "Number of isolates", fill = "", title = "Overall Phenotypic versus Genotypic Sensitivities (Total isolates: 4,726)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/overall_image.png', width = 16, height = 9, dpi = 100)



install.packages("openxlsx", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("plotly", dependencies=TRUE)

library(openxlsx)
library(ggplot2)
library(plotly)

getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_comparison_070519_trainingset.xlsx")
sensitivity_plot<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_comparison_070519_trainingset.xlsx", sheet = "Total95_boxplot" )
colnames(sensitivity_plot)

# grouped boxplot
ggplot(sensitivity_plot, aes(x=Coverage, y=Specificity, color=Coverage)) + geom_point(size=5) +
  theme(legend.title=element_blank(),title = element_text(size = 20),legend.text=element_text(size=15), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15)) + 
  labs(y = "Specificity(%)", x = "Coverage(%)", fill = "", title = "Specificity vs Coverage at 95% Identity") + theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0.80, 1.3, by=0.1), limits=c(0.80,1.3))+scale_x_continuous(breaks = seq(30, 120, by=10))+
  geom_text(aes(label=paste("\n",sprintf("%0.2f", round(sensitivity_plot$Specificity*100, digits = 4)),"% specificity","\n",sensitivity_plot$Coverage,"% COV")), hjust=0.5, vjust=1, check_overlap = TRUE, size=4)
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/specificity_plot.png', width = 16, height = 9, dpi = 100)


# grouped boxplot
ggplot(sensitivity_plot, aes(x=Coverage, y=Sensitivity, color=Coverage)) + geom_point(size=5) +
  theme(legend.title=element_blank(),title = element_text(size = 20),legend.text=element_text(size=15), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15)) + 
  labs(y = "Sensitivity(%)", x = "Coverage(%)", fill = "", title = "Sensitivity vs Coverage at 95% Identity") + theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0.90, 1, by=0.01), limits=c(0.90,1))+ 
  geom_text(aes(label=paste("\n",sprintf("%0.2f", round(sensitivity_plot$Sensitivity*100, digits = 4)),"% sensitivity","\n",sensitivity_plot$Coverage,"% COV")), hjust=0.5, vjust=1, check_overlap = TRUE, size=4)
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/sensitivity_plot.png', width = 16, height = 9, dpi = 100)

## specificity and sensitivity of validation vs training set
training_plot<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_comparison_070519_trainingset.xlsx", sheet = "training" )
validation_plot<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/resistance_workbook_comparison_070519_validationset.xlsx", sheet = "validation" )
trainingvsvalidation <- rbind(data.frame(training_plot,Set="Training"), data.frame(validation_plot,Set="Validation"))
ggplot(trainingvsvalidation, aes(Sensitivity,Specificity,group=Set,col=Set)) + 
  geom_point(size=5) +scale_y_continuous(breaks = seq(0.90, 1, by=0.01), limits=c(0.90,1))

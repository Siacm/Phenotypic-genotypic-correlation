install.packages("BiocManager", dependencies=TRUE)
install.packages("main_gp_tablea.table", dependencies=TRUE)
install.packages("openxlsx", dependencies=TRUE)
install.packages("ggrepel", dependencies=TRUE)
install.packages("forcats", dependencies=TRUE)
install.packages("ggstance", dependencies=TRUE)

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
library(ggstance)

getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx")
all_rgenes<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/Resistance_gene_counts070519.xlsx", sheet = "all_abx_rgene_count" )

all_rgenes$Antibiotic <- factor(all_rgenes$Antibiotic, levels = c("Ampicillin","Cefotaxime", "Meropenem","Azithromycin",
                                                                            "Streptomycin", "Chloramphenicol","Tetracycline","Ciprofloxacin","Sulphathiozole"))
####arranges highest barplot in order of each facet
all_rgenes$n = as.numeric(factor(all_rgenes$Antibiotic))
all_rgenes = ddply(all_rgenes,.(Antibiotic,Resistance.gene),transform, x=paste(c(rep(' ',n-1), Resistance.gene), collapse=''))
all_rgenes$x = factor(all_rgenes$x, levels=all_rgenes[order(all_rgenes$Count), 'x'])

####creates plot, if you want to re-order,make y=-Count
p<-ggplot(all_rgenes, aes(x=x, y=Count, fill=Antibiotic)) + geom_bar(stat = "identity") + #scale_fill_manual(values=c("#de2d26","#fc9272"))+
  theme(title = element_text(size = 15), axis.text.x=element_text(size=15), axis.text=element_text(size=15),legend.position = "none",panel.background = element_rect(fill = "whitesmoke")) + 
  labs(x = NULL,y = "Number of isolates", fill = NULL, title = "") +
  theme(plot.title = element_text(hjust = 0.5))+coord_flip()
p

####adds color for bars in each facet
p<-p+facet_grid(Antibiotic~., scales= "free_y", space= "free_y")+ theme(strip.text.y = element_text(angle = 0, size=15))+
  scale_fill_manual(values= c("#8dd3c7","turquoise","seagreen","#ffffb3","#bebada","#fdb462","#80b1d3","#d9d9d9","#fb8072"))+
  scale_y_continuous(breaks=seq(0, 600, 50),expand=c(0,10))
p

###adds colour to the actual facet
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#8dd3c7","turquoise","seagreen","#ffffb3","#bebada","#fdb462","#80b1d3","#d9d9d9","#fb8072")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

png("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/allrgenes_image.png",width = 1080, height = 1080); plot(g); dev.off()


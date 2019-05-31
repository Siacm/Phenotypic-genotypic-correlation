library(png)
library(igraph)
install.packages("rasterVis")
library(raster)
library(rasterVis)

# To get an image to plot
img_train <- list.files(path = "~/Documents/PhD/Chapter_2/Chapter_2_writeup/co-occurrence pics", pattern = ".png", full.names=TRUE)
img_list <- lapply(img_train, readPNG)

g <- graph.ring(4)
# This is a complex attribute, so supply a list here

V(g)$raster <-img_list

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/all_rasterimages.png", 1000,1000)
plot(g, vertex.shape="raster", vertex.label=NA, vertex.size=160, vertex.size2=160, edge.lty=0)
dev.off()
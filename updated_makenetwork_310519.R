#Datasheet of alleles as column names and no samples names is
#miniCooccur_inputNoSampleNames
install.packages("igraph") 

library(reshape2)
library(ggplot2)
library(igraph)

miniCooccur <- read.csv("~/Documents/PhD/DemoCooccurrence/matched_merged_MDR_forRscript.csv",header=TRUE)
colnames(miniCooccur)<-c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB","blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","catA1","catA2","catB3","cmlA","cmlA1","cmlA5","ere(A)","erm(A)","erm(B)","floR","mef(B)","mph(A)","msr(E)","oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul2","sul3","tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A")

miniCooccur <- as.matrix(miniCooccur)
miniCooccur
miniCooccur_out <- crossprod(miniCooccur)
diag(miniCooccur_out) <- 0
write.csv(x = miniCooccur_out, file="ExampleAmpCooccurrence_dataset.csv")

miniCooccur_melt <- melt(miniCooccur_out)


#PLOTTING HEATMAP
ggplot(data=miniCooccur_melt, aes(x=Var1,y=Var2,fill=value)) + geom_tile()

ggplot(data = miniCooccur_melt, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "red", mid = "white")

ggplot(data = miniCooccur_melt, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "red", mid = "white") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#To plot out the networks with igraph

library(igraph)

mini <- graph.adjacency(miniCooccur_out, weighted=TRUE, mode ='undirected')
mini <- simplify(mini)
# set labels and degrees of vertices
V(mini)$label <- V(mini)$name
V(mini)$degree <- degree(mini)
plot(mini)

plot.igraph(mini,vertex.label=V(mini)$name, edge.color="black",edge.width=E(mini)$weight)
mini.e<-igraph::as_data_frame(mini, what="edges")
mini.v<-igraph::as_data_frame(mini, what="vertices")

hist(mini.e$weight)
mean<-mean(mini.e$weight)
#this calculates the mean of the weights
sd<-sd(mini.e$weight)
#this calculates the sd of the weights


###CUTOFF=NINESD
NINESD<-mean+9*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- NINESD
miniNINESD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniNINESD) 

#If you want to delete all no unconnected nodes
iso <- V(miniNINESD)[degree(miniNINESD)==0]
miniNINESD <- delete.vertices(miniNINESD, iso)

#To add edge weight labels
cut.off <- NINESD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())


plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minNINESD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minNINESD", edge.label = E(g2)$weight, layout=coords)


V(miniNINESD)$color <- ifelse(V(miniNINESD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                               ifelse(V(miniNINESD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                      ifelse(V(miniNINESD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                             ifelse(V(miniNINESD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                    ifelse(V(miniNINESD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                           ifelse(V(miniNINESD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                  ifelse(V(miniNINESD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniNINESD, vertex.label=V(miniNINESD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥427 connections",cex.main=2)
png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/9SD.png", 700,700)

plot.igraph(miniNINESD, vertex.label=V(miniNINESD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥427 connections",cex.main=2)

dev.off()



###CUTOFF=EIGHTSD
EIGHTSD<-mean+8*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- EIGHTSD
miniEIGHTSD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniEIGHTSD) 

#If you want to delete all no unconnected nodes
iso <- V(miniEIGHTSD)[degree(miniEIGHTSD)==0]
miniEIGHTSD <- delete.vertices(miniEIGHTSD, iso)

#To add edge weight labels
cut.off <- EIGHTSD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())


plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minEIGHTSD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minEIGHTSD", edge.label = E(g2)$weight, layout=coords)


V(miniEIGHTSD)$color <- ifelse(V(miniEIGHTSD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                               ifelse(V(miniEIGHTSD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                      ifelse(V(miniEIGHTSD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                             ifelse(V(miniEIGHTSD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                    ifelse(V(miniEIGHTSD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                           ifelse(V(miniEIGHTSD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                  ifelse(V(miniEIGHTSD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniEIGHTSD, vertex.label=V(miniEIGHTSD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50,  edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥381 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/8SD.png", 700,700)

plot.igraph(miniEIGHTSD, vertex.label=V(miniEIGHTSD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥381 connections",cex.main=2)

dev.off()

###CUTOFF=SEVENSD
SEVENSD<-mean+7*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- SEVENSD
miniSEVENSD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniSEVENSD) 

#If you want to delete all no unconnected nodes
iso <- V(miniSEVENSD)[degree(miniSEVENSD)==0]
miniSEVENSD <- delete.vertices(miniSEVENSD, iso)

#To add edge weight labels
cut.off <- SEVENSD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())


plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minSEVENSD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minSEVENSD", edge.label = E(g2)$weight, layout=coords)


V(miniSEVENSD)$color <- ifelse(V(miniSEVENSD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                             ifelse(V(miniSEVENSD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                    ifelse(V(miniSEVENSD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                           ifelse(V(miniSEVENSD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                  ifelse(V(miniSEVENSD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                         ifelse(V(miniSEVENSD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                ifelse(V(miniSEVENSD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniSEVENSD, vertex.label=V(miniSEVENSD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥336 connections",cex.main=2)
png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/7SD.png", 700,700)

plot.igraph(miniSEVENSD, vertex.label=V(miniSEVENSD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥336 connections",cex.main=2)

dev.off()

###CUTOFF=SIXSD
SIXSD<-mean+6*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- SIXSD
miniSIXSD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniSIXSD) 

#If you want to delete all no unconnected nodes
iso <- V(miniSIXSD)[degree(miniSIXSD)==0]
miniSIXSD <- delete.vertices(miniSIXSD, iso)

#To add edge weight labels
cut.off <- SIXSD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())


plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minSIXSD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minSIXSD", edge.label = E(g2)$weight, layout=coords)


V(miniSIXSD)$color <- ifelse(V(miniSIXSD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                              ifelse(V(miniSIXSD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                     ifelse(V(miniSIXSD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                            ifelse(V(miniSIXSD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                   ifelse(V(miniSIXSD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                          ifelse(V(miniSIXSD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                 ifelse(V(miniSIXSD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniSIXSD, vertex.label=V(miniSIXSD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥290 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/6SD.png", 700,700)

plot.igraph(miniSIXSD, vertex.label=V(miniSIXSD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥290 connections",cex.main=2)

dev.off()


###CUTOFF=FIVESD
FIVESD<-mean+5*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- FIVESD
miniFIVESD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniFIVESD) 

#If you want to delete all no unconnected nodes
iso <- V(miniFIVESD)[degree(miniFIVESD)==0]
miniFIVESD <- delete.vertices(miniFIVESD, iso)

#To add edge weight labels
cut.off <- FIVESD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())


plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minFIVESD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minFIVESD", edge.label = E(g2)$weight, layout=coords)


V(miniFIVESD)$color <- ifelse(V(miniFIVESD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                              ifelse(V(miniFIVESD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                     ifelse(V(miniFIVESD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                            ifelse(V(miniFIVESD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                   ifelse(V(miniFIVESD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                          ifelse(V(miniFIVESD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                 ifelse(V(miniFIVESD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniFIVESD, vertex.label=V(miniFIVESD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥244 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/5SD.png", 700,700)

plot.igraph(miniFIVESD, vertex.label=V(miniFIVESD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥244 connections",cex.main=2)


dev.off()
###CUTOFF=FOURSD
FOURSD<-mean+4*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- FOURSD
miniFOURSD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniFOURSD) 

#If you want to delete all no unconnected nodes
iso <- V(miniFOURSD)[degree(miniFOURSD)==0]
miniFOURSD <- delete.vertices(miniFOURSD, iso)

#To add edge weight labels
cut.off <- FOURSD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())


plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minFOURSD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minFOURSD", edge.label = E(g2)$weight, layout=coords)


V(miniFOURSD)$color <- ifelse(V(miniFOURSD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                               ifelse(V(miniFOURSD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                      ifelse(V(miniFOURSD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                             ifelse(V(miniFOURSD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                    ifelse(V(miniFOURSD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                           ifelse(V(miniFOURSD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                  ifelse(V(miniFOURSD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniFOURSD, vertex.label=V(miniFOURSD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥199 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/4SD.png", 700,700)

plot.igraph(miniFOURSD, vertex.label=V(miniFOURSD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥199 connections",cex.main=2)

dev.off()

###CUTOFF=THREESD
THREESD<-mean+3*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- THREESD
miniTHREESD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniTHREESD) 

#If you want to delete all no unconnected nodes
iso <- V(miniTHREESD)[degree(miniTHREESD)==0]
miniTHREESD <- delete.vertices(miniTHREESD, iso)

#To add edge weight labels
cut.off <- THREESD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minTHREESD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minTHREESD", edge.label = E(g2)$weight, layout=coords)


V(miniTHREESD)$color <- ifelse(V(miniTHREESD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                             ifelse(V(miniTHREESD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                    ifelse(V(miniTHREESD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                           ifelse(V(miniTHREESD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                  ifelse(V(miniTHREESD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                         ifelse(V(miniTHREESD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                ifelse(V(miniTHREESD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))

#Grey connections with edge labels
plot.igraph(miniTHREESD, vertex.label=V(miniTHREESD)$name, vertex.size=24,edge.color="grey",edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥153 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/3SD.png", 700,700)

plot.igraph(miniTHREESD, vertex.label=V(miniTHREESD)$name, vertex.size=24,edge.color="black",edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥153 connections",cex.main=2)

dev.off()
###CUTOFF=TWOSD
TWOSD<-mean+2*(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- TWOSD
miniTWOSD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniTWOSD) 

#If you want to delete all no unconnected nodes
iso <- V(miniTWOSD)[degree(miniTWOSD)==0]
miniTWOSD <- delete.vertices(miniTWOSD, iso)

#To add edge weight labels
cut.off <- TWOSD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minTWOSD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minTWOSD", edge.label = E(g2)$weight, layout=coords)


V(miniTWOSD)$color <- ifelse(V(miniTWOSD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                             ifelse(V(miniTWOSD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                    ifelse(V(miniTWOSD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                           ifelse(V(miniTWOSD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                  ifelse(V(miniTWOSD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                         ifelse(V(miniTWOSD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                ifelse(V(miniTWOSD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))
#Grey connections with edge labels
plot.igraph(miniTWOSD, vertex.label=V(miniTWOSD)$name, edge.color="grey",vertex.size=24,edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥108 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/2SD.png", 700,700)

#no grey connection labels
plot.igraph(miniTWOSD, vertex.label=V(miniTWOSD)$name, edge.color="black",vertex.size=24,edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥108 connections",cex.main=2)

dev.off()

###CUTOFF=ONESD
ONESD<-mean+(sd)

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- ONESD
miniONESD <- delete_edges(mini, E(mini)[weight<cut.off])
plot(miniONESD) 

#If you want to delete all no unconnected nodes
iso <- V(miniONESD)[degree(miniONESD)==0]
miniONESD <- delete.vertices(miniONESD, iso)


#To add edge weight labels
cut.off <- ONESD
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minONESD", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minONESD", edge.label = E(g2)$weight, layout=coords)


V(miniONESD)$color <- ifelse(V(miniONESD)$name %in% c("aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","strA","aph(6)-Ic","strB"),"#bebada",
                             ifelse(V(miniONESD)$name %in% c("blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40"),"#8dd3c7",
                                    ifelse(V(miniONESD)$name %in% c("catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR"),"#fdb462",
                                           ifelse(V(miniONESD)$name %in% c("ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)"),"#ffffb3",
                                                  ifelse(V(miniONESD)$name %in% c("oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A"),"#d9d9d9",
                                                         ifelse(V(miniONESD)$name %in% c("sul1","sul2","sul3"),"#fb8072",
                                                                ifelse(V(miniONESD)$name %in% c("tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"),"#80b1d3","white")))))))


#Grey connections with edge labels
plot.igraph(miniONESD, vertex.label=V(miniONESD)$name, edge.color="grey",vertex.size=24,edge.width=E(net.sp)$weight/50,edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥62 connections",cex.main=2)

png ("~/Documents/PhD/Chapter_2/Chapter_2_writeup/1SD.png", 700,700)

#no grey connection labels
plot.igraph(miniONESD, vertex.label=V(miniONESD)$name, edge.color="black",vertex.size=24,edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurrence of MDR Isolates, ≥62 connections",cex.main=2)

dev.off()

###CUTOFF=mean

#This sets the cutoff to the mean value and then gets it to delete any edges with less weight than the cutoff specified
cut.off <- mean
minimean <- delete_edges(mini, E(mini)[weight<cut.off])
plot(minimean) 

#If you want to delete all no unconnected nodes
iso <- V(minimean)[degree(minimean)==0]
minimean <- delete.vertices(minimean, iso)


#To add edge weight labels
cut.off <- mean
net.sp <- delete_edges(mini, E(mini)[weight<cut.off])
iso <- V(net.sp)[degree(net.sp)==0]
g2 <- delete.vertices(net.sp, iso)
coords <- layout_(g2, as_star())
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minmean", edge.label = E(g2)$weight)
plot.igraph(g2, vertex.label=V(g2)$name, edge.color="grey", alpha = 0.5,edge.width=E(g2)$weight/50, main="minmean", edge.label = E(g2)$weight, layout=coords)

#Grey connections with edge labels
plot.igraph(minimean, vertex.label=V(minimean)$name, edge.color="grey",vertex.size=24,edge.width=E(net.sp)$weight/50, edge.label = E(g2)$weight, layout=layout_nicely)
title("Co-occurence of MDR Isolates, ≥16 connections",cex.main=2)

#no grey connection labels
plot.igraph(minimean, vertex.label=V(minimean)$name, edge.color="black",vertex.size=24,edge.width=E(net.sp)$weight/50, layout=layout_nicely)
title("Co-occurence of MDR Isolates, ≥16 connections",cex.main=2)




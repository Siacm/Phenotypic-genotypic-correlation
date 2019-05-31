install.packages("BiocManager", dependencies=TRUE)
install.packages("data.table", dependencies=TRUE)
install.packages("dplyr", dependencies=TRUE)
library("readr")
library("data.table")
library("dplyr")

##reads phenotypic and genotypic files ********sistr has double quotes, remove first (re-save as csv on libreoffice to remove quotes******
genotypic <- read.csv("~/Documents/PhD/Chapter_2/Run_results/abricate_results/95id60cov_mdusamplesinputbig_run070319_metadataAll_simplified.csv")
names(genotypic) <- gsub(x = names(genotypic), pattern = "0_Abricate_ncbi_", replacement = "")
names(genotypic) <- gsub(x = names(genotypic), pattern = "Abricate_ncbi_", replacement = "")
colnames(genotypic) ##to get index numbers below 77:222 to change all genes that had "-" to "." back to "-"
names(genotypic)[77:222]<- c ("A7J11_01233","A7J11_02581","A7J11_04295","A7J11_04464","aac(3)-II","aac(3)-IIa","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","aad9","aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","ant(2'')-Ia","aph(3'')-Ib","aph(3')-IIa","aph(3')-Ia","aph(3')-VI","aph(4)-Ia","aph(6)-Ic","aph(6)-Id","arr","arr-2","arr-3","blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","ble","bleO","catA1","catA2","catB3","cmlA","cmlA1","cmlA5","dfrA1","dfrA10","dfrA12","dfrA14","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","ere(A)","erm(42)","erm(A)","erm(B)","floR","fosA3","fosA4","fosA7","lnu(C)","lnu(F)","lnu(G)","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","mef(B)","mph(A)","mph(E)","msr(E)","oqxA","oqxA2","oqxB","qacE","qacEdelta1","qacG2","qacL","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul2","sul3","tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)")
colnames(genotypic)[1]<- c ("ID")
#genotypic <- read_csv("~/Documents/PhD/Chapter_2/amr_files/resistance_data/Original_files/all_mdu_enterica_amr.csv", trim_ws = TRUE)
phenotypic<- read.csv("~/Documents/PhD/Chapter_2/amr_files/resistance_data/updated_phenotypic_amr.csv")

#convert the blanks to 0 (except column1) to to be able to differentiate from missing geno/pheno
genotypic[-1][sapply(genotypic[-1], is.na)]<-0

##use fread for reading sistr output because it has "", it cant merge.
sistr_genotypic<- fread("~/Documents/PhD/Chapter_2/Sistr_output/All_S_enterica_isolates_MLST_SISTR.csv", quote ="\"")
sistr_genotypic_extra<- fread("~/Documents/PhD/Chapter_2/Sistr_output/Extra_S_enterica_isolates_MLST_SISTR.csv", quote ="\"")
sistr_genotypic_missing<-read.csv("~/Documents/PhD/Chapter_2/Sistr_output/missing_sistr.out.csv", stringsAsFactors = FALSE)
sistr_genotypic_missing_mlst<-read.csv("~/Documents/PhD/Chapter_2/Sistr_output/missing_mlst.csv", stringsAsFactors = FALSE)
ariba_gyrpar_known<-read.csv("~/Documents/PhD/Chapter_2/Run_results/ariba200219_90id100cov/200219_summary_known_var.csv", stringsAsFactors = FALSE)
colnames(ariba_gyrpar_known)[1]<- c ("ID")
ariba_gyrpar_known$ID<-gsub("/report.tsv","",ariba_gyrpar_known$ID)
ariba_gyrpar_known<-subset(ariba_gyrpar_known, select =-c(gyrA.match, gyrB.match, parC.match, parE.match))
ariba_gyrpar_novel<-read.csv("~/Documents/PhD/Chapter_2/Run_results/ariba200219_90id100cov/200219_summary_novel_var.csv", stringsAsFactors = FALSE) #novel ariba variants as per ncbi
colnames(ariba_gyrpar_novel)[1]<- c ("ID")
ariba_gyrpar_novel$ID<-gsub("/report.tsv","",ariba_gyrpar_novel$ID)
#ariba_gyrpar_novel<-subset(ariba_gyrpar_novel, select =c(ID, gyrB.S464F, gyrB.Y421F)) #added two novel variants -to be confirmed -used at 90% id 50% cov
ariba_gyrpar_novel<-subset(ariba_gyrpar_novel, select =c(ID, gyrB.S464F)) #added a novel variants -to be confirmed used at 90%id 60%cov and above
merged_ariba_gyrpar<-merge(ariba_gyrpar_known,ariba_gyrpar_novel, by= "ID", all=TRUE)

#merges all the missing sistr mlst and serovar 
sistr_missing<-merge(sistr_genotypic_missing, sistr_genotypic_missing_mlst, by ="ID", all=TRUE, stringsAsFactors = FALSE)
merged_sistr<-merge(sistr_genotypic[,c(1,11,15)],sistr_genotypic_extra[,c(1,11,15)], by =c("ID","serovar","m_MLST"), all=TRUE, stringsAsFactors = FALSE)
merged_sistr<-merge.data.frame(merged_sistr,sistr_missing[,c(1,6,10)], by =c("ID","serovar","m_MLST"), all=TRUE)

#get QC data for maybe results
genotypic_coverage<-read.csv("~/Documents/PhD/Chapter_2/Run_results/abricate_results/95id75cov_mdusamplesinput010319_metadataAll.csv")
names(genotypic_coverage) <- gsub(x = names(genotypic_coverage), pattern = "0_Abricate_ncbi_", replacement = "")
names(genotypic_coverage) <- gsub(x = names(genotypic_coverage), pattern = "Abricate_ncbi_", replacement = "")
colnames(genotypic_coverage)[24:232]<- c ("A7J11_01233","A7J11_02581","A7J11_02581.1","A7J11_02581.2","A7J11_04295","A7J11_04295.1","A7J11_04295.2","A7J11_04464","A7J11_04464.1","A7J11_04464.2","A7J11_04464.3","A7J11_04464.4","aac(3)-II","aac(3)-IIa","aac(3)-IIa.1","aac(3)-IIa.2","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","aad9","aadA1","aadA1.1","aadA1.2","aadA1.3","aadA1.4","aadA11","aadA12","aadA12.1","aadA15","aadA15.1","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","ant(2'')-Ia","aph(3'')-Ib","aph(3'')-Ib.1","aph(3'')-Ib.2","aph(3')-IIa","aph(3')-Ia","aph(3')-Ia.1","aph(3')-Ia.2","aph(3')-VI","aph(4)-Ia","aph(6)-Ic","aph(6)-Id","aph(6)-Id.1","aph(6)-Id.2","aph(6)-Id.3","aph(6)-Id.4","arr","arr-2","arr-3","blaBIL-1","blaBIL-1.1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-10.1","blaCTX-M-10.2","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-1.1","blaTEM-1.2","blaTEM-1.3","blaTEM-1.4","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-135.1","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","ble","bleO","catA1","catA2","catA2.1","catA2.2","catB3","catB3.1","catB3.2","cmlA","cmlA1","cmlA5","dfrA1","dfrA1.1","dfrA10","dfrA12","dfrA14","dfrA14.1","dfrA14.2","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","ere(A)","ere(A).1","erm(42)","erm(42).1","erm(A)","erm(B)","floR","floR.1","fosA3","fosA4","fosA7","fosA7.1","lnu(C)","lnu(F)","lnu(G)","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","mef(B)","mef(B).1","mef(B).2","mph(A)","mph(E)","msr(E)","oqxA","oqxA2","oqxB","qacE","qacE.1","qacEdelta1","qacEdelta1.1","qacG2","qacL","qacL.1","qacL.2","qacL.3","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul1.1","sul1.2","sul1.3","sul1.4","sul2","sul2.1","sul2.2","sul2.3","sul2.4","sul3","tet(A)","tet(A).1","tet(A).2","tet(A).3","tet(A).4","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(M).1","tet(X)")
genotypic_coverage<-genotypic_coverage[,c(1,24:232)]
colnames(genotypic_coverage)[1]<- c ("ID")
genotypic_qc<-fread("~/Documents/PhD/Chapter_2/Run_results/fa_qcrun_161218/qc_mdusamples.txt", quote ="\"")
colnames(genotypic_qc)[1]<- c ("ID")
genotypic_qc<-subset(genotypic_qc, select =c(ID, no, bp, N50)) #only these columns
genotypic_qc_cover<-merge(genotypic_qc,genotypic_coverage,by ="ID",all=TRUE)
genotypic_qc_cover<-merge(genotypic_qc_cover,genotypic, by =c("ID","A7J11_01233","A7J11_02581","A7J11_04295","A7J11_04464","aac(3)-II","aac(3)-IIa","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","aad9","aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","ant(2'')-Ia","aph(3'')-Ib","aph(3')-IIa","aph(3')-Ia","aph(3')-VI","aph(4)-Ia","aph(6)-Ic","aph(6)-Id","arr","arr-2","arr-3","blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","ble","bleO","catA1","catA2","catB3","cmlA","cmlA1","cmlA5","dfrA1","dfrA10","dfrA12","dfrA14","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","ere(A)","erm(42)","erm(A)","erm(B)","floR","fosA3","fosA4","fosA7","lnu(C)","lnu(F)","lnu(G)","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","mef(B)","mph(A)","mph(E)","msr(E)","oqxA","oqxA2","oqxB","qacE","qacEdelta1","qacG2","qacL","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul2","sul3","tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)"), all=TRUE)
genotypic_qc_cover<-merge(genotypic_qc_cover,merged_sistr, by = "ID", all=TRUE)
#genotypic_qc_cover$new <- apply(genotypic_qc_cover[,1:319], 1, function(x) any(x %in% c("maybe")))
#genotypic_qc_cover <- genotypic_qc_cover[(genotypic_qc_cover$new==TRUE),]
#genotypic_qc_cover$new <- NULL
genotypic_qc_cover$all_maybe <-  apply(genotypic_qc_cover, 1, function(x) paste("[", paste(colnames(genotypic_qc_cover)[grepl("maybe", x, ignore.case = T)], collapse = ","), "]", sep = "")) #makes a list with all genes that say maybe
genotypic_qc_cover<-genotypic_qc_cover[,c("ID", "serovar", "m_MLST", "all_maybe", "no", "bp","N50","A7J11_01233","A7J11_02581","A7J11_02581.1","A7J11_02581.2","A7J11_04295","A7J11_04295.1","A7J11_04295.2","A7J11_04464","A7J11_04464.1","A7J11_04464.2","A7J11_04464.3","A7J11_04464.4","aac(3)-II","aac(3)-IIa","aac(3)-IIa.1","aac(3)-IIa.2","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","aad9","aadA1","aadA1.1","aadA1.2","aadA1.3","aadA1.4","aadA11","aadA12","aadA12.1","aadA15","aadA15.1","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","ant(2'')-Ia","aph(3'')-Ib","aph(3'')-Ib.1","aph(3'')-Ib.2","aph(3')-IIa","aph(3')-Ia","aph(3')-Ia.1","aph(3')-Ia.2","aph(3')-VI","aph(4)-Ia","aph(6)-Ic","aph(6)-Id","aph(6)-Id.1","aph(6)-Id.2","aph(6)-Id.3","aph(6)-Id.4","arr","arr-2","arr-3","blaBIL-1","blaBIL-1.1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-10.1","blaCTX-M-10.2","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-1.1","blaTEM-1.2","blaTEM-1.3","blaTEM-1.4","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-135.1","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","ble","bleO","catA1","catA2","catA2.1","catA2.2","catB3","catB3.1","catB3.2","cmlA","cmlA1","cmlA5","dfrA1","dfrA1.1","dfrA10","dfrA12","dfrA14","dfrA14.1","dfrA14.2","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","ere(A)","ere(A).1","erm(42)","erm(42).1","erm(A)","erm(B)","floR","floR.1","fosA3","fosA4","fosA7","fosA7.1","lnu(C)","lnu(F)","lnu(G)","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","mef(B)","mef(B).1","mef(B).2","mph(A)","mph(E)","msr(E)","oqxA","oqxA2","oqxB","qacE","qacE.1","qacEdelta1","qacEdelta1.1","qacG2","qacL","qacL.1","qacL.2","qacL.3","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul1.1","sul1.2","sul1.3","sul1.4","sul2","sul2.1","sul2.2","sul2.3","sul2.4","sul3","tet(A)","tet(A).1","tet(A).2","tet(A).3","tet(A).4","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(M).1","tet(X)")]
write.csv(genotypic_qc_cover, file="~/Documents/PhD/Chapter_2/amr_files/r_output/genotypic_qccover_maybe.csv",row.names=FALSE, quote=TRUE)


#merges all the geno/pheno data (check to make that mdu_typed column says YES..remove all those that say no)
all_merged_gp_files<-merge(genotypic,merged_ariba_gyrpar, by ="ID", all=TRUE) 
all_merged_gp_files<-merge(genotypic[,-c(2)],phenotypic, by ="ID", all=TRUE)#removes all_genes column to create a new all_genes column that will have gyra/parc genes included (see below)
#all_merged_gp_files$all_genes <-  apply(all_merged_gp_files, 1, function(x) paste("[", paste(colnames(all_merged_gp_files)[x == "yes"], collapse = ","), "]", sep = ""))
all_merged_gp_files$all_genes <-  apply(all_merged_gp_files, 1, function(x) paste("[", paste(colnames(all_merged_gp_files)[grepl("yes", x, ignore.case = T)], collapse = ","), "]", sep = ""))
all_merged_gp_files<-merge(merged_sistr, all_merged_gp_files, by ="ID", all=TRUE)
matched_merged_gp_files<-merge(genotypic,phenotypic, by.x = "ID", all.x=TRUE)
matched_merged_gp_files<-merge(matched_merged_gp_files[,-c(2)], merged_ariba_gyrpar, by ="ID", all=TRUE)#removes all_genes column to create a new all_genes column that will have gyra/parc genes included (see below)
matched_merged_gp_files$all_genes <-  apply(matched_merged_gp_files, 1, function(x) paste("[", paste(colnames(matched_merged_gp_files)[grepl("yes", x, ignore.case = T)], collapse = ","), "]", sep = ""))
matched_merged_gp_files<-merge(matched_merged_gp_files, merged_sistr, by ="ID", all=FALSE)

#print 
write.csv(all_merged_gp_files, file="~/Documents/PhD/Chapter_2/amr_files/r_output/all_merged_gp_files.csv",row.names=FALSE, quote=TRUE)
matched_merged_gp_files<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","m_MLST","all_genes","genes_confirmed","genes_unconfirmed","A7J11_01233","A7J11_02581","A7J11_04295","A7J11_04464","aac(3)-II","aac(3)-IIa","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","aad9","aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","ant(2'')-Ia","aph(3'')-Ib","aph(3')-IIa","aph(3')-Ia","aph(3')-VI","aph(4)-Ia","aph(6)-Ic","aph(6)-Id","arr","arr-2","arr-3","blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","ble","bleO","catA1","catA2","catB3","cmlA","cmlA1","cmlA5","dfrA1","dfrA10","dfrA12","dfrA14","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","ere(A)","erm(42)","erm(A)","erm(B)","floR","fosA3","fosA4","fosA7","lnu(C)","lnu(F)","lnu(G)","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","mef(B)","mph(A)","mph(E)","msr(E)","oqxA","oqxA2","oqxB","qacE","qacEdelta1","qacG2","qacL","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul2","sul3","tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A","phage_type","ampicillin","streptomycin","tetracycline","chloramphenicol","sulphathiozole","trimethoprim","cotrimoxazole","kanamycin","nalidixic_acid","spectinomycin","gentamicin","ciprofloxacin","azithromycin","cefotaxime","meropenem","mdu_typed","card_type","site","collection_date")]
#matched_merged_gp_files<-matched_merged_gp_files[,c("ID", "MLST_Scheme","ST","m_MLST","all_genes","genes_confirmed","genes_unconfirmed","ampicillin","streptomycin","tetracycline","chloramphenicol","sulphathiozole","trimethoprim","cotrimoxazole","kanamycin","nalidixic_acid","spectinomycin","gentamicin","ciprofloxacin","azithromycin","cefotaxime","meropenem","A7J11_01233","A7J11_02581","A7J11_04295","A7J11_04464","aac(3)-II","aac(3)-IIa","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","aad9","aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","ant(2'')-Ia","aph(3'')-Ib","aph(3')-IIa","aph(3')-Ia","aph(3')-VI","aph(4)-Ia","aph(6)-Ic","aph(6)-Id","arr","arr-2","arr-3","blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","ble","bleO","catA1","catA2","catB3","cmlA","cmlA1","cmlA5","dfrA1","dfrA10","dfrA12","dfrA14","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","ere(A)","erm(42)","erm(A)","erm(B)","floR","fosA3","fosA4","fosA7","lnu(C)","lnu(F)","lnu(G)","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","mef(B)","mph(A)","mph(E)","msr(E)","oqxA","oqxA2","oqxB","qacE","qacEdelta1","qacG2","qacL","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","sul1","sul2","sul3","tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrA.Y100H","gyrB.E466D","parC_1.E84G","parC_1.S57T","parC_1.S80I","parC_1.S80R","parE_1.S458A","site","collection_date")]
colnames(matched_merged_gp_files)[4]<- c ("ST")
write.csv(matched_merged_gp_files, file="~/Documents/PhD/Chapter_2/amr_files/r_output/matched_merged_gp_files.csv", row.names=FALSE, quote=TRUE)

#merge to only training set numbers
training_set_numbers<- read.csv("~/Documents/PhD/Chapter_2/amr_files/resistance_data/Training_ids.csv")
matched_merged_gp_files<- merge(training_set_numbers,matched_merged_gp_files, by ="ID",all=FALSE)

#get the index number of each antibiotic and resistance gene
colnames(matched_merged_gp_files)

#sort samples resistance genes columns acdg to abx
gentamicin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","gentamicin","aac(3)-II","aac(3)-IIa","aac(3)-IId","aac(3)-IIe","aac(3)-IVa","aac(3)-Id","aac(3)-VIa","aac(6')-IIc","aac(6')-Il","ant(2'')-Ia","site","collection_date")]
gentamicin_resistance$gentamicin<-gsub("M|P","I", gentamicin_resistance$gentamicin)
kanamycin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","kanamycin","aac(6')-IIc","aac(6')-Ib-cr","aac(6')-If","aac(6')-Il","ant(2'')-Ia","aph(3')-IIa","aph(3')-Ia","aph(3')-VI","site","collection_date")]
kanamycin_resistance$kanamycin<-gsub("M|P","I", kanamycin_resistance$kanamycin)
streptomycin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","streptomycin","aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","aph(3'')-Ib","aph(6)-Ic","aph(6)-Id","site","collection_date")]
streptomycin_resistance$streptomycin<-gsub("M|P","I", streptomycin_resistance$streptomycin)
spectinomycin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","spectinomycin","aad9","aadA1","aadA11","aadA12","aadA15","aadA16","aadA2","aadA22","aadA4","aadA5","aadA7","site","collection_date")]
spectinomycin_resistance$spectinomycin<-gsub("M|P","I", spectinomycin_resistance$spectinomycin)
ampicillin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","ampicillin","blaBIL-1","blaCARB-2","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4","blaLAP-2","blaNDM-1","blaOXA-1","blaOXA-10","blaOXA-9","blaSHV-12","blaTEM-1","blaTEM-106","blaTEM-12","blaTEM-135","blaTEM-150","blaTEM-168","blaTEM-176","blaTEM-215","blaTEM-40","site","collection_date")]
ampicillin_resistance$ampicillin<-gsub("M|P","I", ampicillin_resistance$ampicillin)
cefotaxime_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","cefotaxime","blaBIL-1","blaCMY-107","blaCMY-132","blaCMY-16","blaCMY-2","blaCMY-40","blaCMY-7","blaCTX-M-1","blaCTX-M-10","blaCTX-M-101","blaCTX-M-14","blaCTX-M-144","blaCTX-M-15","blaCTX-M-183","blaCTX-M-211","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-65","blaCTX-M-9","blaDHA-1","blaIMP-4", "blaNDM-1", "blaOXA-10","blaSHV-12","blaTEM-106","blaTEM-12","blaTEM-168","site","collection_date")]
cefotaxime_resistance$cefotaxime<-gsub("M|P","I", cefotaxime_resistance$cefotaxime)
meropenem_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","meropenem","blaIMP-4","blaNDM-1","site","collection_date")]
meropenem_resistance$meropenem<-gsub("M|P","I", meropenem_resistance$meropenem)
azithromycin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","azithromycin","ere(A)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)","site","collection_date")]
azithromycin_resistance$azithromycin<-gsub("M|P","I", azithromycin_resistance$azithromycin)
azithromycin_resistance_novel<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","azithromycin","ere(A)","erm(42)","erm(A)","erm(B)","mef(B)","mph(A)","msr(E)","site","collection_date")]
azithromycin_resistance_novel$azithromycin<-gsub("M|P","I", azithromycin_resistance_novel$azithromycin)
chloramphenicol_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","chloramphenicol","catA1","catA2","catB3","cmlA","cmlA1","cmlA5","floR","site","collection_date")]
chloramphenicol_resistance$chloramphenicol<-gsub("M|P","I", chloramphenicol_resistance$chloramphenicol)
tetracycline_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","tetracycline","tet(A)","tet(B)","tet(C)","tet(D)","tet(G)","tet(H)","tet(M)","tet(X)","site","collection_date")]
tetracycline_resistance$tetracycline<-gsub("M|P","I", tetracycline_resistance$tetracycline)
ciprofloxacin_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","ciprofloxacin","aac(6')-Ib-cr","oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A","site","collection_date")]
ciprofloxacin_resistance$ciprofloxacin<-gsub("I|P","D", ciprofloxacin_resistance$ciprofloxacin)
nalidixicacid_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","nalidixic_acid","aac(6')-Ib-cr","oqxA","oqxA2","oqxB","qepA2","qnrA1","qnrB1","qnrB19","qnrB2","qnrB4","qnrB6","qnrB7","qnrS1","qnrS13","qnrS2","gyrA.D87G","gyrA.D87N","gyrA.D87Y","gyrA.S83A","gyrA.S83F","gyrA.S83I","gyrA.S83Y","gyrB.E466D","parC.E84G","parC.S80I","parC.S80R","parC.T57S","parE.S458A","site","collection_date")]
nalidixicacid_resistance$nalidixic_acid<-gsub("I|P","D", nalidixicacid_resistance$nalidixic_acid)
nalidixicacid_resistance_novel<-merge(nalidixicacid_resistance,ariba_gyrpar_novel, by = "ID", all=TRUE)
sulfonamides_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","sulphathiozole","sul1","sul2","sul3","site","collection_date")]
sulfonamides_resistance$sulphathiozole<-gsub("M|P","I", sulfonamides_resistance$sulphathiozole)
trimethoprim_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","trimethoprim","dfrA1","dfrA10","dfrA12","dfrA14","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","site","collection_date")]
trimethoprim_resistance$trimethoprim<-gsub("M|P","I", trimethoprim_resistance$trimethoprim)
cotrimoxazole_resistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","cotrimoxazole","dfrA1","dfrA10","dfrA12","dfrA14","dfrA15","dfrA16","dfrA17","dfrA19","dfrA23","dfrA27","dfrA32","dfrA5","dfrA7","dfrB4","sul1","sul2","sul3","site","collection_date")]
cotrimoxazole_resistance$cotrimoxazole<-gsub("M|P|I","D", cotrimoxazole_resistance$cotrimoxazole)
rifamycin_genoresistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","arr","arr-2","arr-3","site","collection_date")]
bleomycin_genoresistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","ble","bleO","site","collection_date")]
colistin_genoresistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","mcr-1.1","mcr-3.1","mcr-3.11","mcr-3.2","mcr-5.1","site","collection_date")]
lincosamide_genoresistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","lnu(C)","lnu(F)","lnu(G)","site","collection_date")]
fosfomycin_genoresistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","fosA3","fosA4","fosA7","site","collection_date")]
qac_genoresistance<-matched_merged_gp_files[,c("ID", "MLST_Scheme","serovar","ST","all_genes","genes_confirmed","genes_unconfirmed","qacE","qacEdelta1","qacG2","qacL","site","collection_date")]
                                                    
#print
write.csv(gentamicin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/gentamicin_resistance.csv",row.names=FALSE, quote=TRUE)
write.csv(kanamycin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/kanamycin_resistance.csv",row.names=FALSE, quote=TRUE)
write.csv(spectinomycin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/spectinomycin_resistance.csv",row.names=FALSE, quote=TRUE)
write.csv(streptomycin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/streptomycin_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(ampicillin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/ampicillin_resistance.csv",row.names=FALSE, quote=TRUE)
write.csv(cefotaxime_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/cefotaxime_resistance.csv",row.names=FALSE, quote=TRUE)
write.csv(meropenem_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/meropenem_resistance.csv",row.names=FALSE, quote=TRUE)
write.csv(azithromycin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/azithromycin_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(azithromycin_resistance_novel, file="~/Documents/PhD/Chapter_2/amr_files/r_output/azithromycin_resistance_novel.csv", row.names=FALSE, quote=TRUE)
write.csv(tetracycline_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/tetracycline_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(chloramphenicol_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/chloramphenicol_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(ciprofloxacin_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/ciprofloxacin_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(ciprofloxacin_resistance_novel, file="~/Documents/PhD/Chapter_2/amr_files/r_output/ciprofloxacin_resistance_novel.csv", row.names=FALSE, quote=TRUE)
write.csv(nalidixicacid_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/nalidixicacid_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(nalidixicacid_resistance_novel, file="~/Documents/PhD/Chapter_2/amr_files/r_output/nalidixicacid_resistance_novel.csv", row.names=FALSE, quote=TRUE)
write.csv(sulfonamides_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/sulfonamides_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(trimethoprim_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/trimethoprim_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(cotrimoxazole_resistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/cotrimoxazole_resistance.csv", row.names=FALSE, quote=TRUE)
write.csv(qac_genoresistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/qac_genoresistance.csv", row.names=FALSE, quote=TRUE)
write.csv(fosfomycin_genoresistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/fosfomycin_genoresistance.csv", row.names=FALSE, quote=TRUE)
write.csv(rifamycin_genoresistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/rifamycin_genoresistance.csv", row.names=FALSE, quote=TRUE)
write.csv(colistin_genoresistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/colistin_genoresistance.csv",row.names=FALSE, quote=TRUE)
write.csv(lincosamide_genoresistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/lincosamide_genoresistance.csv", row.names=FALSE, quote=TRUE)
write.csv(bleomycin_genoresistance, file="~/Documents/PhD/Chapter_2/amr_files/r_output/bleomycin_genoresistance.csv", row.names=FALSE, quote=TRUE)
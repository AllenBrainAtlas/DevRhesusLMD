## Code for generating figures 1d, 4c, and 5; extended data figures 5 and 6, and supplemental table 8 in the NHP paper
# Init workspace
if (!grepl("reports$", getwd())) setwd("reports")

##############################################################################################
##############################################################################################
print("Read in and properly format the data.")

outputFolder = "../reports/"
setwd(outputFolder)
load("../cache/nhp_PrePost_StartingData.RData")
codeFolder   = "../src/" 
inputFolder  = "../cache/inputFiles/"

# LIBRARIES
library(WGCNA);
library(reshape);
library(gplots);
library(flashClust);

# Files with other functions
source(paste(codeFolder,"functionsForMacaqueAnalysis.r",sep=""))

useThese   = probes$is_best_probeid&probes$is_passing_probeid  # (This is all genes, but is kept in for compatibility)
bestProbes = probes$probeid[useThese]
bestGenes  = probes$macaque_genesymbol[useThese]
humanGenes = probes$human_genesymbol[useThese]
humanGenes[humanGenes==bestGenes] = NA
bestEntrez = probes$macaque_entrezid[useThese]
exprPP     = exprPrePost[useThese,]
samplePP   = samplePrePost
exprPPG    = exprPP; 
rownames(exprPP) = bestProbes;
medianVals = 2^apply(exprPP,1,median) 
ages       = levels(samplePP$age)

print("-----For areas with replicates, take the average of the replicates as the value for that area.")
regLayDonor = paste(samplePP$struc_age,samplePP$donor_name,sep="_")
exprPPGrep  = findFromGroups(t(exprPPG),regLayDonor)
kpRep       = rep(TRUE, length(regLayDonor))
for (i in 2:length(regLayDonor))
  if(length(intersect(regLayDonor[i],regLayDonor[1:(i-1)])>0)) kpRep[i] = FALSE
exprPPGrep  = exprPPGrep[,regLayDonor[kpRep]]
samplePPrep = samplePP[kpRep,]
colnames(exprPPGrep) = rownames(samplePPrep)
ord         = order(samplePPrep$age_log2pcd,samplePPrep$toplevel,samplePPrep$region,samplePPrep$layer)
exprPPGrep  = exprPPGrep[,ord]
samplePPrep = samplePPrep[ord,] 


##############################################################################################
##############################################################################################
print("Code for Figure 1d.")

layerPositionsS2 = read.csv(paste(inputFolder,"layerRectanglePositions_small2.csv",sep=""),row.names=1)
agePositionsS    = read.csv(paste(inputFolder,"ageRectanglePositions_small.csv",sep=""),row.names=1)
layerPlot        = as.character(samplePPrep$layer)
layerPlot[layerPlot=="L2-3"] = "CPo"
names(layerPlot) = rownames(samplePPrep)
fig1dGenes       = c("PAX6","EOMES","DCX","SYT1","GFAP","MOBP")
kpV1             = (samplePPrep$subregion=="V1")&(!is.element(layerPlot,c("TMZ","ICD")))
# Only plot data from V1 for figure 1d, and exclude a few samples from minor fiber layers.

pdf("Figure1d_GenePlots.pdf",width=7,height=7)
for (a in c(fig1dGenes,"IL7"))
  plotMacaqueCortexSmall(exprPPGrep[a,][kpV1],layerPlot[kpV1],samplePPrep$age[kpV1],layerPositionsS2,
    agePositionsS,plotTitle=paste("Mean expr for M01, V1 -",a),quantileScale=c(0.05,0.95),isLog2=TRUE,legendPos = c(8,-13,-10))
dev.off()
# IL7 is only plotted to get the legend for the figure


print("Code for Extended Data Figure 2 -- These are expanded versions of the plots for Figure 1d.")

layerPositions   = read.csv(paste(inputFolder,"layerRectanglePositions.csv",sep=""),row.names=1)
regionPositions  = read.csv(paste(inputFolder,"regionRectanglePositions.csv",sep=""),row.names=1)
ageOffsets       = read.csv(paste(inputFolder,"ageRectangleOffsets.csv",sep=""),row.names=1)

pdf("ExtendedDataFigure2_GenePlots.pdf",width=18,height=9)
for (g in fig1dGenes) {
 i  = which(bestGenes==g);
 hg = humanGenes[i];   
 if (!is.na(hg)) g = paste(g,hg,sep="__")
 p  = bestProbes[i];
 e  = bestEntrez[i];
 fn = paste(outputFolder,"fig_1d/",g,"__",e,"__",p,".pdf",sep="")
 plotMacaqueCortex(exprPP[p,],NULL,samplePP$layer,samplePP$subregion,samplePP$age,
    layerPositions,regionPositions,ageOffsets,paste(g,e,p,sep=" - "),
    quantileScale=c(0.1,0.95),medianVals=medianVals)
}
dev.off()


print("NOTE: Uncomment this section of the code to generate plots like the one in Extended Data Figure 2 for ALL passing genes.")
#print("Make plots like the one in Extended Data Figure 2 for ALL passing genes used in the analysis.")

#dir.create(paste(outputFolder,"allGenePlots/",sep=""))
#for (g in bestGenes) {
# i  = which(bestGenes==g);
# hg = humanGenes[i];   
# if (!is.na(hg)) g = paste(g,hg,sep="__")
# p  = bestProbes[i];
# e  = bestEntrez[i];
# fn = paste(outputFolder,"allGenePlots/",g,"__",e,"__",p,".pdf",sep="")
# pdf(fn,width=18,height=9)
# plotMacaqueCortex(exprPP[p,],NULL,samplePP$layer,samplePP$subregion,samplePP$age,
#    layerPositions,regionPositions,ageOffsets,paste(g,e,p,sep=" - "),
#    quantileScale=c(0.1,0.95),medianVals=medianVals)
# dev.off()
#}


##############################################################################################
##############################################################################################
print("Code for Figure 3c and extended data figure 5d.")

gnsF3  = c("MOBP","OPALIN")
gnsED5 = c("MOG","ERMN","MAL","ASPA")
gns    = c(gnsF3,gnsED5)
wmAcro = as.character(samplePPrep$struc_acro)
wmAcro[is.element(wmAcro,c("GPi","GPo"))]   = "GP"
wmAcro[is.element(wmAcro,c("V1wm","V1iz"))] = "V1wm/iz"
wmAcro[!is.element(wmAcro,c("V1wm/iz","GP","ic","CA1or"))] = "other"

regs    = c("CA1or","V1wm/iz","ic","GP")
datIn   = list()
for (r in regs)
 datIn[[r]] = findFromGroups(t(exprPPGrep[gns,wmAcro==r]),samplePPrep$age[wmAcro==r])
ages    = levels(samplePPrep$age)
agePlot = ages2pcd(ages)
cols2   = c("cyan","green","purple","red")
cols    = c("blue","green","purple","red")

pdf(paste(outputFolder,"Figure3c_OligodendrocyteMarkerGenes.pdf",sep=""),height=6,width=6);
for(gn in gnsF3){      
 plot(0,0,col="white",xlim=c(1,10),ylim=c(2,14),ylab="Gene Expression",
   xlab="Age (log2(pcd))",main=paste(gn,"across WM regions"))
 for (i in 1:length(regs)){   r = regs[i]
   lines(1:10,datIn[[r]][gn,],col=cols2[i],lwd=3);
   text(8+0.3*i,5,r,col=cols[i],srt=90)
 }
 abline(v=6.5); points(1:10,rep(2.5,10),pch=19,cex=0.5)
}
dev.off()

pdf(paste(outputFolder,"ExtendedDataFigure5d_OligodendrocyteMarkerGenes.pdf",sep=""),height=6,width=6);
for(gn in gnsED5){      
 plot(0,0,col="white",xlim=c(1,10),ylim=c(2,14),ylab="Gene Expression",
   xlab="Age (log2(pcd))",main=paste(gn,"across WM regions"))
 for (i in 1:length(regs)){   r = regs[i]
   lines(1:10,datIn[[r]][gn,],col=cols2[i],lwd=3);
   text(8+0.3*i,5,r,col=cols[i],srt=90)
 }
 abline(v=6.5); points(1:10,rep(2.5,10),pch=19,cex=0.5)
}
dev.off()


##############################################################################################
##############################################################################################
print("Code for Figure 4a ==> comparison of layer patterning between each pair of ages")

# Load pre-generated laminar marker genes (from "nhp_laminar_analysis.R")

orderedLaminarGenes = read.table(paste(inputFolder,"V1(L1-WM)_max_laminar_genes_byage.txt",sep=""),header=TRUE)
colnames(orderedLaminarGenes) = gsub("X","",colnames(orderedLaminarGenes))

agesLA   = c("E80","E90","E120","0M","3M","12M","48M");
layersLA = c("L1","L2_3","L4","L5","L6","WM")  
ageLA    = rep(agesLA,length(layersLA));
layerLA  = rep(layersLA[1],length(agesLA))
for (i in 2:length(layersLA))
  layerLA  = c(layerLA, rep(layersLA[i],length(agesLA)))
ageLayerLA = paste(ageLA,layerLA,sep="_")
orderedLaminarGenes = orderedLaminarGenes[,ageLayerLA]

kpV1     = (samplePP$subregion=="V1")&(samplePP$age_log2pcd>6)&(samplePP$layer!="CP")&
           is.element(samplePP$layer_dev,c("L2","L3","L4","L5","L6"))
whMax    = function(x,nm) return(nm[which.max(x)]) 
meanV1   = list()
layerTmp = as.character(samplePP$layer_dev)
layerTmp[is.element(layerTmp,c("L2","L3"))] = "L2_3"
for (a in agesLA){       
 kpV            = kpV1&(samplePP$age==a);
 meanV1[[a]]    = findFromGroups2(t(exprPPG[,kpV]),layerTmp[kpV])
}

N = 1244 # Find the top 10% (1244) laminar genes at each time point.
distAges = matrix(0,nrow=7,ncol=7);  rownames(distAges) <- colnames(distAges) <- agesLA
meanCorsV2 <- distAgesN <- distAgesE <- distAgesC <- distAgesT <- distAges
for (a1 in agesLA) for (a2 in agesLA) { 
 lamGenes = NULL
 for (l in layersLA) {
  nm = paste(a2,l,sep="_")
  lamGenes = c(lamGenes,as.character(orderedLaminarGenes[1:round(N/length(layersLA)),nm]))
 }
 isLaminarV        = is.element(rownames(meanV1[[a1]]),lamGenes)
 meanCorsV2[a1,a2] = mean(diag(cor(t(meanV1[[a1]][isLaminarV,]),t(meanV1[[a2]][isLaminarV,]))))
 distAges[a1,a2]   = paste(a1,a2,sep="_")
 distAgesN[a1,a2]  = abs(ages2pcd(a1)-ages2pcd(a2))
 distAgesC[a1,a2]  = abs(which(ages==a1)-which(ages==a2))
 distAgesE[a1,a2]  = abs(ages2inv(a1)-ages2inv(a2))
 distAgesT[a1,a2]  = abs(convertAge(ages2pcd(a1),speciesOut="eventScore")- 
                                      convertAge(ages2pcd(a2),speciesOut="eventScore"))
}

pdf("Figure4a_LaminarSimilarityBetweenAges_top10percent.pdf",width=6.5,height=7);
t  = as.vector(distAges);
y  = as.vector(meanCorsV2);
k  = y<0.99
xC = as.vector(distAgesC);
xN = as.vector(distAgesN);
xE = as.vector(distAgesE);
xT = as.vector(distAgesT);
verboseScatterplot(xE[k],y[k],xlim=range(distAgesE),ylim=c(0,1),
  main="",ylab="Laminar expression correlation",abline=TRUE,
  xlab="Age difference (delta pcd -1)",col="white");
tmpCols = colorRampPalette(c("blue","red"))(8)[2:8];  names(tmpCols) = agesLA
for (i in 1:length(t[k])){
 tmpAges = strsplit(t[k][i],"_")[[1]]
 points(xE[k][i],y[k][i],cex=3,pch=19,col=tmpCols[tmpAges[2]])
 points(xE[k][i],y[k][i],cex=3,pch=1,col="black")
 points(xE[k][i],y[k][i],cex=1.5,pch=19,col=tmpCols[tmpAges[1]])
}
points(rep(0.009,7),(0.5+(0:6)*0.05),cex=3,pch=19,col=tmpCols[7:1])
text(rep(0.0102,7),(0.5+(0:6)*0.05),names(tmpCols)[7:1],cex=1.5)
points(rep(0.009,2),(0.4+(0:1)*0.05),cex=c(1.5,3),pch=19)
text(rep(0.0102,2),(0.4+(0:1)*0.05),c("Comparison","Source list"),cex=0.8)
dev.off()


print("Make the plot for figure 4c (layer 6 marker genes)")
print("-----Note: Only part of each of these plots were used in making figure 4e.")

layerPositionsS = read.csv(paste(inputFolder,"layerRectanglePositions_small.csv",sep=""),row.names=1)
agePositionsS   = read.csv(paste(inputFolder,"ageRectanglePositions_small.csv",sep=""),row.names=1)

layer6gns  = c("SLC26A7","LGR5","CA12","SYT6")
kpV1       = samplePP$subregion=="V1"
agePlot    = as.character(samplePP$age)[kpV1]
layerPlot  = as.character(samplePP$layer)[kpV1];     layerPlot[layerPlot=="L2-3"] = "CPo"
layerPlot[is.element(layerPlot,c("L4Cb","L4Ca","L4B","L4A"))] = "L4"
kpACG      = samplePP$subregion=="ACG"
agePlotA   = as.character(samplePP$age)[kpACG]
layerPlotA = as.character(samplePP$layer)[kpACG];     layerPlotA[layerPlotA=="L2-3"] = "CPo"
layerPlotA[is.element(layerPlotA,c("L4Cb","L4Ca","L4B","L4A"))] = "L4"
legendPos  = c(8,-11.5,-8.5);

pdf("Figure4c_MarkerGenesForLayer4.pdf",width=7,height=7)
for (g in layer6gns) {
 i  = which(bestGenes==g);
 hg = humanGenes[i];
 if (!is.na(hg)) g = paste(g,hg,sep="__")
 p  = bestProbes[i];
 e  = bestEntrez[i];
 plotMacaqueCortexSmall(exprPP[p,kpV1],layerPlot,agePlot,layerPositionsS,agePositionsS,
  paste(g,e,p,"V1",sep=" - "),isLog2=TRUE,combineFn="meanNA",quantileScale=c(0.1,0.95),
  linearOrLog="linear",bgPar="grey95",displayLayers=FALSE,legendPos=legendPos)
 abline(v=c(0,6,10),lwd=2);   abline(h=c(0,-8,-12),lwd=2)
}
dev.off();


##############################################################################################
##############################################################################################
print("Code for Figure 4d ==> Marker genes for V1 and ACG at different time points.")
print("Note: The remainder of this figure combines replicate samples as in figure 1d and 4c.")
print("Note: All plots in the rest of this figure have V1 as some shade of blue and ACG as some.")
print("      shade of red, which are then matched in photoshop.")
print("Note: The line plots in Figure e-h are generated in a different code document.")

# Prepare the data for this figure
layKp3   = c("L2","L3","L2-3","L5","L6","CPi","CPo","VZ","VZi","VZo","SZ","SZi","SZo")
kpSamp3  = is.element(samplePPrep$subregion,c("ACG","V1","S1"))&is.element(samplePPrep$layer,layKp3)
age3     = as.character(samplePPrep$age)[kpSamp3];
layAll3  = as.character(samplePPrep$layer)[kpSamp3];
layCol3  = rep("L2/L3",sum(kpSamp3))
layCol3[is.element(layAll3,c("L5","L6","CPi"))] = "L5/L6"
layCol3[is.element(layAll3,c("SZ","SZi","SZo"))] = "SZ"
layCol3[is.element(layAll3,c("VZ","VZi","VZo"))] = "VZ"
region3  = as.character(samplePPrep$subregion)[kpSamp3];
laAgDon3 = paste(layCol3,age3,as.character(samplePPrep$donor)[kpSamp3])
laAgDons = intersect(laAgDon3[region3=="ACG"],laAgDon3[region3=="V1"])
laAgDons = setdiff(laAgDons,"L2/L3 E70 DAM35650") # This sample is of a different format from the rest, so it is removed
exprG3   = exprPPGrep[,kpSamp3]
exprA3   = findFromGroups(t(exprG3[,region3=="ACG"]),laAgDon3[region3=="ACG"])
exprA3   = exprA3[,laAgDons]
exprV3   = findFromGroups(t(exprG3[,region3=="V1"]),laAgDon3[region3=="V1"])
exprV3   = exprV3[,laAgDons]
exprS3   = findFromGroups(t(exprG3[,region3=="S1"]),laAgDon3[region3=="S1"])
exprS3   = exprS3[,intersect(colnames(exprS3),laAgDons)]
age3av   = age3[match(laAgDons,laAgDon3)]
layAll3av= layAll3[match(laAgDons,laAgDon3)]
layCol3av= layCol3[match(laAgDons,laAgDon3)]
lays3    = c("L2/L3","L5/L6","SZ","VZ")

# Run all of the t-tests (note that the p-values are uncorrected for multiple comparisons)
fc3 <- pval3 <- list()
for (l in lays3){
 print(paste("---------- Working on layer",l))
 pval3[[l]] = matrix(NA,nrow=dim(exprG3)[1],ncol=length(ages))
 colnames(pval3[[l]]) = ages;  rownames(pval3[[l]]) = rownames(exprG3)
 fc3[[l]]   = pval3[[l]]
 for (a in ages){
  kp = (age3av==a)&(layCol3av==l)
  if (sum(kp)>=3){
   pval3[[l]][,a] = apply(cbind(exprA3[,kp],exprV3[,kp]),1, t.test.l.paired2)
   fc3[[l]][,a]   = rowMeans(exprA3[,kp])-rowMeans(exprV3[,kp])
  }
 }
}


# Make the plot for figure 4d
pdf(paste(outputFolder,"Figure4d_NumberOfV1andACGgenes_perAgeAndLayer.pdf",sep=""),height=15,width=7);
par(mfrow=c(4,1))
pThresh  = 0.05  # Uncorrected p-value
fcThresh = 1     # But also a fold change correction
la       = length(ages)
mx       = 1000
for (l in lays3){
 numV = colSums((pval3[[l]]<pThresh)&(fc3[[l]]< -fcThresh))
 numA = colSums((pval3[[l]]<pThresh)&(fc3[[l]]>  fcThresh))
 plot(1:la,numV,ylim=c(-0.2*mx,mx),main= paste("Number of DEX genes in",l),xlim=c(0.5,la+0.5),
  xlab="",ylab="Number of genes",col="blue",type="l",lwd=3);  abline(h=0);
 lines(1:la,numA,col="red",lwd=3)
 text(1:10,rep(-0.1*mx,10),ages,cex=2)
} 
dev.off()


print("Output Supplemental Table 8, which includes all of these significant genes")
gnR <- layR <- ageR <- avR <- fcR <- pR <- NULL
for (l in lays3) for (a in ages[!is.na(pval3[[l]][1,])]){
 f    = as.numeric(fc3[[l]][,a]);     p   = as.numeric(pval3[[l]][,a])
 kpV  = (p<pThresh)&(f< -fcThresh);   kpA = (p<pThresh)&(f>fcThresh);  
 gnR  = c(gnR,bestGenes[kpV],bestGenes[kpA])
 layR = c(layR,rep(l,sum(kpV)+sum(kpA)))
 ageR = c(ageR,rep(a,sum(kpV)+sum(kpA)))
 avR  = c(avR,rep("V1",sum(kpV)),rep("ACG",sum(kpA)))
 fcR  = c(fcR,-f[kpV],f[kpA])
 pR   = c(pR,p[kpV],p[kpA])
}
regOutZZ = data.frame(Gene=gnR,Layer=layR,Age=ageR,Region=avR,log2FoldChange=fcR,Pvalue=pR)
regOutZZ = regOutZZ[order(layR,factor(ageR,levels=ages),avR,pR,gnR),]
write.csv(regOutZZ,"SupplementalTable8_V1vsACG_markerGenes.csv",row.names=FALSE)


##############################################################################################
##############################################################################################
print("Code for Extended Data Figure 6a ==> Compare expression levels in S1 with V1 and S1.")
print("-----Specifically, find genes higher in ACG in at least 2 of E70, E80, and E90, then check pattern in S1.")
print("-----Some of this is now part of a Supplemental Figure.")
print("Note: code for Extended Data Figure 5b is found in a separate R document.")

gnA <- list()
pdf("ExtendedDataFigure6a_V1andACGvsS1.pdf",height=6,width=17)
par(mfrow=c(1,3))
for (l in c("L5/L6","SZ","VZ")){
  gnA[[l]] = sort(rownames(pval3[[l]])[rowSums((pval3[[l]][,3:5]<pThresh)&(fc3[[l]][,3:5]>fcThresh))>=2])
  kpAV =  colnames(exprA3)[is.element(age3av,c("E70","E80","E90"))&(layCol3av==l)]
  kpS  =  intersect(kpAV,colnames(exprS3))
  reg  =  factor(c(rep("ACG",length(kpAV)),rep("S1",length(kpS)),rep("V1",length(kpAV))),levels=c("ACG","S1","V1"))
  exprASV = findFromGroups(t(cbind(exprA3[,kpAV],exprS3[,kpS],exprV3[,kpAV])[gnA[[l]],]-rowMeans(exprV3[gnA[[l]],kpAV])),reg) 
  reg  = factor(c(rep("ACG",dim(exprASV)[1]),rep("S1",dim(exprASV)[1])),levels=c("ACG","S1"))
  dat  = c(exprASV[,1],exprASV[,2])
  verboseBoxplot(dat,reg,notch=FALSE,ylab="log2FC (region vs. V1)",main=paste("Layer",l),ylim=c(-1,5))
  abline(h=0)
  verboseScatterplot(exprASV[,1],exprASV[,2],xlab="log2FC (ACG vs. V1)", 
    ylab="log2FC (S1 vs. V1)",pch=19,cex=0.5,xlim=c(0,5),ylim=c(-1,4))
  abline(0,1); abline(0,0.5); abline(h=0)
  points(median(exprASV[,1]),median(exprASV[,2]),cex=3.2,pch=19)
  verboseScatterplot(exprASV[,1],exprASV[,2],xlab="log2FC (ACG vs. V1)", 
    ylab="log2FC (S1 vs. V1)",col="white",xlim=c(0,5),ylim=c(-1,4))
  abline(0,1,col="grey"); abline(0,0.5,col="grey"); abline(h=0,col="grey")
  text(exprASV[,1],exprASV[,2],gnA[[l]],cex=0.5)
}
dev.off()




##############################################################################################
##############################################################################################
print("Code for Figure 4e-h ==> Perform (and plot) enrichment analyses for cell type marker genes.")
print("Note: Marker genes for inhibitory and excitatory neurons were identified from Tasic et. al 2016.  Markers for")
print("      cell cycle (likely reflecting progenitor cells herein) were identified from Bar Joseph et. al 2008.")
print("      Markers for astrocytes were identified from Cahoy et. al 2008.  All markers are included in")
print("      userInputMarkers_updatedTasic.csv.  Marker assignment strategies are described in Supplemental Methods.")
print("Note: The legend was generated by comparing the size of the boxes in these plots to the p-values output")
print("      from userListEnrichment.")

enrichs3V <- enrichs3A <- NULL
genes3   = rownames(pval3[[l]])
fnIn     = paste(inputFolder,"userInputMarkers_updatedTasic.csv",sep="") 
 
for (l in lays3) for (a in ages){
 lala   = paste(l,a)
 label3 = rep("grey",length(genes3))
 label3[(pval3[[l]][,a]<pThresh)&(fc3[[l]][,a]< -fcThresh)] = lala
 if(sum(label3==lala)>=5){
  enrichsTmp = userListEnrichment(genes3,label3,catNmIn="input",fnIn=fnIn)
  enrichsTmp = enrichsTmp$sigOverlaps[,c(1,2,4)]
  enrichsTmp = as.matrix(enrichsTmp)
  enrichs3V  = rbind(enrichs3V,enrichsTmp)
 }
}
enrichs3V = as.data.frame(enrichs3V)
enrichs3V$CorrectedPvalues = as.numeric(as.character(enrichs3V$CorrectedPvalues))
for (l in lays3) for (a in ages){
 lala   = paste(l,a)
 label3 = rep("grey",length(genes3))
 label3[(pval3[[l]][,a]<pThresh)&(fc3[[l]][,a]>fcThresh)] = lala
 if(sum(label3==lala)>=5){
  enrichsTmp = userListEnrichment(genes3,label3,catNmIn="input",fnIn=fnIn)
  enrichsTmp = enrichsTmp$sigOverlaps[,c(1,2,4)]
  enrichsTmp = as.matrix(enrichsTmp)
  enrichs3A  = rbind(enrichs3A,enrichsTmp)
 }
}
enrichs3A = as.data.frame(enrichs3A)
enrichs3A$CorrectedPvalues = as.numeric(as.character(enrichs3A$CorrectedPvalues))
omitX = c(1,2,3,1,7,8,9,10,7,8,9,10)
omitY = c(3,3,3,2,1,1,1,1,0,0,0,0)

# Scale results for appropriate multiple comparisons  # OMIT THIS SECTION FOR THE TF ANALYSIS
pScale   = 242/4  # ALL p-values are scaled by this value to adjust for 238 additional tests not used in this paper or included in the code
enrichs3A[,3] = enrichs3A[,3]*pScale
enrichs3A = enrichs3A[enrichs3A[,3]<pThresh,]
enrichs3V[,3] = enrichs3V[,3]*pScale
enrichs3V = enrichs3V[enrichs3V[,3]<pThresh,]


# Plot enrichment analyses for figure 4e-h
pdf("Figure4efgh_CellTypeEnrichmentPlots_ACGvsV1.pdf",width=6,height=3.5)
for (cat in sort(unique(rbind(enrichs3V,enrichs3A)[,2]))){
 plot(0,0,col="white",xlim=c(-0.5,length(ages)+0.5),ylim=c(-1.5,3.5),main=cat,
   axes=FALSE,xlab="",ylab="")
 enrichTmpA = enrichs3A[enrichs3A[,2]==cat,]
 enrichTmpV = enrichs3V[enrichs3V[,2]==cat,]
 abline(v=1:length(ages),lty="dotted");  abline(h=0:3,lty="dotted"); 
 text(1:length(ages),rep(-1,length(ages)),ages)
 text(rep(0,length(lays3)),4-(1:length(lays3)), lays3)
 points(omitX,omitY,pch=19)
 for (l in lays3) for (a in ages){
  lala   = paste(l,a)
  kpV    = which(enrichTmpV$InputCategories==lala)
  kpA    = which(enrichTmpA$InputCategories==lala)
  if(length(kpV)>0){
   cexes = sqrt(as.numeric(-log10(enrichTmpV$CorrectedPvalues)))[kpV]
   if(length(cexes)>0){
    x = which(ages==a)+0.1
    y = 4-which(lays3==l)
    points(x,y,pch=15,col="blue",cex=cexes)
   }
  }
  if(length(kpA)>0){
   cexes = sqrt(as.numeric(-log10(enrichTmpA$CorrectedPvalues)))[kpA]
   if(length(cexes)>0){
    x = which(ages==a)-0.1
    y = 4-which(lays3==l)
    points(x,y,pch=15,col="red",cex=cexes)
}}}}
dev.off()


file.remove("enrichment.csv")
## End of code.
##############################################################################################
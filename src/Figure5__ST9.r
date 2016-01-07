## Code for generating figure 5 and supplemental table 9 in the NHP paper

##############################################################################################
##############################################################################################
# Init workspace
if (!grepl("reports$", getwd())) setwd("reports")

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


###################################################################################################
###################################################################################################
print("Run 10 WGCNA networks, one for each age, on replicate-averaged data sets from each age.")
print("=== These networks only include samples from V1 and ACG, the two cortical regions with data at each age.")
print(" ")
print("***** NOTE: THIS ANALYSIS TAKES A FEW HOURS AND REQUIRES A LOT OF RAM. *****")
print("I ran it on a Windows desktop, 64-bit OS, 24GB of RAM.  My *guess* is 12GB of RAM would be enough.")

# Run the actual networks
allColors  = matrix(0,nrow=length(bestGenes),ncol=10)
AGES = c("E40","E50","E70","E80","E90","E120","0M","3M","12M","48M")[10:1]
rownames(allColors) = bestGenes
colnames(allColors) = AGES
omitLabels   = NULL
matchModules = NA
kpRegion     = is.element(samplePPrep$subregion,c("ACG","V1"))
for (AGE in AGES){
 print(paste("Starting",AGE))
 print(date())
 kpSubset  = samplePPrep$age==AGE
 datExprRo = exprPPGrep[,kpSubset&kpRegion]
 if (ages2pcd(AGE)>7){
  print("===Scale V1 and ACG to minimize potential batch effects in postnatal data.")
  subVA = as.character(samplePPrep$subregion[kpSubset&kpRegion])
  useLay = is.element(as.character(samplePPrep$layer_dev)[kpSubset&kpRegion],c("L2","L3","L5","L6"))
  scaleVal = findFromGroups(t(datExprRo[,useLay]),subVA[useLay],function(x) return(mean(sort(x)[4:9])))
  scaleVal = rowMeans(scaleVal)-scaleVal
  scaleVal = scaleVal[,subVA]
  scaleVal = scaleVal[,subVA]
  datExprRo = datExprRo + scaleVal
 }
 dirName  = paste("NCX_WGCNA",AGE,sep="_")
 source(paste(codeFolder,"WGCNA_automated_code_Block.r",sep=""))
 allColors[names(moduleColors),AGE] = moduleColors
 omitLabels   = sort(unique(c(omitLabels,moduleColors)))
 matchModules = moduleColors
 print(paste("Done with",AGE))
}


print("Reorder and merge a few modules to line them up across ages as closely as possible.")
AGES           = c("E40","E50","E70","E80","E90","E120","0M","3M","12M","48M")[10:1]
omitLabels     = NULL
matchModules   = NA
modTable       = read.csv(paste(inputFolder,"modSwap.csv",sep=""))
modSwap        = as.character(modTable[,1])
names(modSwap) = modTable[,2]
for (AGE in AGES){
 print(paste("Starting",AGE))
 kpRegion  = is.element(samplePPrep$subregion,c("ACG","V1"))
 kpSubset  = samplePPrep$age==AGE
 datExprRo = exprPPGrep[,kpSubset&kpRegion]
 if (ages2pcd(AGE)>7){
  print("Scale V1 and ACG to minimize potential batch effects in postnatal data.")
  subVA = as.character(samplePPrep$subregion[kpSubset&kpRegion])
  useLay = is.element(as.character(samplePPrep$layer_dev)[kpSubset&kpRegion],c("L2","L3","L5","L6"))
  scaleVal = findFromGroups(t(datExprRo[,useLay]),subVA[useLay],function(x) return(mean(sort(x)[4:9])))
  scaleVal = rowMeans(scaleVal)-scaleVal
  scaleVal = scaleVal[,subVA]
  datExprRo = datExprRo + scaleVal
 }
 dirName  = paste("NCX_WGCNA",AGE,sep="_")
 moduleSwap = modSwap[grep(AGE,modSwap)]
 source(paste(codeFolder,"WGCNA_reorderAndMergeModules.r",sep=""))
 allColors[names(moduleColors),AGE] = moduleColors
 omitLabels   = sort(unique(c(omitLabels,moduleColors)))
 matchModules = moduleColors
 print(paste("Done with",AGE))
}


print("Collect the module assignments into a single matrix.")
AGES = c("E40","E50","E70","E80","E90","E120","0M","3M","12M","48M")
MEsL <- list()
moduleColorsL = NULL
for (AGE in AGES){
 load(paste("NCX_WGCNA_",AGE,"/WGCNA_Robject.RData",sep=""))
 MEsL[[AGE]] = MEs
 MEsL[[AGE]] = MEsL[[AGE]][,setdiff(colnames(MEsL[[AGE]]),ignoreMods)]
 moduleColorsL = cbind(moduleColorsL,moduleColors)
}
rownames(moduleColorsL) = names(moduleColors)
colnames(moduleColorsL) = AGES


print("Supplemental Table 9 ==> Output module assignments to a table.")
write.csv(moduleColorsL,"SupplementalTable9_ModuleAssignments.csv")


###################################################################################################
###################################################################################################
print("Identify and plot the module lineages (i.e., which modules have overlapping genes at adjacent ages).")
mods = unique(c(sort(moduleColorsL[,1]), sort(moduleColorsL[,2]), sort(moduleColorsL[,3]), 
       sort(moduleColorsL[,4]), sort(moduleColorsL[,5]), sort(moduleColorsL[,6]),   
	   sort(moduleColorsL[,7]), sort(moduleColorsL[,8]), sort(moduleColorsL[,9]), sort(moduleColorsL[,10])))
mods = mods[substr(mods,nchar(mods)-2,nchar(mods))!="M00"]
ageMods = substr(mods,1,nchar(mods)-4)
xMods = match(ageMods,AGES)
names(xMods) <- names(ageMods) <- mods
yMods <- countMods <- xMods*0
for (a in AGES){
 kpMods = mods[ageMods==a]
 l = length(kpMods)+1
 yMods[kpMods] = 1+9*(1:(l-1))/l
 countMods[kpMods] = table(moduleColorsL[,a])[kpMods]
}
yMods = 11-yMods
shortMods = as.character(lapply(as.list(mods),function(x) return(strsplit(x,"_")[[1]][2])))
shortMods = gsub("M0","",shortMods)
shortMods = gsub("M","",shortMods)

segX1 <- segX2 <- segY1 <- segY2 <- segW <- NULL
pThresh=50
ignoreMods = paste(AGES,"M00",sep="_")
for (a in 2:10){
 a1      = AGES[a-1];  a2 = AGES[a]
 ovTable = overlapTable(moduleColorsL[,a1],moduleColorsL[,a2],ignore=ignoreMods)
 pTable  = -log10(ovTable$pTable+10^(-99));
 pTable2 = pmin(pTable,100)
 pTable2[pTable2<pThresh] = NA
 pTable2 = sqrt(pTable2)/4
 for (m1 in rownames(pTable2)) for (m2 in colnames(pTable2)) if(!is.na(pTable2[m1,m2])){
  segX1 = c(segX1,xMods[m1])
  segY1 = c(segY1,yMods[m1])
  segX2 = c(segX2,xMods[m2])
  segY2 = c(segY2,yMods[m2])
  segW  = c(segW,pTable2[m1,m2])
 }
}

strucCol = rep("",length(mods))
names(strucCol) <- mods
strucCol[mods]  = "grey"

pdf("Figure5a_ModuleSizesAndGeneOverlap.pdf",height=9,width=7)
plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="")
text(1:10,rep(0,10),AGES,srt=90)
points(xMods,yMods,pch=1,cex=sqrt(countMods)/10,col="grey")#"green")
segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
text(xMods,yMods,shortMods,cex=1)
dev.off()
# This is the baseline plot for Figure 5. The size of each circle repersents the number of genes 
# in that module and connected modules have a highly significant number of overlapping genes.


###################################################################################################
###################################################################################################
print("Determine which modules have higher expression in cortical plate layers and WM/progenitor layers.")

pdf("Figure5a_PostmitoticNeuronVsOtherLayers.pdf",height=9,width=7)
strucCol2  = strucCol
strucCol2[names(strucCol2)] = "cyan"
fcThresh  = log2(1.0)  # 0.5
pThresh   = 10
layerPlot = as.character(samplePPrep$layer)
layerPlot[is.element(layerPlot,c("L2_3","CPo","L2","L3","L2-3","CP"))] = "L2/3"
layerPlot[is.element(layerPlot,c("L4A","L4B","L4Ca","L4Cb"))] = "L4"
layerPlot[is.element(layerPlot,c("CPi","L5","L6"))] = "L5/6"
names(layerPlot) = rownames(samplePPrep)
cpLayer = c("CPi","L6","L5","CPo","CP","L4","L3","L2","L4Cb","L4Ca","L4B","L4A","L2/3","L5/6","SP")
for (a in ages) {
 kp       = is.element(colnames(exprPPGrep),rownames(MEsL[[a]]))
 meanExpr = findFromGroups(t(exprPPGrep[,kp]),is.element(layerPlot[kp],cpLayer))
 laysTmp  = colnames(meanExpr)
 for (m in mods[ageMods==a]){
  gnM   = rownames(moduleColorsL)[moduleColorsL[,a]==m]
  kpTmp = diff(colMeans(meanExpr[gnM,]))
  if(kpTmp<0) strucCol2[m] = "orange"
  if(abs(kpTmp)<fcThresh) strucCol2[m] = "grey"
 }
}
plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="",main="CP enrichment")
text(1:10,rep(0,10),AGES,srt=90)
points(xMods,yMods,pch=19,cex=5,col=strucCol2)
segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
dev.off()


###################################################################################################
###################################################################################################
print("Plot enrichment for each layer.")

pdf("Figure5a_MaximumLayerExpression.pdf",height=9,width=7)
layerPlot = as.character(samplePPrep$layer)
layerPlot[is.element(layerPlot,c("L2_3","CPo","L2","L3","L2-3","CP"))] = "L2/3"
layerPlot[is.element(layerPlot,c("L4A","L4B","L4Ca","L4Cb"))] = "L4"
layerPlot[is.element(layerPlot,c("CPi","L5","L6"))] = "L5/6"
names(layerPlot) = rownames(samplePPrep)
useLayers = c("VZi","VZo","VZ","SZi","SZo","SZ","WM","IZ","SP","L5/6","L4","L2/3","L1","MZ")
topLay = mods
names(topLay) = mods
topLay[mods] = " "
fcThresh = 0
for (a in ages) {
 kp       = is.element(colnames(exprPPGrep),rownames(MEsL[[a]]))&(is.element(layerPlot,useLayers))
 meanExpr = findFromGroups(t(exprPPGrep[,kp]),layerPlot[kp])
 laysTmp  = colnames(meanExpr)
 for (m in mods[ageMods==a]){
  gnM   = rownames(moduleColorsL)[moduleColorsL[,a]==m]
  exprTmp = colMeans(meanExpr[gnM,])
  diffFC  = diff(sort(-exprTmp)[1:2])
  if(diffFC>fcThresh) topLay[m] = laysTmp[which.max(exprTmp)[1]]
 }
}
plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="",main="Maximum Layer (split VZ/SZ; all)")
text(1:10,rep(0,10),AGES,srt=90)
points(xMods,yMods,pch=19,cex=5,col=strucCol2)
segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
text(xMods,yMods,topLay)

layerPlot[is.element(layerPlot,c("VZi","VZo"))] = "VZ"
layerPlot[is.element(layerPlot,c("SZi","SZo"))] = "SZ"
topLay = mods
names(topLay) = mods
topLay[mods] = " "
fcThresh = 0.25
for (a in ages) {
 kp       = is.element(colnames(exprPPGrep),rownames(MEsL[[a]]))&(is.element(layerPlot,useLayers))
 meanExpr = findFromGroups(t(exprPPGrep[,kp]),layerPlot[kp])
 laysTmp  = colnames(meanExpr)
 for (m in mods[ageMods==a]){
  gnM   = rownames(moduleColorsL)[moduleColorsL[,a]==m]
  exprTmp = colMeans(meanExpr[gnM,])
  diffFC  = diff(sort(-exprTmp)[1:2])
  if(diffFC>fcThresh) topLay[m] = laysTmp[which.max(exprTmp)[1]]
 }
}
plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="",main="Maximum Layer (mean log2FC>0.25)")
text(1:10,rep(0,10),AGES,srt=90)
points(xMods,yMods,pch=19,cex=5,col=strucCol2)
segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
text(xMods,yMods,topLay)
dev.off()
# This is used for annotating the modules in figure 2a for prenatal laminar enrichment


###################################################################################################
###################################################################################################
print("Plot the enrichments for the laminar markers presented in Figure 4.")

laminarMarkers = read.csv(paste(inputFolder,"laminar.summary.all_L1-WM_dev.csv",sep=""),header=TRUE)
laminarMarkers = laminarMarkers[laminarMarkers$region=="NCX_Ocx_V1",]
agesLA         = c("E80","E90","E120","0M","3M","12M","48M");
layersLA       = c("L1","L2_3","L4","L5","L6","WM")  
colsLA         = c("pink","#ABC5DC","#00D0FF","blue","#0093FF","red")
names(colsLA)  = layersLA
ageLA          = rep(agesLA,length(layersLA));
layerLA        = rep(layersLA[1],length(agesLA))
for (i in 2:length(layersLA))
  layerLA      = c(layerLA, rep(layersLA[i],length(agesLA)))
ageLayerLA     = paste(ageLA,layerLA,sep="_")
pValThresh     = 30  # Only plot highly significant laminar enrichments (p<10^-30)

pdf("Figure5a_PostnatalLaminarEnrichments.pdf",height=9,width=7)
reg = "NCX_Ocx_V1"
for(lay in layersLA){
 plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="",main=lay)
 text(1:10,rep(0,10),AGES,srt=90)
 for (a in agesLA){
  kp = (laminarMarkers$region==reg)&(laminarMarkers$age==a)&(laminarMarkers$layer==lay)
  gn = as.character(laminarMarkers$gene)[kp]
  for (m in mods[ageMods==a]){
   gnM  = rownames(moduleColorsL)[moduleColorsL[,a]==m]
   val  = phyper3(length(bestGenes),length(gn),length(gnM),length(intersect(gn,gnM)))
   val  = min(val*length(agesLA)*length(layersLA),1)
   val  = -log10(val)
   if (val>pValThresh)  points(xMods[m],yMods[m],cex=5,pch=19,col=colsLA[lay])
  }
 }
 points(xMods,yMods,pch=1,cex=sqrt(countMods)/10,col="grey")
 segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
}
dev.off()


###################################################################################################
###################################################################################################
print("Plot enrichments for cell types, using the same lists as in Figure 4.")

cellGenes = read.csv(paste(inputFolder,"userInputMarkers_updatedTasic.csv",sep=""))
useCats = as.character(unique(cellGenes[,2]))

# Get the actual p-values
pvCat <- gnDis <- list()
for (cat in useCats)  gnDis[[cat]] = intersect(bestGenes,as.character(cellGenes[cellGenes[,2]==cat,1]))
for (cat in useCats){  
 pvCat[[cat]] = rep(1,length(mods));
 names(pvCat[[cat]]) = mods
 for (m in mods){
  gnIn =  bestGenes[rowSums(moduleColorsL==m)>0]
  pvCat[[cat]][m] = phyper3(length(bestGenes),length(gnIn),length(gnDis[[cat]]),length(intersect(gnIn,gnDis[[cat]])))
 }
}
pvalThresh     = 10^(-15)
useCols        = c("grey","#FF413A","grey","grey","#FE5A00")
names(useCols) = useCats

pdf("Figure5a_ModuleEnrichmentsForCellTypes.pdf",height=9,width=7)
for(cat in useCats){
 plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="",
   main=paste(substr(cat,1,30),"- Pvalue <",pvalThresh))
 text(1:10,rep(0,10),AGES,srt=90)
 m = pvCat[[cat]]<pvalThresh
 points(xMods[m],yMods[m],cex=5,pch=19,col=useCols[cat])  
 points(xMods,yMods,pch=1,cex=sqrt(countMods)/10,col=strucCol)
 segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
}
dev.off()
 


###################################################################################################
###################################################################################################
print("Code for Figure 5b ==>  Identify module enrichments for disease.")
print("Run a permutaton analysis to explicitly model emperical Pvalue for module disease gene significance.")

diseaseGenes = read.csv(paste(inputFolder,"diseaseGenesFromGeneticStudies.csv",sep=""))
useCats = as.character(unique(diseaseGenes[,2]))

# Get the actual p-values
pvCat <- gnDis <- list()
for (cat in useCats)  gnDis[[cat]] = intersect(bestGenes,as.character(diseaseGenes[diseaseGenes[,2]==cat,1]))
for (cat in useCats){  
 pvCat[[cat]] = rep(1,length(mods));
 names(pvCat[[cat]]) = mods
 for (m in mods){
  gnIn =  bestGenes[rowSums(moduleColorsL==m)>0]
  pvCat[[cat]][m] = phyper3(length(bestGenes),length(gnIn),length(gnDis[[cat]]),length(intersect(gnIn,gnDis[[cat]])))
 }
}

# Get the permuation p-values
pvDis    = rep(1,length(useCats));  names(pvDis) = useCats
numPerms = 999
pvRand   = list();  for (cat in useCats) pvRand[[cat]] = rep(1,length(numPerms))
for (p in 1:numPerms){
 set.seed(p)
 mcRand = moduleColorsL
 for (a in colnames(mcRand))  mcRand[,a] = sample(mcRand[,a],length(mcRand[,a]))
 for (cat in useCats){ 
  pvTmp = NULL 
  for (m in mods){
   gnIn =  bestGenes[rowSums(mcRand==m)>0]
   pvTmp = c(pvTmp,phyper3(length(bestGenes),length(gnIn),length(gnDis[[cat]]),length(intersect(gnIn,gnDis[[cat]]))))
  }
  pvRand[[cat]][p] = min(pvTmp)
 }
}

pdf("Figure5b_ModuleEnrichmentsForGeneticDiseaseGeneLists.pdf",height=9,width=7)
empPvalue = 0.1
for(cat in useCats){
 pvDis[cat] = sort(pvRand[[cat]])[floor((numPerms+1)*empPvalue)]
 plot(0,0,col="white",axes=FALSE,ylim=c(0,11),xlim=c(0,11),xlab="",ylab="",
   main=paste(substr(cat,1,30),"- Emp. Pvalue =",empPvalue))
 text(1:10,rep(0,10),AGES,srt=90)
 m = pvCat[[cat]]<pvDis[[cat]]
 points(xMods[m],yMods[m],cex=5,pch=19,col="red")  # Modules passing empPvalue
 points(xMods,yMods,pch=1,cex=sqrt(countMods)/10,col=strucCol)
 segments(segX1,segY1,segX2,segY2,lwd=segW,col="grey")
}
dev.off()


###################################################################################################
###################################################################################################
print("Plot the expression of all disease genes (present in at least two modules) in V1 and ACG.")

layerPositionsS2 = read.csv(paste(inputFolder,"layerRectanglePositions_small2.csv",sep=""),row.names=1)
agePositionsS    = read.csv(paste(inputFolder,"ageRectanglePositions_small.csv",sep=""),row.names=1)
layerPlot        = as.character(samplePPrep$layer)
layerPlot[layerPlot=="L2-3"] = "CPo"
rl    = as.character(samplePPrep$subregion)
names(layerPlot) <- names(rl) <- rownames(samplePPrep)
kpV1  = (samplePPrep$subregion=="V1")&(!is.element(layerPlot,c("TMZ","ICD")))
kpACG = (samplePPrep$subregion=="ACG")&(!is.element(layerPlot,c("TMZ","ICD")))
gnAll = list()
empPvalue      = 0.1
modCountThresh = 2
for (cat in useCats){
 pvDis[cat] = sort(pvRand[[cat]])[floor((numPerms+1)*empPvalue)]
 modsKp = mods[pvCat[[cat]]<pvDis[[cat]]]
 gnTmp  = NA
 for (m in modsKp)
  gnTmp = c(gnTmp,intersect(gnDis[[cat]],bestGenes[rowSums(moduleColorsL==m)>0]))
 gnAll[[cat]] = table(gnTmp)
}
for (cat in useCats) print(paste(cat,sum(gnAll[[cat]]>=1),"/",length(intersect(gnDis[[cat]],bestGenes))))

datScale = colMeans(exprPPGrep)-mean(exprPPGrep)
# Scale gene expression to the average expression in each sample across all genes

pdf("Figure5b_AverageExpressionLevelOfDiseaseGenesInNCX.pdf",width=7,height=7)
for (cat in useCats){
 disGn  = names(gnAll[[cat]])[gnAll[[cat]]>=modCountThresh]
 datGn  = colMeans(exprPPGrep[disGn,])-datScale
 ctGn   = length(disGn)
 plotMacaqueCortexSmall(datGn[kpV1],layerPlot[kpV1],samplePPrep$age[kpV1],layerPositionsS2,
    agePositionsS,plotTitle=paste(ctGn,cat,"genes in V1"),quantileScale=c(0.05,0.95),isLog2=TRUE,legendPos = c(8,-13,-10))
 plotMacaqueCortexSmall(datGn[kpACG],layerPlot[kpACG],samplePPrep$age[kpACG],layerPositionsS2,
    agePositionsS,plotTitle=paste(ctGn,cat,"genes in ACG"),quantileScale=c(0.05,0.95),isLog2=TRUE,legendPos = c(8,-13,-10))
}
dev.off()


###################################################################################################
###################################################################################################
print("Plot the expression of select ID genes in ACG")

pdf("Figure5e_expressionOfIDgenesInNCX.pdf",width=7,height=7)
idGenes = c("CPLX1","ASPM","ASXL3","RTTN")
for (gn in idGenes){
 datGn  = exprPPGrep[gn,]-datScale
 plotMacaqueCortexSmall(datGn[kpV1],layerPlot[kpV1],samplePPrep$age[kpV1],layerPositionsS2,
    agePositionsS,plotTitle=paste(gn,"in V1"),quantileScale=c(0.05,0.95),isLog2=TRUE,legendPos = c(8,-13,-10))
 plotMacaqueCortexSmall(datGn[kpACG],layerPlot[kpACG],samplePPrep$age[kpACG],layerPositionsS2,
    agePositionsS,plotTitle=paste(gn,"in ACG"),quantileScale=c(0.05,0.95),isLog2=TRUE,legendPos = c(8,-13,-10))
}
dev.off()

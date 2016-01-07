## This code block runs the WGCNA.  It is NOT a function call, but rather a separate code block from the
## Figure2_ST3.r file, which is separated for convenience.


print("Set up all of the variables and parameters in advance.")
minKMEtoStay = 0.4   # Genes with lower kME are assigned to the grey module
threshM      = 0.1   # Iteratively merge modules until the distance of all modules is greater than threshM
minModSize   = 50    # Smallest allowable module for the initial network generation
N            = 9331  # Use 75% of top variable genes in calculation.  12441 # USE ALL GENES
power        = 14
deepSplit    = 1
maxModCount  = 15    # THIS IS DONE SEPARATELY FOR EACH BLOCK
maxBlockSize = 17500 # ONLY INCLUDE 1 BLOCK
if(!exists("kpSubset"))     kpSubset     = samplePPrep$age==AGE
if(!exists("matchModules")) matchModules = NA
if(!exists("omitLabels"))   omitLabels   = NULL
if(!exists("maxModCount"))  maxModCount  = 1000
if(!exists("kpRegion"))     kpRegion     = rep(TRUE,length(kpSubset))   
if(!exists("dirName"))      dirName      = paste("NCX_WGCNA",AGE,sep="_")
if(!exists("datExprRo"))    datExprRo    = exprPPGrep[,kpSubset&kpRegion]
dir.create(dirName)


print(paste("Run automated WGCNA using the",N,"most variable genes."))
saveFn     = paste(dirName,"WGCNA_Robject.RData",sep="/")
varR0      = apply(datExprRo,1,var)
topN       = names(-sort(-varR0))[1:N]
datExprRun = t(datExprRo[topN,])
datExprRun = as.data.frame(datExprRun)
ignoreMods = paste(AGES,"M00",sep="_")

# THIS IS THE MAIN FUNCTION CALL FOR THE WGCNA
blockRun = blockwiseModules(datExprRun, checkMissingData = TRUE, maxBlockSize = maxBlockSize,
  power = power, networkType = "signed", deepSplit = deepSplit, minModuleSize = minModSize, 
  minCoreKMESize = minModSize/3, minKMEtoStay = minKMEtoStay, mergeCutHeight = threshM, 
  numericLabels = TRUE, verbose = 1)
unmergedLabels = blockRun$colors
names(unmergedLabels) <- colnames(datExprRun)


print("Merge modules if there are more than the specfied amount at any age.")
mergedLabels = unmergedLabels
modCount = length(unique(mergedLabels))
while(modCount>maxModCount){
 mergedLabels = as.numeric(mergedLabels)
 mergedMEs = moduleEigengenes(datExprRun, colors=mergedLabels, verbose=0,excludeGrey=TRUE)$eigengenes
 corMEdiss = 1-cor(mergedMEs);  diag(corMEdiss)=1
 mLab = rownames(which(corMEdiss==min(corMEdiss),arr.ind=TRUE))
 mLab = substr(mLab,3,nchar(mLab))
 mergedLabels[mergedLabels == mLab[1]] = mLab[2]
 modCount = length(unique(mergedLabels))
 print(paste("===Merged",mLab[1],"and",mLab[2],"- Distance =",signif(min(corMEdiss),3),"-",modCount,"modules remain."))
} 


print("Assign each gene to a final module based on kME and output some representative plots and data.")
colorMergedDeleted = mergedLabels
names(colorMergedDeleted) = topN
colorMergedDeleted = paste("M",colorMergedDeleted,sep="");    
a = nchar(colorMergedDeleted)==2
colorMergedDeleted[a] = paste("M0",substr(colorMergedDeleted[a],2,2),sep="")
colorMergedDeleted = paste(AGE,colorMergedDeleted,sep="_")
ignoreMods   = paste(AGES,"M00",sep="_")

MEs    = moduleEigengenes(datExprRun,colorMergedDeleted,excludeGrey=TRUE,verbose=0,grey=paste(AGE,"M00",sep="_"))$eigengenes
colnames(MEs) = substr(colnames(MEs),3,nchar(colnames(MEs)))
rownames(MEs) = colnames(datExprRo)
kME    = cor(t(datExprRo),MEs,use="p")
kME[is.na(kME)] = 0
mods   = c(colnames(kME),paste(AGE,"M00",sep="_"))
maxKme = apply(kME,1,max)
moduleColors = apply(kME,1,function(x,y) return(y[which.max(x)[1]]),mods)
moduleColors[maxKme<minKMEtoStay] = paste(AGE,"M00",sep="_")
moduleColorsOld1  = moduleColors


print("Relabel modules and order based on the percent neuronal content of module (Cahoy et al 2008).")
data(BrainLists)
bl = c("Astrocyte_probable__Cahoy","Neuron_probable__Cahoy","Oligodendrocyte_probable__Cahoy")
comparisons = matrix(FALSE,nrow=length(colorMergedDeleted),ncol=length(bl))
rownames(comparisons) = names(unmergedLabels);         colnames(comparisons) = bl
for (b in bl)
  comparisons[is.element(rownames(comparisons),BrainLists[BrainLists[,2]==b,1]),b] = TRUE
colorAssignedF = factor(colorMergedDeleted,levels = names(table(colorMergedDeleted)))
names(colorAssignedF) = names(unmergedLabels)
out = table(colorAssignedF)
for (b in bl) out = cbind(out, table(colorAssignedF[comparisons[,b]]))
colnames(out) = c("TotalGenes",colnames(comparisons))

out2 = out[,c(1:4,1)]
colnames(out2) = c("TotalGenes","Astrocyte","Neuron","Oligodendrocyte","OtherGenes")
out2 = apply(out2,1,function(x) { x[5] = x[1]-(x[2]+x[3]+x[4]);  x = 100*x/x[1];  x[5] = x[5];  return(x)})
out2 = -(out2[3,]-out2[2,]-out2[4,])
out2 = out2[setdiff(names(out2),ignoreMods)]
ord      = order(out2)
labels   = c(colnames(MEs),paste(AGE,"M00",sep="_"))
names(labels) = c(colnames(MEs)[ord],paste(AGE,"M00",sep="_"))


if(!is.na(matchModules[1])){
 print("Reassign module colors to roughly match module ordering from adjacent later age.")
 ovTable = overlapTable(moduleColors,matchModules,ignore=ignoreMods)
 countTable = ovTable$countTable;  
 pTable  = -log10(ovTable$pTable+10^(-299))+countTable/1000;
 pTable[grep("_M",rownames(pTable)),grep("_C",colnames(pTable))] = 0
 pTable[grep("_C",rownames(pTable)),grep("_M",colnames(pTable))] = 0
 ord     = order(apply(pTable,1,which.max))
 seqs    = paste("M",c(paste("0",0:9,sep=""),as.character(10:99)),sep="")
 labels  = seqs[1:(length(rownames(pTable))+1)]
 labels  = paste(AGE,labels,sep="_")
 names(labels)      = c(paste(AGE,"M00",sep="_"),rownames(pTable)[ord])
 colorMergedDeleted = labels[colorMergedDeleted]
 moduleColorsOld2   = moduleColors
 moduleColors       = labels[moduleColors]
 moduleColorsOld3   = moduleColors
 
 ovTable = overlapTable(moduleColors,matchModules,ignore=ignoreMods)
 pTable  = -log10(ovTable$pTable+10^(-299));    
 pTable2 = pmin(pTable,100);     countTable = ovTable$countTable;   
 xLabels = colnames(pTable);     yLabels    = rownames(pTable)
 
 pdf(paste(dirName,"moduleComparisonWithPreviousAge.pdf",sep="/"), height=9,width=15)  
 par(mfrow=c(1,1));   par(cex = 1.0);   par(mar=c(8, 12.4, 2.7, 1)+0.3);
 labeledHeatmap(Matrix = pTable2, xLabels = paste(" ",xLabels), yLabels = paste(" ",yLabels), 
  colorLabels = TRUE, ySymbols = yLabels, main = " ", xSymbols = xLabels, textMatrix = countTable, 
  colors = greenWhiteRed(100)[50:100], setStdMargins = FALSE, cex.text = 1.0, cex.lab = 1.0);
 dev.off()
} else {
 colorMergedDeleted = labels[colorMergedDeleted]
}


print("Finalize and save the modules.")
MEs    = moduleEigengenes(datExprRun,colorMergedDeleted,excludeGrey=TRUE,verbose=0,grey=paste(AGE,"M00",sep="_"))$eigengenes
colnames(MEs) = substr(colnames(MEs),3,nchar(colnames(MEs)))
rownames(MEs) = colnames(datExprRo)
kME    = cor(t(datExprRo),MEs,use="p")
kME[is.na(kME)] = 0
mods   = colnames(kME)
maxKme = apply(kME,1,max)
moduleColors = apply(kME,1,function(x,y) return(y[which.max(x)[1]]),mods)
moduleColors[maxKme<minKMEtoStay] = paste(AGE,"M00",sep="_")
kMEtable = data.frame(gene=rownames(kME),kME=maxKme,module=moduleColors)
kMEtable = kMEtable[order(factor(moduleColors,levels=mods),-maxKme),]
kME      = kME[rownames(kMEtable),mods]

write.csv(kMEtable,paste(dirName,"WGCNA_MEcorrelationTable_kME.csv",sep="/"),row.names=FALSE)
save(MEs,mergedLabels,colorMergedDeleted,blockRun,moduleColors,file=saveFn)

pdf(paste(dirName,"WGCNA_dendrogram.pdf",sep="/"),height=8,width=15)
mName  = c("Intial Modules","Final Assignments")
cs     = c("grey",standardColors())
mColor = cbind(labels2colors(mergedLabels),labels2colors(moduleColors[topN],colorSeq=cs))
plotDendroAndColors(blockRun$dendrograms[[1]],mColor,mName,dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

# Remove variables that need to be reset
rm(kpSubset,datExprRo)

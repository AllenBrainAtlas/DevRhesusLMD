# moduleSwap needs to be a vector where the VALUES are the DESIRED module names and names(moduleSwap) are the CURRENT module names.
# This code block is not a function, but rather a section of code cut out of the main Figure2_ST3.r for convenience.  Some parts
#   of the code are redundant with the autoamted_WGCNA_code_Block.r and may slow things down, but are included so that files
#   that have already been written out are overwritten with the correct module labels.

print("Read in the orginal network information and the reassigned module labels.")
saveFn  = paste(dirName,"WGCNA_Robject.RData",sep="/")
load(saveFn)

cn = names(moduleSwap)
moduleSwap = c(moduleSwap,paste(AGE,"M00",sep="_"))
names(moduleSwap) = c(cn,paste(AGE,"M00",sep="_"))
colorMergedDeleted_old = colorMergedDeleted
colorMergedDeleted = moduleSwap[as.character(colorMergedDeleted)]

N            = 9331  # Use 75% of top variable genes in calculation.  12441 # USE ALL GENES
minKMEtoStay = 0.4   # Genes with lower kME are assigned to the grey module
varR0        = apply(datExprRo,1,var)
topN         = names(-sort(-varR0))[1:N]
datExprRun   = t(datExprRo[topN,])
datExprRun   = as.data.frame(datExprRun)
ignoreMods   = paste(AGES,"M00",sep="_")


print("Reassign genes to modules - this should have no effect, except in networks with merged modules.")
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
save(MEs,mergedLabels,colorMergedDeleted,colorMergedDeleted_old,blockRun,moduleColors,file=saveFn)


if(!is.na(matchModules[1])){
 library(gplots)
 ovTable    = overlapTable(moduleColors,matchModules,ignore=ignoreMods)
 pTable     = -log10(ovTable$pTable+10^(-299));    
 pTable2 = pmin(pTable,100);     countTable = ovTable$countTable;   
 xLabels = colnames(pTable);     yLabels    = rownames(pTable)
 
 pdf(paste(dirName,"moduleComparisonWithPreviousAge.pdf",sep="/"), height=9,width=15)  
 par(mfrow=c(1,1));   par(cex = 1.0);   par(mar=c(8, 12.4, 2.7, 1)+0.3);
 labeledHeatmap(Matrix = pTable2, xLabels = paste(" ",xLabels), yLabels = paste(" ",yLabels), 
  colorLabels = TRUE, ySymbols = yLabels, main = " ", xSymbols = xLabels, textMatrix = countTable, 
  colors = greenWhiteRed(100)[50:100], setStdMargins = FALSE, cex.text = 1.0, cex.lab = 1.0);
 dev.off()
} 

rm(kpSubset,datExprRo)

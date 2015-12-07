## Code for generating figures 1d, 3e, 4a,b,e, and 5a-f in the NHP paper

##############################################################################################
##############################################################################################
print("Read in and properly format the data.")

load("\\\\aibsdata/humancelltypes/manuscripts/DevMacaqueLMD/Analysis_code/TB/cache/nhp_PrePost_StartingData.RData")
inputFolder  = "\\\\aibsdata/humancelltypes/JeremyM/rhesusMonkey_prenatal_LCM/ScienceAnalyses/inputFiles/"
outputFolder = "\\\\aibsdata/humancelltypes/JeremyM/rhesusMonkey_prenatal_LCM/ScienceAnalyses/outputFiles/"
dir.create(outputFolder)
setwd(outputFolder)

# LIBRARIES
library(WGCNA);
library(reshape);
library(gplots)

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

medianVals = 2^apply(exprPP,1,median)  # THESE ARE IN LINEAR SPACE  # I don't remember why I need this value.... XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#sa=scan("\\\\Aibsdata/humancelltypes/Analyses/DevMacaqueLMD/data/struc_acro.txt",what="character")  # ALREADY INCLUDED IN SAMPLEINFO
#samplePP$struc_acro = sa


##############################################################################################
##############################################################################################
print("Code for Figure 1d.")
print("-----Note: Only part of each of these plots were used in making figure 1d.")


source(paste(inputFolder,"functionsForMacaqueAnalysis.r",sep=""))
source(paste(inputFolder,"plotMacaqueCortex.r",sep=""))
layerPositions  = read.csv(paste(inputFolder,"layerRectanglePositions.csv",sep=""),row.names=1)
regionPositions = read.csv(paste(inputFolder,"regionRectanglePositions.csv",sep=""),row.names=1)
ageOffsets      = read.csv(paste(inputFolder,"ageRectangleOffsets.csv",sep=""),row.names=1)

fig1dGenes = c("PAX6","EOMES","DCX","SYT1","GFAP","MOBP")
fig1dGenes = which(is.element(bestGenes,fig1dGenes))
dir.create(paste(outputFolder,"fig_1e/",sep=""))
for (i in fig1dGenes) {
 g  = bestGenes[i];
 hg = humanGenes[i];   
 if (!is.na(hg)) g = paste(g,hg,sep="__")
 p  = bestProbes[i];
 e  = bestEntrez[i];
 fn = paste(outputFolder,"fig_1e/",g,"__",e,"__",p,".pdf",sep="")
 pdf(fn,width=18,height=9)
 plotMacaqueCortex(exprPP[p,],NULL,samplePP$layer,samplePP$subregion,samplePP$age,
    layerPositions,regionPositions,ageOffsets,paste(g,e,p,sep=" - "),
    quantileScale=c(0.1,0.95),medianVals=medianVals)
 dev.off()
}


##############################################################################################
##############################################################################################
print("Code for Figure 3e.")

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

gns    = c("MOG","MOBP","ERMN","MAL","ASPA","OPALIN")
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

pdf(paste(outputFolder,"Oligodendrocyte_select_genes_with_time_Fig3e.pdf",sep=""),height=6,width=6);
for(gn in gns){      
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
print("Code for Figure 4a, b, and e.")

# Trygve's pre-generated laminar marker genes.  I do not have the code for making these gene lists...

laminarMarkers      = read.csv(paste(inputFolder,"laminar.summary.all_L1-WM_dev.csv",sep=""),header=TRUE)
laminarMarkers      = laminarMarkers[laminarMarkers$region=="NCX_Ocx_V1",]
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


# For each layer, find the top 10% (1244) laminar genes at each time point.
N=1244 # 
laminarGenesLA = list();
for (i in 1:length(ageLayerLA)) {
  l = layerLA[i];  
  a = ageLA[i];  
  laminarGenesLA[[ageLayerLA[i]]] = as.character(laminarMarkers$gene)[(laminarMarkers$age==a)&(laminarMarkers$layer==l)]
}

percentOverlap = matrix(0,nrow=length(layersLA),ncol=length(agesLA));
rownames(percentOverlap) <- layersLA;  colnames(percentOverlap) <- agesLA
countSignificant <- percentSignificant <- percentOverlap
for (l in layersLA) for (a in agesLA){   nm = paste(a,l,sep="_")
 percentOverlap[l,a]     = 100*length(intersect(orderedLaminarGenes[1:N,nm], 
                           orderedLaminarGenes[1:N,paste("48M",l,sep="_")]))/N
 percentSignificant[l,a] = 100*length(intersect(laminarGenesLA[[nm]], 
                           orderedLaminarGenes[1:N,paste("48M",l,sep="_")]))/N
 countSignificant[l,a]   = length(laminarGenesLA[[nm]])
}
percentOverlap     = percentOverlap[,agesLA];
percentSignificant = percentSignificant[,agesLA];    

print("-----Make the plot for figure 4a")
xs        = 8.4-0.33*(1:10)
layersLA2 = c("L2_3","L4","L5","L6"); 
agesLA2   = agesLA[1:6]
x         = 1:6;
cols      = c("purple","cyan","green","red")
ymax      = max(percentOverlap[layersLA2,agesLA2],na.rm=TRUE)+7

pdf(paste(outputFolder,"topLaminarGenes_newLists_linePlots_Fig4a.pdf",sep=""),height=6,width=6);
plot(x,colMeans(percentOverlap[layersLA2,agesLA2]),main=paste("% of top",N,"laminar genes in common with 48M"),
  type="l",lwd=5,xlab="Age (log2[pcd])",ylab="% in common with 48M",ylim=c(-10,ymax)); 
abline(h=0);
text(x,rep(-5,length(x)),agesLA2)
for (i in 1:length(layersLA2))
 lines(x,percentOverlap[layersLA2,agesLA2][i,],type="l",lwd=4,col=cols[i])
for (i in 1:length(x)) segments(x[i],0,x[i],-1)
abline(h=10,lty="dashed",lwd=1.5);
abline(v=3.5,lwd=1.5)
dev.off()


# Now find which (if any) genes are significant at ALL time points  # CODE NOT USED???
#laminarGenesAllL = NULL;  for (l in layersLA) laminarGenesAllL[[l]] = bestGenes
#for (i in 1:length(ageLayerLA))   laminarGenesAllL[[layerLA[i]]] = 
#   sort(intersect(laminarGenesAllL[[layerLA[i]]], orderedLaminarGenes[1:N,i]))

print("-----Make the plot for figure 4b")

ages2inv <- function(ages){
 convert = c(40,50,70,80,90,120,170,255,530,1625,2000)
 names(convert) = c("E40","E50","E70","E80","E90","E120","0M","3M","12M","48M","Adult")
 return(1-1/convert[ages])
}
  
kpV1  = (samplePP$subregion=="V1")&(samplePP$age_log2pcd>6)&(samplePP$layer!="CP")&
         is.element(samplePP$layer_dev,c("L2","L3","L4","L5","L6"))
whMax = function(x,nm) return(nm[which.max(x)]) 

meanV1   = list()
layerTmp = as.character(samplePP$layer_dev)
layerTmp[is.element(layerTmp,c("L2","L3"))] = "L2_3"
for (a in agesLA){       
 kpV            = kpV1&(samplePP$age==a);
 meanV1[[a]]    = findFromGroups2(t(exprPPG[,kpV]),layerTmp[kpV])
}

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

pdf("LaminarSimilarityBetweenAges_top10perc_dots_Fig4b.pdf",width=6.5,height=7);
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


print("-----Make the plot for figure 4e")
print("-----Note: Only part of each of these plots were used in making figure 4e.")

layerPositionsS = read.csv(paste(inputFolder,"layerRectanglePositions_small.csv",sep=""),row.names=1)
agePositionsS   = read.csv(paste(inputFolder,"ageRectanglePositions_small.csv",sep=""),row.names=1)

l23Genes   = c("CNN3","DNAH5","PRR15","MLF1")
l4Genes    = c("POLQ","CALB1","CPA6","RORA")
l5Genes    = c("FAM84A","LRRK1","EMILIN1","GRM8")
l6Genes    = c("SLC26A7","LGR5","CA12","SYT6")
fig4eGenes = c(l23Genes,l4Genes,l5Genes,l6Genes)
fig4eGenes = which(is.element(bestGenes,fig4eGenes))
fd = paste(outputFolder,"fig_4e/",sep="")
dir.create(fd)

kpV1       = samplePP$subregion=="V1"
agePlot    = as.character(samplePP$age)[kpV1]
layerPlot  = as.character(samplePP$layer)[kpV1];     layerPlot[layerPlot=="L2-3"] = "CPo"
layerPlot[is.element(layerPlot,c("L4Cb","L4Ca","L4B","L4A"))] = "L4"
kpACG      = samplePP$subregion=="ACG"
agePlotA   = as.character(samplePP$age)[kpACG]
layerPlotA = as.character(samplePP$layer)[kpACG];     layerPlotA[layerPlotA=="L2-3"] = "CPo"
layerPlotA[is.element(layerPlotA,c("L4Cb","L4Ca","L4B","L4A"))] = "L4"
legendPos  = c(8,-11.5,-8.5);

for (i in fig4eGenes) {
 g  = bestGenes[i];
 hg = humanGenes[i];
 if (!is.na(hg)) g = paste(g,hg,sep="__")
 p  = bestProbes[i];
 e  = bestEntrez[i];
 fn = paste(fd,g,"__",e,"__",p,".pdf",sep="")
 
 pdf(fn,width=7,height=7)
 plotMacaqueCortexSmall(exprPP[p,kpV1],layerPlot,agePlot,layerPositionsS,agePositionsS,
  paste(g,e,p,"V1",sep=" - "),isLog2=TRUE,combineFn="meanNA",quantileScale=c(0.1,0.95),
  linearOrLog="linear",bgPar="grey95",displayLayers=FALSE,legendPos=legendPos)
 abline(v=c(0,6,10),lwd=2);   abline(h=c(0,-8,-12),lwd=2)
 dev.off();
}


##############################################################################################
##############################################################################################
print("Code for Figure 5.")
print("-----This combines replicate samples as in figure 3e to avoid artificially improving statistical power.")
print("-----All plots in this figure have V1 as some shade of blue and ACG as some shade of red.")

# Prepare the data for this figure
layKp3   = c("L2","L3","L2-3","L5","L6","CPi","CPo","VZ","VZi","VZo","SZ","SZi","SZo")
kpSamp3  = is.element(samplePPrep$subregion,c("ACG","V1"))&is.element(samplePPrep$layer,layKp3)
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
age3av   = age3[match(laAgDons,laAgDon3)]
layAll3av= layAll3[match(laAgDons,laAgDon3)]
layCol3av= layCol3[match(laAgDons,laAgDon3)]
lays3    = c("L2/L3","L5/L6","SZ","VZ")

# Run all of the t-tests (note that the p-values are uncorrected for multiple comparisons)
t.test.l.paired2 <- function(x){  
  l  = length(x)/2
  tt = t.test(x[1:l],x[(l+1):(2*l)],paired=TRUE)
  return(tt$p.val)
}

fc3 <- pval3 <- list()
for (l in lays3){
 pval3[[l]] = matrix(NA,nrow=dim(exprG3)[1],ncol=length(ages))
 colnames(pval3[[l]]) = ages;  rownames(pval3[[l]]) = rownames(exprG3)
 fc3[[l]]   = pval3[[l]]
 for (a in ages){  kp = (age3av==a)&(layCol3av==l)
  if (sum(kp)>=3){
   pval3[[l]][,a] = apply(cbind(exprA3[,kp],exprV3[,kp]),1, t.test.l.paired2)
   fc3[[l]][,a]   = rowMeans(exprA3[,kp])-rowMeans(exprV3[,kp])
  }
 }
}


print("-----Make the plot for figure 5a")

pdf(paste(outputFolder,"NumberOfDexGenes_perAge_fig5a.pdf",sep=""),height=15,width=7);
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


print("-----Make a supplemental table with all of these significant genes")

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
write.csv(regOutZZ,"SupplementalTable_V1vsACG_markerGenes.csv",row.names=FALSE)


print("-----Perform (and plot) enrichment analyses for figure 5b-e")
print("----------Note: code for generated the line plots is elsewhere.")

enrichs3V <- enrichs3A <- NULL
genes3   = rownames(pval3[[l]])
fnIn     = paste(inputFolder,"userInputMarkers.csv",sep="")  
  
for (l in lays3) for (a in ages){
 lala   = paste(l,a)
 label3 = rep("grey",length(genes3))
 label3[(pval3[[l]][,a]<pThresh)&(fc3[[l]][,a]< -fcThresh)] = lala
 if(sum(label3==lala)>=5){
  enrichsTmp = userListEnrichment(genes3,label3,catNmIn="input",fnIn=fnIn)
  enrichsTmp = enrichsTmp$sigOverlaps[,c(1,2,4)]
  enrichs3V  = rbind(enrichs3V,enrichsTmp)
 }
}
for (l in lays3) for (a in ages){
 lala   = paste(l,a)
 label3 = rep("grey",length(genes3))
 label3[(pval3[[l]][,a]<pThresh)&(fc3[[l]][,a]>fcThresh)] = lala
 if(sum(label3==lala)>=5){
  enrichsTmp = userListEnrichment(genes3,label3,catNmIn="input",fnIn=fnIn)
  enrichsTmp = enrichsTmp$sigOverlaps[,c(1,2,4)]
  enrichs3A  = rbind(enrichs3A,enrichsTmp)
 }
}
omitX = c(1,2,3,1,7,8,9,10,7,8,9,10)
omitY = c(3,3,3,2,1,1,1,1,0,0,0,0)

# Scale results for appropriate multiple comparisons
pScale   = 242/4  # ALL p-values are scaled by this value to adjust for 238 additional tests not used in this paper or included in the code
enrichs3A[,3] = enrichs3A[,3]*pScale
enrichs3A = enrichs3A[enrichs3A[,3]<pThresh,]
enrichs3V[,3] = enrichs3V[,3]*pScale
enrichs3V = enrichs3V[enrichs3V[,3]<pThresh,]


print("-----Plot enrichment analyses for figure 5b-e")
print("----------Note that the legend was generated by comparing the size of the boxes")
print("----------  in these plots to the p-values output from userListEnrichment.")

pdf(paste(outputFolder,"EnrichmentPlots_ACGvsV1_Fig5bcde.pdf",sep=""),width=6,height=3.5)
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


print("-----Plot bar plot for fig 5f (page 1) and the inset (page 2)")
print("----------Note that the gene line plots are plotted separately.")

countV = matrix(0,ncol=length(lays3),nrow=length(rownames(pval3[[l]])));     
colnames(countV) = lays3;
rownames(countV) = rownames(pval3[[l]]);
countA <- countV

for (l in lays3) for (a in ages){
 label3 = genes3[(pval3[[l]][,a]<pThresh)&(fc3[[l]][,a]< -fcThresh)]
 if (sum(is.na(label3))==0)  countV[label3,l] = countV[label3,l]+1
 label3 = genes3[(pval3[[l]][,a]<pThresh)&(fc3[[l]][,a]>fcThresh)]
 if (sum(is.na(label3))==0)  countA[label3,l] = countA[label3,l]+1
}
topV = countV;  topA = countA
topV[,1] = floor(topV[,1]/5);
topV[,2] = floor(topV[,2]/6);
topV[,3] = floor(topV[,3]/4);
topV[,4] = floor(topV[,4]/4)
topA[,1] = floor(topA[,1]/5);
topA[,2] = floor(topA[,2]/6);
topA[,3] = floor(topA[,3]/4);
topA[,4] = floor(topA[,4]/4)
countsAV = colSums(cbind(countV>0,countA>0))[c(1,5,2,6,3,7,4,8)]
topsAV   = cbind(topV,topA)[,c(1,5,2,6,3,7,4,8)]
perAV    = (100*colSums(topsAV)/countsAV)
perAV    = rbind(perAV,100-perAV)
nm       = c("L23_V1","L23_ACG","L56_V1","L56_ACG","SZ_V1","SZ_ACG","VZ_V1","VZ_ACG")
colnames(perAV) = paste(colnames(perAV),c("V1","ACG"),sep=",")

pdf(paste(outputFolder,"barplotOfConsistentMarkers_Fig5f.pdf",sep=""),height=7.5,width=9.5)
barplot(perAV,ylab="# of marker genes",col="white",border="white",ylim=c(0,45))
x = -0.5+(1:8)*1.2;  y = c(1:100);  cols = c("blue","red","blue","red","blue","red","blue","red")
for (i in 1:dim(perAV)[2]){
 txt = sort(rownames(pval3[[1]])[topsAV[,i]>0],decreasing=TRUE)
 txt = setdiff(txt,txt[grep("LOC",substr(txt,1,3))])
 l   = length(txt)
 text(rep(x[i],l),y[1:l],txt,cex=1,col=cols[i])
}
abline(h=0)
par(mfrow=c(2,2))
a = length(bestGenes);  b = sum(rowSums(countV+countA)>0);  d = sum(rowSums(topV+topA)>0)
valAV2 = cbind(a-b-d, b-d, d);   names(valAV2) = c("Non-regional","Briefly regional","Persistent regional")
pie(valAV2,col=c("grey","purple","black"))
dev.off()










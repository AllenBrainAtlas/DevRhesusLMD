convertAge <- function(ageIn,speciesIn="human",speciesOut="macaque"){
 # speciesOut can also be "eventScore"
 if(speciesIn=="human")     { I1 = 3.167; I2 = 3.72;  }
 if(speciesIn=="macaque")   { I1 = 3.27;  I2 = 2.413; }
 if(speciesIn=="mouse")     { I1 = 2.145; I2 = 1.894; }
 eventScore = (log(2^ageIn)-I1)/I2
 if(speciesOut=="eventScore") return(eventScore)
 if(speciesOut=="human")    { R1 = 3.167; R2 = 3.72;  }
 if(speciesOut=="macaque")  { R1 = 3.27;  R2 = 2.413; }
 if(speciesOut=="mouse")    { R1 = 2.145; R2 = 1.894; }
 ageOut = log2(exp(eventScore*R2+R1))
 return(ageOut)
}

t.test.l.paired2 <- function(x){  
  l  = length(x)/2
  tt = t.test(x[1:l],x[(l+1):(2*l)],paired=TRUE)
  return(tt$p.val)
}

ages2inv <- function(ages){
 convert = c(40,50,70,80,90,120,170,255,530,1625,2000)
 names(convert) = c("E40","E50","E70","E80","E90","E120","0M","3M","12M","48M","Adult")
 return(1-1/convert[ages])
}
  
ages2pcd <- function(ages,log2=TRUE){   convert = c(40,50,70,80,90,120,170,255,530,1625,2000)
 names(convert) = c("E40","E50","E70","E80","E90","E120","0M","3M","12M","48M","Adult")
 if (log2) return(log2(convert[ages]));   return(convert[ages])  }

ovTest <- function(grp1,grp2,all,pval=TRUE,writeOv=FALSE,twoTailed=TRUE){
  all  = unique(all);    grp1 = intersect(grp1,all);    grp2 = intersect(grp2,all)
  lenA = length(all);    len1 = length(grp1);           len2 = length(grp2);
  ov   = intersect(grp1,grp2);                          lenO = length(ov);
  OE   = lenO/(len1*len2/lenA)
  if(writeOv) write(ov,"")
  if(twoTailed&pval) if(OE<1) return(-phyper3(lenA,len1,len2,lenO,lt=FALSE))
  if(pval) return(phyper3(lenA,len1,len2,lenO))
  return(OE)
}

phyper3 <- function(total, group1, group2, overlap, ...){  if(overlap<=2) return(1); 
 return(max(min(1,phyper2(total, group1, group2, overlap, ...)),10^-310)) }

integer.base.b <- function(x, b=2){
 xi <- as.integer(x);  N <- length(x);  xMax <- max(x)
 ndigits <- (floor(logb(xMax, base=2))+1)
 Base.b  <- array(NA, dim=c(N, ndigits))
 for(i in 1:ndigits){	Base.b[, ndigits-i+1] <- (x %% b)
   				x <- (x %/% b)  }
 if(N==1) Base.b[1, ] else Base.b
}


getTvalue  <- function(x) {    # T TEST FOR USE WITH APPLY
  x1 = x[var[[1]]];
  x2 = x[var[[2]]]
  return((mean(x1)-mean(x2))/sqrt(var(x1)/length(x1)+var(x2)/length(x2)))
}


getAnovaPvalforApply <- function(x){
  if (!exists("varWeights")) varWeights=NULL
  anovadat  = as.data.frame(cbind(varLabels,x))
  aov.out   = summary(aov(as.numeric(anovadat[,2])~anovadat[,1],data=anovadat,weights=varWeights))
  return(aov.out[[1]]$'Pr(>F)'[1])
}


getAnovaFstatforApply <- function(x){
  if (!exists("varWeights")) varWeights=NULL
  anovadat  = as.data.frame(cbind(varLabels,x))
  aov.out   = summary(aov(as.numeric(anovadat[,2])~anovadat[,1],data=anovadat,weights=varWeights))
  return(aov.out[[1]]$'F value'[1])
}


decimals <- function(x,dec=2){
  scale = 10^dec
  return(round(x*scale)/scale)
}


memory.usage <- function(o, num=10){
        sizeO = NULL
        for (i in 1:length(o))
        sizeO = c(sizeO, object.size(get(o[i])))
        names(sizeO) = o
        sizeO = sort(sizeO,decreasing=TRUE)
        return(sizeO[1:num])
}


phyper2 <- function (total, group1, group2, overlap, verySig=TRUE ,lt=TRUE){
  # This function is the same is phyper, just allows for more sensible input values
  q = overlap
  m = group1
  n = total-group1
  k = group2
  prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
  if (verySig) return(-prob)
  return(1-prob)
}


mouse2human2 <- function (mouse, m2h){
 # Convert mouse to human symbols
 rownames(m2h) = m2h$Mou
 noHumn = which(!(mouse%in%m2h$Mouse_Symbol))
 humn = which((mouse%in%m2h$Mouse_Symbol))
 mouse[humn] = as.character(m2h[mouse[humn],1])
 mouse[noHumn] = toupper(mouse[noHumn])
 return(mouse)
}


separabilityScore <- function (x1, x2, group){
# Returns the mean within-group distance scaled by the mean between-group distance.  Scores close to ZERO are the most separable. NOTE: requires library(fields)
 x = cbind(x1,x2); 
 d = as.matrix(dist(x))
 groups = unique(group)
 within <- midpt <- NULL
 for (g in groups){
  dg     = d[group==g,group==g]
  within = c(within,dg[lower.tri(dg)])
  midpt  = rbind(midpt,c(mean(x1[group==g]),mean(x2[group==g])))
 }
 within  = mean(within)
 between = as.matrix(dist(midpt))
 between = mean(between[lower.tri(between)]) 
 return(within/between)
}


t.test.l.paired <- function(x){
  tt = t.test(x[var[[1]]],x[var[[2]]],paired=TRUE)
  return(c(tt$est,tt$stat,tt$p.val))
}


t.test.l = function(x){
 tt = t.test(x[var[[1]]],x[var[[2]]])
 return(c(tt$est,sd(x[var[[1]]]),sd(x[var[[2]]]),tt$p.val))
}


cor.test.l = function(x){
 ct = cor.test(x,var)
 return(c(ct$est,ct$p.val))
}


plotSamples <- function(data,value,GorP,phenoOrder=1:dim(data)[2],filter=1:dim(data)[2],
  filename="plotSample.pdf",las=2,xlab="",ylab="Expression",col=NA, mainIn=NA, maxToPlot=8,
  cex.main=4,...){
  data      = data[,filter]
  if(!is.na(col[1])) col = col[filter]
  if(!is.na(phenoOrder[1])){
    phenoOrder = phenoOrder[filter]
    sampOrder  = order(phenoOrder,colnames(data))
    data       = data[,sampOrder]
    if(!is.na(col[1])) col = col[sampOrder]
  }
  plotData  = data[GorP==value,]
  if(is.na(col[1])) col = "grey"

  if (length(plotData)==0){ write("Invalid gene/probe.",""); return(NULL)}
  if (length(plotData)==dim(data)[2]){ 
    plotData2 = list(as.numeric(plotData))
    names(plotData2) = value
  } else {
    if (dim(plotData)[1]>maxToPlot) plotData = plotData[1:maxToPlot,]
    plotData2 = list()
    for (i in 1:dim(plotData)[1]) plotData2[[i]] = as.numeric(plotData[i,])
    names(plotData2) = paste(value,"-",rownames(plotData))
  }
  numValues = length(plotData2)
  sampNames = colnames(data)
  if (!is.na(filename))
    pdf(filename,height=(numValues+2)*6,width=round(10+dim(data)[2]/10))
  par(mfrow=c(numValues*2,1))
  for (i in 1:numValues){
    if (is.na(mainIn)) {
     main = names(plotData2)[i]
    } else { main = mainIn }
    q = round(quantile(plotData2[[i]],c(0,.5,1)))
    main = paste(main,":  min=",q[1],"; median=",q[2],"; max=",q[3],sep="")
    barplot(plotData2[[i]],main=main,names.arg=sampNames,las=las,
      xlab=xlab,ylab=ylab,col=col,cex.main=cex.main,...)
    plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
 }
  if (!is.na(filename)) dev.off()
}


probe2Gene <- function(probes, probeList, geneList) {
 probes<-sub("X", "", probes)
 names(geneList) = probeList
 return(geneList[probes])
}


write.geneList = function(PG, filename, allProbes=0, allGenes=0, probe="g") {
 gene = PG
 if (probe!="g") {
  gene = probe2Gene(PG,allProbes,allGenes)
 }
 write(gene,filename,sep="\n")
}


findFromGroups <- function(datExpr,groupVector,fn="mean"){
# Performs a function at the group level (default is "mean") and returns a
# matrix where the output ROWS are the same as the COLUMNS of datExpr (typically
# genes or probes) and the output columns are the components of the group (i.e.,
# regions of the brain). datExpr is expression data (genes = COLUMNS, 
# samples = ROWS).  groupVector is a vector or factor corresponding to group 
# with one element per sample
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[groupVector==groups[i],]
    if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}

findFromGroups2 <- function(datExpr,groupVector,fn="mean",...){
# Performs a function at the group level (default is "mean") and returns a
# matrix where the output ROWS are the same as the COLUMNS of datExpr (typically
# genes or probes) and the output columns are the components of the group (i.e.,
# regions of the brain). datExpr is expression data (genes = COLUMNS, 
# samples = ROWS).  groupVector is a vector or factor corresponding to group 
# with one element per sample
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[groupVector==groups[i],]
    if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn,...)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}



findFromNonGroups <- function(datExpr,groupVector,largeGroupVector=
  rep("all",length(groupVector)),fn="mean"){
# Same as findFromGroups, except that the function is applied to all samples
# NOT in the relevant group.  If largeGroupVector is specified, it will apply
# the function to samples both NOT in the group but IN the largeGroup.
  groups  <- nonGroups <- names(table(groupVector))
  fn       = match.fun(fn)
  for (i in 1:length(groups))
    nonGroups[i] = largeGroupVector[which(groupVector==groups[i])[1]]
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[(groupVector!=groups[i])&(largeGroupVector==nonGroups[i]),]
    if (length(datIn)>0) if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}


remapModuleColors <- function(colorh,colOut=standardColors()) {
  ## This function remaps one set of colors to another.
  ## Note: colOut returns colors from greatest to fewest genes
  ##   in the module, excluding grey
  ## Requires library(moduleColor)
  colIn = names(sort(table(colorh), decreasing=TRUE))
  colIn = c("grey",colIn[colIn!="grey"])
  colOut = c("grey",colOut[colOut!="grey"])
  colorhOut = colorh
  for (i in 1:length(colIn)) colorhOut[colorh==colIn[i]]=colOut[i]
  return(colorhOut)  
}


visantPrepOverall <- function(colorFinal2, moduleColor, datExprrest, genes, numInt, power, signed=FALSE)
{
## This file returns the overall strongest connections within a module for one group of subjects

## USER inputs
# colorFinal2 = vector with module association for each probe
# moduleColor = color of the module... should be a member of colorFinal2
# datExprrest = expression data for the genes corresponding to cIndex
# genes = list of genes that correspond to the probes
# numInt = number of interactions to output to the visant file

cIndex = (colorFinal2==moduleColor)
datExpr=datExprrest[,cIndex]
if (signed){AdjMat =((1+cor(datExpr, use="p"))/2)^power
        } else {AdjMat =abs(cor(datExpr, use="p"))^power}
diag(AdjMat)=0
Degree <- apply(AdjMat,1,sum)
Degree = Degree/max(Degree)
meanExpr=apply(datExpr,2,mean)
ProbeSet=colnames(datExpr)
GeneSet=genes[cIndex];
Module=rep(moduleColor,length(meanExpr))

## This file summarizes intramodular connectivity and expression for each gene in each group:
fn=paste(moduleColor,"_connectivityOverall.csv",sep="")
datConn = cbind(ProbeSet,GeneSet,meanExpr,Degree,Module)
datConn = datConn[order(Degree, decreasing=TRUE),]
write.table(datConn,file=fn,sep=",",
            row.names=F, col.names= c("probes","genes","meanExpr","kin","module"))
write(paste(fn, "written."),"")

## TOM matrix
distTOM <- TOMdist1(AdjMat)
simTOM = 1-distTOM
diag(simTOM)=0
simTOMcutoff = simTOM

## Correlation matrix
Pearson = cor(datExpr ,use="p")
diag(Pearson)=0

## Dynamically determine the appropriate cutoff
cutoff = 0.24
len    = 10000
dir    = "increase"
loops  = 0
split  = 0.01
numInt = numInt+100
while(len>100){
 loops = loops+1
 if (dir == "increase") { cutoff = cutoff+split; }
 if (dir == "decrease") { cutoff = cutoff-split; }
 indices = (simTOMcutoff>cutoff)
 len = sum(sum(indices))
 if (len < numInt) {dir = "decrease";}
 if (len >=numInt) {dir = "increase";}
 len = abs(len-numInt)
 if (loops>500){ len=0;}
 if ((loops%%100)==0){ split = split/2; }
}
write(c(loops,cutoff,len),"")

## Output using cutoffs:
indices = (simTOMcutoff[1,]>cutoff)
datout=cbind(
       rep(GeneSet[1], length(ProbeSet[indices])),
       GeneSet[indices],rep(0, length(ProbeSet[indices])),
       rep("M1002", length(ProbeSet[indices])),
       simTOM[1,][indices], Pearson[1,][indices])
colnames(datout) = c("gene1","gene2","zero","color","TO","Correlation")
for(i in seq(2,length(ProbeSet),by=1)){
 indices = (simTOMcutoff[i,]>cutoff)
 datout=rbind(datout,cbind(
        rep(GeneSet[i], length(ProbeSet[indices])),
        GeneSet[indices],rep(0, length(ProbeSet[indices])),
        rep("M1002", length(ProbeSet[indices])),
        simTOM[i,][indices], Pearson[i,][indices]))
}
datout = datout[order(datout[,5],decreasing=TRUE),]
fn = paste(moduleColor,"_visantOverall.csv",sep="")
write.table(datout,file=fn,sep=",",row.names=F, col.names=c("gene1","gene2",
            "zero","color","TO","Correlation"))
write(paste(fn, "written."),"")
}

isFirst = function(x){
 out = rep(TRUE,length(x));
 for (i in 2:length(x)) if(is.element(x[i],x[1:(i-1)])) out[i]=FALSE
 return(out)
}

####################################################################################################
####################################################################################################
####################################################################################################
#### The remainder of the functions are used for generating the heat map gene expression plots. ####
####################################################################################################
####################################################################################################
####################################################################################################


meanNA   = function(x) return(mean(x,na.rm=TRUE))     # FUNCTION FOR CALCULATING MEAN WITH NAs
sdNA     = function(x) return(sd(x,na.rm=TRUE))       # FUNCTION FOR CALCULATING SD WITH NAs
medianNA = function(x) return(median(x,na.rm=TRUE))   # FUNCTION FOR CALCULATING MEDIAN WITH NAs

plotMacaqueCortex2 <- function(inputPPA,layer,region,age,layerPositions,regionPositions,ageOffsets,
  plotTitle="CortexPlot",scaleA=FALSE,isLog2=TRUE,combineFn="medianNA",quantileScale = c(0,1),
  linearOrLog="linear",suppressAge=FALSE, bgPar="white", naBoxFile=NA, naBoxCol="lightyellow"){
  # This is a wrapper for the main function that allows all data to be included in a single vector
  # and has slightly different defaults (it assumes the Merck data is already scaled).
  inputPP = inputPPA[age!="Adult"]
  inputA  = inputPPA[age=="Adult"]
  plotMacaqueCortex(inputPP,inputA,layer,region,age,layerPositions,regionPositions,ageOffsets,
  plotTitle,scaleA,isLog2,combineFn,quantileScale,linearOrLog,suppressAge,bgPar,naBoxFile)
  
}
  
plotMacaqueCortex <- function(inputPP,inputA,layer,region,age,layerPositions,regionPositions,ageOffsets,
  plotTitle="CortexPlot",scaleA=FALSE,isLog2=TRUE,combineFn="medianNA",quantileScale = c(0,1),
  linearOrLog="linear",suppressAge=FALSE, bgPar="white",naBoxFile=NA, naBoxCol="lightgrey",
  showAdultDLG = FALSE, medianVals=NA){
  
  ## Format and subset the data
  
  inputPP  = as.numeric(inputPP);      inputA = as.numeric(inputA)
  layer    = as.character(layer);      layer[is.na(layer)]   = "none"
  age      = as.character(age);        age[is.na(age)]       = "none"
  region   = as.character(region);     region[is.na(region)] = "none"
  regInitial = region
  region[substr(age,1,1)=="E"] = paste(region[substr(age,1,1)=="E"],"_prenatal",sep="")
  region[substr(age,nchar(age),nchar(age))=="M"] = paste(region[substr(age,nchar(age),nchar(age))=="M"],"_postnatal",sep="")
  region[age=="Adult"] = paste(region[age=="Adult"],"_adult",sep="")  
  
  if(!showAdultDLG) {
   regionPositions = regionPositions[rownames(regionPositions)!="DLG_adult",] }
  
  if((!isLog2)&(linearOrLog!="linear")){ 
   inputPP = log2(inputPP);
   inputA  = log2(inputA);
  }
  if((isLog2)&(linearOrLog=="linear")){ 
   inputPP = 2^inputPP;
   inputA  = 2^inputA;
  }
  if(scaleA[1]!=FALSE) {
   regLay  = paste(regInitial,layer)
   kpA     = is.element(regLay[age=="Adult"],regLay[age!="Adult"][is.element(age[age!="Adult"],scaleA)])
   kpPP    = is.element(regLay[age!="Adult"],regLay[age=="Adult"])&is.element(age[age!="Adult"],scaleA)
   inputA  = inputA - mean(inputA[kpA],na.rm=TRUE) + mean(inputPP[kpPP],na.rm=TRUE)
  }
  inputPPA = c(inputPP,inputA)
  
  kpLayer  = is.element(layer,rownames(layerPositions))
  kpRegion = is.element(region,rownames(regionPositions))
  kpAge    = is.element(age,rownames(ageOffsets))
  kp       = kpLayer&kpRegion&kpAge
  inputPPA = inputPPA[kp]
  layer    = layer[kp]
  region   = region[kp]
  age      = age[kp]
  
  
  ## Compine all replicate samples (within each age/layer/region) using the input function
  
  ageLayReg  = paste(age,layer,region,sep="%")
  inputPPA   = cbind(inputPPA,inputPPA)
  inputPPA2  = findFromGroups(inputPPA,ageLayReg, match.fun(combineFn))
  ageLayReg  = colnames(inputPPA2)
  ageLayReg2 = strsplit(ageLayReg,"%")
  inputPPA2  = as.numeric(inputPPA2[1,])
  layer2 <- age2 <- region2 <- rep("A",length(inputPPA2))
  for (i in 1:length(age2)){
   age2[i]    = ageLayReg2[[i]][1]
   layer2[i]  = ageLayReg2[[i]][2]  
   region2[i] = ageLayReg2[[i]][3]
  }
  
  
  ## Quantile scale the data as requested by user input

  inputPPA4 <- inputPPA3 <- inputPPA2
  qS = as.numeric(quantile(inputPPA2,quantileScale,na.rm=TRUE))
  inputPPA4[inputPPA4<qS[1]] = qS[1]
  inputPPA4[inputPPA4>qS[2]] = qS[2]
  
  
  ## Prepare data to plot the expression levels
  
  rectBT = layerPositions[layer2,]
  rectLR = regionPositions[region2,3:4] + ageOffsets[age2,]
  rectB  = as.numeric(rectBT[,1])
  rectT  = as.numeric(rectBT[,2])
  rectL  = as.numeric(rectLR[,1])
  rectR  = as.numeric(rectLR[,2])
  textX  = (rectL+rectR)/2
  textY  = (rectB+rectT)/2
  if(suppressAge) for (a in unique(age2))
    inputPPA4[age2==a] = inputPPA4[age2==a]/max(inputPPA4[age2==a])
  
  rectCol   = numbers2colors(inputPPA4,signed=FALSE,centered=TRUE)
  
  
  ## The main plot
  
  par(bg=bgPar);
  plot(0,0,col="white",xlim=c(ageOffsets["E40",1]-2.5,max(rectLR)+0.1),ylim=c(min(rectBT)-2,max(rectBT)+1.5),
      axes=FALSE,xlab="",ylab="",main=plotTitle)
  rect(rectL,rectB,rectR,rectT,col=rectCol,border="black")
  text(textX,textY,layer2,cex=0.7)
  
  # Plot the NA data, if provided
  if (!is.na(naBoxFile)){
    ageLayRegNA = read.csv(naBoxFile) # Four columns: Age, Layer, Region, Text of missing data
	ageLayRegNA = as.matrix(ageLayRegNA)
	regTmp = ageLayRegNA[,3]
	ageTmp = ageLayRegNA[,1]
    regTmp[substr(ageTmp,1,1)=="E"] = paste(regTmp[substr(ageTmp,1,1)=="E"],"_prenatal",sep="")
    regTmp[substr(ageTmp,nchar(ageTmp),nchar(ageTmp))=="M"] = paste(regTmp[substr(ageTmp,nchar(ageTmp),nchar(ageTmp))=="M"],"_postnatal",sep="")
    regTmp[ageTmp=="Adult"] = paste(regTmp[ageTmp=="Adult"],"_adult",sep="") 
  	rectBTna = layerPositions[ageLayRegNA[,2],]
    rectLRna = regionPositions[regTmp,3:4] + ageOffsets[ageTmp,]
    rectBna  = as.numeric(rectBTna[,1])
    rectTna  = as.numeric(rectBTna[,2])
    rectLna  = as.numeric(rectLRna[,1])
    rectRna  = as.numeric(rectLRna[,2])
    textXna  = (rectLna+rectRna)/2
    textYna  = (rectBna+rectTna)/2
	rect(rectLna,rectBna,rectRna,rectTna,col=rep(naBoxCol,length(rectTna)),border="black")
    text(textXna,textYna,ageLayRegNA[,4],cex=0.7)
  }   # End plot the NA data
  
  xOff = as.numeric(ageOffsets[,1]);  xOff = c(xOff,max(rectR));  
  xOff = (xOff[1:(length(xOff)-1)] + xOff[2:length(xOff)])/2
  text(xOff,max(rectBT)+1,rownames(ageOffsets))
  abline(h=max(rectBT)+0.5)
  abline(v=ageOffsets[,1]-0.05); 
  abline(v=max(rectR)+0.05)
  

  for (a in unique(age2)){
    reg = unique(region2[age2==a])
    text(regionPositions[reg,5]+ageOffsets[a,],regionPositions[reg,6],regionPositions[reg,1],srt=90)
  }

  
  ## Clean up the plot
  
  rect(-1000,-1000,1000,layerPositions["Bottom",1],col="white",border="white")
  rect(ageOffsets["0M",1],-1000,ageOffsets["Adult",1]-0.1,regionPositions["V1_postnatal",6]-1,
    col="white",border="white")
  segments(ageOffsets["0M",1]-0.05,regionPositions["V1_postnatal",6]-1,ageOffsets["Adult",1]-0.05,
    regionPositions["V1_postnatal",6]-1,col="black")
  segments(ageOffsets["Adult",1]-0.05,regionPositions["V1_adult",6]-1,1000,
   regionPositions["V1_adult",6]-1,col="black")
  segments(-1000,layerPositions["Bottom",1],1000,layerPositions["Bottom",1],col="black")
  rect(max(rectR)+0.06,-1000,1000,1000,col="white",border="white")
  segments(-1000,layerPositions["Bottom",1],ageOffsets["0M",1]-0.05,layerPositions["Bottom",1])
  segments(-1000,layerPositions["Bottom",2],ageOffsets["0M",1]-0.05,layerPositions["Bottom",2])
	  
  xlab = paste(c("Min","25%","Median","75%","Max"),"=",signif(quantile(inputPPA4,na.rm=TRUE),2))
  xlab = paste(xlab,collapse=";  ")
  if (!is.na(medianVals[1])){
    perQ = round(100*mean(median(inputPPA2,na.rm=TRUE)>medianVals))
	totQ = paste(c("25%","50%","75%","Max"),"=",signif(quantile(medianVals,na.rm=TRUE),2)[2:5],sep="")
	totQ = paste(totQ,collapse=", ")
	xlab = paste(xlab,"        ===>        Expressed higher than ",perQ,"% of genes.  Distribution: ",totQ,sep="")
  }
  text((max(rectLR)+ageOffsets["E40",1]-2.5)/2,layerPositions["Bottom",1]-0.7,xlab,cex=1.1)
  
 
  ## Add a side plot summarizing expression by LAYER
  
  # ------- Update this plot to include postmitotic layers and germinal layers separately
  # ------- Appropriately align the lines to the middle of the box
  
  text(ageOffsets["E40",1]-1.25,regionPositions["V1_prenatal",6],"Layer Summary",cex=1.3)
  yBin  = 1/50
  yMean = max(rectBT)-sort((1:((max(rectBT)-min(rectBT))/yBin)*yBin))
  xMean = yMean*NA
  for (i in (2:length(yMean)))
    xMean[i] = meanNA(inputPPA2[(rectB<yMean[i])&(rectT>=yMean[i])])
  xMean = xMean-min(0,min(xMean,na.rm=TRUE))
  xMean[is.na(xMean)] = 0
  xMean = ageOffsets["E40",1]-0.2-2*(xMean)/max(xMean,na.rm=TRUE)
  segments(xMean,yMean,ageOffsets["E40",1]-0.2,yMean,lwd=2,col="blue")
  kp = (age2=="E90")&(region2=="V1_prenatal")
  text(ageOffsets["E40",1]-2.5,textY[kp],layer2[kp])
  text(ageOffsets["E40",1]-2.5,textY[(age2=="E70")&(region2=="CGE_prenatal")],"GE")

  
  ## Add a side plot summarizing expression by REGION (only show ACG and V1)
  
  ## Boundaries of boxes (for the remaining plots)
  minPlots    = layerPositions["Bottom",1]
  maxAgePlots = regionPositions["V1_postnatal",6]-1
  maxRCPlots  = regionPositions["V1_adult",6]-1
  ltAgePlots  = ageOffsets["0M",1]
  midPlots    = ageOffsets["Adult",1]
  rtRCPlots   = max(rectR)

  ## Function for plotting line plot in box
  plotInBox = function(means, sems, t, b, l, r, lab, main="",lrBord = 0.3, tbBord = 0.4, textOffset=0.6){
     l = l + lrBord;  r = r-lrBord;  t = t-tbBord;   b = b+tbBord;   ln = length(means)
	 xPos = (0:(ln-1))*(r-l)/(ln-1)+l;             val = round(max(means+sems))
	 text(xPos,b,lab,cex=0.6)
	 text((l+r)/2,t,main,cex=1.2)
	 b = b+textOffset;  t = t-textOffset
	 sc = (t-b)/max(means+sems)
	 means = means*sc+b;  sems = sems*sc
	 lines(xPos,means+sems,col="grey")
	 lines(xPos,pmax(means-sems,rep(b,length(means))),col="grey")
	 lines(xPos,means,lwd=3)
	 segments(l,b,r,b); segments(l,t,r,t); segments(l,b,l,t); segments(r,b,r,t)
	 l = l-lrBord/2
	 text(l,b,0,srt=90);  text(l,t-(tbBord/2)*nchar(as.character(val)),val,srt=90);
  } 
  
  regEg  = c("ACG","S1","V1")
  regEp  = c("ACG","V1")
  regA   = c("ACG","OG","dlPF","RG","V2","V1")
  
  omitLay= c(unique(layer2)[grep("4",unique(layer2))],"WM","L1","OFZ","IFZ","TMZ","ICD","L","C","Lcx","M")
  isGerm = rectT< -9.1
  kpEp   = (substr(age2,1,1)=="E")&(!is.element(layer2,omitLay))&(!isGerm)
  kpEg   = (substr(age2,1,1)=="E")&(!is.element(layer2,omitLay))&(isGerm)
  kpA    = (substr(age2,nchar(age2),nchar(age2))=="M")&(!is.element(layer2,omitLay))
  yMeanEp= findFromGroups(cbind(inputPPA2[kpEp],inputPPA2[kpEp]),regionPositions[region2[kpEp],1],meanNA)[1,regEp]
  ySDEp  = findFromGroups(cbind(inputPPA2[kpEp],inputPPA2[kpEp]),regionPositions[region2[kpEp],1],sdNA)[1,regEp]
  yMeanEg= findFromGroups(cbind(inputPPA2[kpEg],inputPPA2[kpEg]),regionPositions[region2[kpEg],1],meanNA)[1,regEg]
  ySDEg  = findFromGroups(cbind(inputPPA2[kpEg],inputPPA2[kpEg]),regionPositions[region2[kpEg],1],sdNA)[1,regEg]
  yMeanA = findFromGroups(cbind(inputPPA2[kpA],inputPPA2[kpA]),regionPositions[region2[kpA],1],meanNA)[1,regA]
  ySDA   = findFromGroups(cbind(inputPPA2[kpA],inputPPA2[kpA]),regionPositions[region2[kpA],1],sdNA)[1,regA]
  #ySDEp  = ySDEp/sqrt(sum(kpEp));    ySDEg  = ySDEg/sqrt(sum(kpEg));    ySDA  = ySDA/sqrt(sum(kpA)); # To use SEM, uncomment line
  
  # halfRCtb = (minPlots+maxRCPlots)/2;    halfRClr = (midPlots+rtRCPlots)/2
  # plotInBox(yMeanA,ySDA,halfRCtb,minPlots,midPlots,rtRCPlots,regA,"Postnatal")
  # plotInBox(yMeanEg,ySDEg,maxRCPlots,halfRCtb,midPlots,halfRClr,regEg,"Germinal")
  # plotInBox(yMeanEp,ySDEp,maxRCPlots,halfRCtb,halfRClr,rtRCPlots,regEp,"Postmitotic")
  
  
  ## Add a side plot summarizing expression by AGE
  
  # agep   = rownames(ageOffsets);   agee = agep[1:6]
  # yMeanEp= findFromGroups(cbind(inputPPA2[!isGerm],inputPPA2[!isGerm]),age2[!isGerm],meanNA)[1,agep]
  # ySDEp  = findFromGroups(cbind(inputPPA2[!isGerm],inputPPA2[!isGerm]),age2[!isGerm],sdNA)[1,agep]
  # yMeanEg= as.data.frame(findFromGroups(cbind(inputPPA2[isGerm],inputPPA2[isGerm]),age2[isGerm],meanNA))[1,agee]
  # ySDEg  = as.data.frame(findFromGroups(cbind(inputPPA2[isGerm],inputPPA2[isGerm]),age2[isGerm],sdNA))[1,agee]
  # yMeanEg= c(as.numeric(yMeanEg),rep(0,length(agep)-length(agee)))
  # ySDEg  = c(as.numeric(ySDEg),rep(0,length(agep)-length(agee)))
  
  # halfAgelr = (midPlots+ltAgePlots)/2
  # plotInBox(yMeanEg,ySDEg,maxAgePlots,minPlots,ltAgePlots,halfAgelr,agep,"Germinal")
  # plotInBox(yMeanEp,ySDEp,maxAgePlots,minPlots,halfAgelr,midPlots,agep,"Postmitotic")
  # text(ltAgePlots+(halfAgelr-ltAgePlots)*0.75,(minPlots+maxAgePlots)/2,"(no data)")
  
  ## Create a legend for plotting in the upper left corner
  
  legendVal = min(inputPPA4,na.rm=TRUE)+(0:5)*(max(inputPPA4,na.rm=TRUE)-min(inputPPA4,na.rm=TRUE))/5
  legendCol = numbers2colors(legendVal,signed=FALSE,centered=TRUE)
  legendVal = round(legendVal)
  xOff = ageOffsets["E40",1]-2.5
  xOff = xOff+(0:(length(legendVal)-1)*(2.5/length(legendVal)))
  points(xOff,rep(max(rectBT)+1.5,length(xOff)),pch=15,cex=4,col=legendCol)
  text(xOff,rep(max(rectBT)+1.5,length(xOff)),legendVal,srt=90)
  
}


plotMacaqueCortexSmall <- function(inputPP,layer,age,layerPositionsS,agePositionsS,
  plotTitle="CortexPlot",isLog2=TRUE,combineFn="medianNA",quantileScale = c(0,1),
  linearOrLog="linear",bgPar="white",displayLayers=FALSE,legendPos=NULL,outputDataOnly=FALSE){
  
  ## Format and subset the data
  inputPP  = as.numeric(inputPP);      
  layer    = as.character(layer);      layer[is.na(layer)]   = "none"
  age      = as.character(age);        age[is.na(age)]       = "none"
  if((!isLog2)&(linearOrLog!="linear"))   inputPP = log2(inputPP);
  if((isLog2)&(linearOrLog=="linear"))    inputPP = 2^inputPP;
  kpLayer  = is.element(layer,rownames(layerPositionsS))
  kpAge    = is.element(age,rownames(agePositionsS))
  kp       = kpLayer&kpAge
  inputPP  = inputPP[kp]
  layer    = layer[kp]
  age      = age[kp]
  
  ## Combine all replicate samples (within each age/layer/region) using the input function
  ageLayReg  = paste(age,layer,sep="%")
  inputPP    = cbind(inputPP,inputPP)
  inputPP2   = findFromGroups(inputPP,ageLayReg, match.fun(combineFn))
  ageLayReg  = colnames(inputPP2)
  ageLayReg2 = strsplit(ageLayReg,"%")
  inputPP2   = as.numeric(inputPP2[1,])
  if (outputDataOnly) { names(inputPP2) = ageLayReg;  return(inputPP2); }
  layer2 <- age2 <- rep("A",length(inputPP2))
  for (i in 1:length(age2)){
   age2[i]    = ageLayReg2[[i]][1]
   layer2[i]  = ageLayReg2[[i]][2]  
  }
  agePositionsS = as.matrix(agePositionsS)
  layerPositionsS = as.matrix(layerPositionsS)
  positions  = cbind(agePositionsS[age2,1:2],layerPositionsS[layer2,1:2])
  labelX     = c(agePositionsS[,3],layerPositionsS[,3])
  labelY     = c(agePositionsS[,4],layerPositionsS[,4])
  labelText  = c(rownames(agePositionsS),rownames(layerPositionsS))
  textData   = layer2
  if(!displayLayers) textData = rep("",length(textData))
  
  ## Plot the expression levels in the appropriate positions on the plot
  par(bg=bgPar);
  qS      = as.numeric(quantile(inputPP2,quantileScale,na.rm=TRUE))
  plotRectangle(inputPP2, textData, positions, quantileScale=quantileScale, main=plotTitle, 
  colorRange = qS, legendPos=legendPos, labelX=labelX, labelY=labelY, labelText = labelText,
   colVector = blueWhiteRed(100)[51:100],signed=FALSE, numDecimals=0)
}



## Create a less complicated summary plot with rectangles
plotRectangle <- function(colorData, textData, positions, quantileScale=c(0,1), main="plot", 
  colorRange = range(colorData), legendPos = c(8.5,-13.5,-10), labelX=NULL, labelY=NULL,
  labelText = NULL, colVector = blueWhiteRed(100)[51:100],signed=FALSE,numDecimals=1,rectBorder="black") {

  qS      = as.numeric(quantile(colorData,quantileScale,na.rm=TRUE))
  colorData[colorData<qS[1]] = qS[1];    colorData[colorData>qS[2]] = qS[2]
  rectB   = as.numeric(positions[,3]);   rectT  = as.numeric(positions[,4])
  rectL   = as.numeric(positions[,1]);   rectR  = as.numeric(positions[,2])
  textX   = (rectL+rectR)/2;             textY  = (rectB+rectT)/2
  colTmp  = c(colorData,colorRange) 
  rectTmp = numbers2colors(colTmp,signed=signed,centered=TRUE,colors=colVector)
  rectCol = rectTmp[1:(length(rectTmp)-2)]
  
  ## The main plot
  plot(0,0,col="white",xlim=range(rectL,rectR,labelX),ylim=range(rectT,rectB,labelY),
    axes=FALSE,xlab="",ylab="",main=main)
  rect(rectL,rectB,rectR,rectT,col=rectCol,border=rectBorder)
  text(textX,textY,textData);
  
  ## Plot the legend?
  if(!is.null(legendPos[1])){
   legendVal = min(colTmp,na.rm=TRUE)+(0:5)*(max(colTmp,na.rm=TRUE)-min(colTmp,na.rm=TRUE))/5
   legendCol = numbers2colors(legendVal,signed=signed,centered=TRUE,colors=colVector)
   legX      = rep(legendPos[1],6);
   legY      = quantile(legendPos[2:3],probs=(0:5)/5)
   points(legX,legY,pch=15,cex=4,col=legendCol)
   text(legX,legY,round((10^numDecimals)*legendVal)/(10^numDecimals))
  }
  
  ## Plot the labels
  if (!is.null(labelX[1])) text(labelX,labelY,labelText)
  
}



plotSDvsLength <- function(varDat,enrichGene,allGenes,geneLengths,bins=20,col="green",
 xlab="Samples sorted by gene length",ylab="Number of genes",...){

 if(!is.null(dim(varDat)[1]))  varDat  = apply(varDat,1,sd)
 isEnr   = is.element(allGenes,enrichGene)
 isLong  = is.element(allGenes,names(sort(-geneLengths))[1:sum(isEnr)])
 names(isEnr) <- names(isLong) <- allGenes
 perBin  = round(length(allGenes)/bins)
 groups  = rep(bins,length(allGenes))
 for (i in 1:(bins-1)) groups[(perBin*(i-1)+1):(perBin*i)] = i 

 enrDat  = isEnr[names(sort(varDat))]
 lenOrd  = isLong[names(sort(varDat))]
 y       = as.numeric(table(factor(groups)[enrDat]))
 r       = as.numeric(table(factor(groups)[lenOrd]))
 plot(1:bins,y,pch=19,col=col,ylim=c(0,max(c(y,r))),xlab=xlab,ylab=ylab,...)
 lines(lowess(1:bins,y),col=col,lwd=2)
 points(1:bins,r,pch=19,col="grey")
 lines(lowess(1:bins,r),col="grey")
}




plotMacaqueCortexStats <- function(inputPP,layer,region,age,layerPositions,regionPositions,ageOffsets,
  plotTitle="CortexPlot",isLog2=FALSE,combineFn="medianNA",quantileScale = c(0,1),
  linearOrLog="linear",suppressAge=FALSE, bgPar="white",showAdultDLG = FALSE, medianVals=NA, boxText=layer){
  
  ## Format and subset the data
  
  inputPP  = as.numeric(inputPP);      
  layer    = as.character(layer);      layer[is.na(layer)]   = "none"
  age      = as.character(age);        age[is.na(age)]       = "none"
  region   = as.character(region);     region[is.na(region)] = "none"
  regInitial = region
  region[substr(age,1,1)=="E"] = paste(region[substr(age,1,1)=="E"],"_prenatal",sep="")
  region[substr(age,nchar(age),nchar(age))=="M"] = paste(region[substr(age,nchar(age),nchar(age))=="M"],"_postnatal",sep="")
  region[age=="Adult"] = paste(region[age=="Adult"],"_adult",sep="")  
  
  if(!showAdultDLG) {
   regionPositions = regionPositions[rownames(regionPositions)!="DLG_adult",] }
  
  if((!isLog2)&(linearOrLog!="linear")){ 
   inputPP = log2(inputPP);  }
  if((isLog2)&(linearOrLog=="linear")){ 
   inputPP = 2^inputPP;  }
   
  kpLayer  = is.element(layer,rownames(layerPositions))
  kpRegion = is.element(region,rownames(regionPositions))
  kpAge    = is.element(age,rownames(ageOffsets))
  kp       = kpLayer&kpRegion&kpAge
  inputPP  = inputPP[kp]
  layer2   = layer[kp]
  region2  = region[kp]
  age2     = age[kp]
  boxText  = boxText[kp]
  
  ## Quantile scale the data as requested by user input

  inputPPA4 <- inputPPA2 <- inputPP
  qS = as.numeric(quantile(inputPP,quantileScale,na.rm=TRUE))
  inputPPA4[inputPPA4<qS[1]] = qS[1]
  inputPPA4[inputPPA4>qS[2]] = qS[2]
  
  
  ## Prepare data to plot the expression levels
  
  rectBT = layerPositions[layer2,]
  rectLR = regionPositions[region2,3:4] + ageOffsets[age2,]
  rectB  = as.numeric(rectBT[,1])
  rectT  = as.numeric(rectBT[,2])
  rectL  = as.numeric(rectLR[,1])
  rectR  = as.numeric(rectLR[,2])
  textX  = (rectL+rectR)/2
  textY  = (rectB+rectT)/2
  if(suppressAge) for (a in unique(age2))
    inputPPA4[age2==a] = inputPPA4[age2==a]/max(inputPPA4[age2==a])
  
  rectCol   = numbers2colors(inputPPA4,signed=FALSE,centered=TRUE)
  
  
  ## The main plot
  
  par(bg=bgPar);
  plot(0,0,col="white",xlim=c(ageOffsets["E40",1]-2.5,max(rectLR)+0.1),ylim=c(min(rectBT)-2,max(rectBT)+1.5),
      axes=FALSE,xlab="",ylab="",main=plotTitle)
  rect(rectL,rectB,rectR,rectT,col=rectCol,border="black")
  text(textX,textY,boxText,cex=0.7)
  
  xOff = as.numeric(ageOffsets[,1]);  xOff = c(xOff,max(rectR));  
  xOff = (xOff[1:(length(xOff)-1)] + xOff[2:length(xOff)])/2
  text(xOff,max(rectBT)+1,rownames(ageOffsets))
  abline(h=max(rectBT)+0.5)
  abline(v=ageOffsets[,1]-0.05); 
  abline(v=max(rectR)+0.05)
  

  for (a in unique(age2)){
    reg = unique(region2[age2==a])
    text(regionPositions[reg,5]+ageOffsets[a,],regionPositions[reg,6],regionPositions[reg,1],srt=90)
  }

  
  ## Clean up the plot
  
  rect(-1000,-1000,1000,layerPositions["Bottom",1],col="white",border="white")
  rect(ageOffsets["0M",1],-1000,ageOffsets["Adult",1]-0.1,regionPositions["V1_postnatal",6]-1,
    col="white",border="white")
  segments(ageOffsets["0M",1]-0.05,regionPositions["V1_postnatal",6]-1,ageOffsets["Adult",1]-0.05,
    regionPositions["V1_postnatal",6]-1,col="black")
  segments(ageOffsets["Adult",1]-0.05,regionPositions["V1_adult",6]-1,1000,
   regionPositions["V1_adult",6]-1,col="black")
  segments(-1000,layerPositions["Bottom",1],1000,layerPositions["Bottom",1],col="black")
  rect(max(rectR)+0.06,-1000,1000,1000,col="white",border="white")
  segments(-1000,layerPositions["Bottom",1],ageOffsets["0M",1]-0.05,layerPositions["Bottom",1])
  segments(-1000,layerPositions["Bottom",2],ageOffsets["0M",1]-0.05,layerPositions["Bottom",2])
	  
  xlab = paste(c("Min","25%","Median","75%","Max"),"=",signif(quantile(inputPPA4,na.rm=TRUE),2))
  xlab = paste(xlab,collapse=";  ")
  if (!is.na(medianVals[1])){
    perQ = round(100*mean(median(inputPPA2,na.rm=TRUE)>medianVals))
	totQ = paste(c("25%","50%","75%","Max"),"=",signif(quantile(medianVals,na.rm=TRUE),2)[2:5],sep="")
	totQ = paste(totQ,collapse=", ")
	xlab = paste(xlab,"        ===>        Expressed higher than ",perQ,"% of genes.  Distribution: ",totQ,sep="")
  }
  text((max(rectLR)+ageOffsets["E40",1]-2.5)/2,layerPositions["Bottom",1]-0.7,xlab,cex=1.1)
  
  ## Add a side plot summarizing expression by LAYER
  
  # ------- Update this plot to include postmitotic layers and germinal layers separately
  # ------- Appropriately align the lines to the middle of the box
  
  text(ageOffsets["E40",1]-1.25,regionPositions["V1_prenatal",6],"Layer Summary",cex=1.3)
  yBin  = 1/50
  yMean = max(rectBT)-sort((1:((max(rectBT)-min(rectBT))/yBin)*yBin))
  xMean = yMean*NA
  for (i in (2:length(yMean)))
    xMean[i] = meanNA(inputPPA2[(rectB<yMean[i])&(rectT>=yMean[i])])
  xMean = xMean-min(0,min(xMean,na.rm=TRUE))
  xMean[is.na(xMean)] = 0
  xMean = ageOffsets["E40",1]-0.2-2*(xMean)/max(xMean,na.rm=TRUE)
  segments(xMean,yMean,ageOffsets["E40",1]-0.2,yMean,lwd=2,col="blue")
  kp = (age2=="E90")&(region2=="V1_prenatal")
  text(ageOffsets["E40",1]-2.5,textY[kp],layer2[kp])
  text(ageOffsets["E40",1]-2.5,textY[(age2=="E70")&(region2=="CGE_prenatal")],"GE")

  
  ## Create a legend for plotting in the upper left corner
  
  legendVal = min(inputPPA4,na.rm=TRUE)+(0:5)*(max(inputPPA4,na.rm=TRUE)-min(inputPPA4,na.rm=TRUE))/5
  legendCol = numbers2colors(legendVal,signed=FALSE,centered=TRUE)
  legendVal = round(legendVal)
  xOff = ageOffsets["E40",1]-2.5
  xOff = xOff+(0:(length(legendVal)-1)*(2.5/length(legendVal)))
  points(xOff,rep(max(rectBT)+1.5,length(xOff)),pch=15,cex=4,col=legendCol)
  text(xOff,rep(max(rectBT)+1.5,length(xOff)),legendVal,srt=90)
  
}



plotTwoGroupsSmoothed <- function (expr, age, groupVector, groups=unique(groupVector)){




}

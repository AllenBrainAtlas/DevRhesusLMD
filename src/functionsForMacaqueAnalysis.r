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
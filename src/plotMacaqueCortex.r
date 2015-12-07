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

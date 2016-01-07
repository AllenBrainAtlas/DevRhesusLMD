# Script to perform Multivariate Distance Matrix Regression (MDMR) in R

# Based in part on methods presented in:
# McArdle, B.H. and Anderson, M.J. (2001).
# "Fitting multivariate models to community data:
# A comment on distance-based redundancy analysis". Ecology, 82, 290-297.
#
# Salem RM, O'Connor DT, Schork NJ.
# "Curve-Based Multivariate Distance Matrix Regression (CMDMR) Analysis: 
# Application to Genetic Association Analyses Involving Repeated Measures". Physiol Genomics 2010; 42: 236-247 


# Written by Rany M. Salem

# Version 1.0.1 (BETA)
# Date 07.01.10 (Lastest Update)
# Code should run on all version of R, written and tested on version 2.6.0 & 2.7.0
# Program requires the following libraries: (must be downloaded and installed)
#    No Dependencies

# Key Features & Notes
#    1. Input files consists of square distance matrix (NxN) and predictor variables file (N x Vars)
#    2. Program checks if matrix is square and that number of subjects in dmatrix and predictor file are the same,
#	   Warning printed and program STOPS if mismatches found
#    3. Variable header row in predictor matrix is required and used in results output
#    4. Accomodates missing data in predictor file (excluding subjects with missing values), though not in distance matrix
#    5. Dummy variable coding of predictor variables can be specified (into categorical/factor variables), 
#	  with "cat" in 2nd line of predictor input file for variable into a custom set

#NOTE this is wrong/updated version of code does not use cat/con optiona
#    6. All variables can be analyzed as categorical ("cat") or continuous ("con") using DUMMY.VAR option, 2nd line 
#	  of predictor file must then contain predictor data

#    7. Performs Univariate (Marginal) test for each predictor
#    8. Performs Forward Regression (Conditional), forward selection adds variables to model until fit does not improved
#    9. Performs Multiple Regression, analysis of each variables adjusted for all other variables in the model
#    10. Allows for user specified custom model, selecting set of predictors, covariates and interaction terms, specified 
#	   variable names are checked against names in predictor file, error message printed if match not found

# List of Parameters: file name, file location, legend title, and range of variables for correlation matrix
#	DMATRIX			- Distance Matrix (NxN)
#	DATA			- DATA Matrix of Predictor Variables (N by # Wariables): 1st line variable header,
#	PERMS			- Number of permutations (for P-Value calculations), Default = 100

#	FORMULA			- Formula describing variables (similar format as other R functions)

#	UNIVARIATE		- Requests univariate tests (TRUE/FALSE), Default = FALSE
#	FOWARD.REG		- Requests forward regression (conditional) tests (TRUE/FALSE), Default = FALSE
#	MULTIVARIABLE		- Requests multivariable tests (TRUE/FALSE), each predictor is  adjusted for all other variables in model, Default = FALSE
#	USER.CUSTOM		- Requests user defined custom model (TRUE/FALSE), utilizing USER.PREDICT, USER.COVARS & USER.INTERACTION options. 
#				  Note: Program STOPS and Warning printed if variables specified are not in PREDICT file, Default = FALSE
#	USER.PREDICT  		- Specifies a set of variables for analysis, each variable is analyzed 1 at a time, with covariates 
#				  or interaction terms in model, Default = NULL, ignored if USER.CUSTOM not specified
#	USER.COVARS 	  	- Specifies set of variables to use as covariates/adjustment for each variable specified in USER.PREDCT. 
#			    	  Default = NULL, option ignored if USER.CUSTOM not specified
#  	USER.INTERACTION 	- Requests and calculates interaction term for Single CONTINUOUS variable with each variable listed in USER.PREDICT. 
#				  Internally generates model with Main Effect (predictor variable), interaction variable (USER.INTERACTION) and interaction 
#				  term, results only reported for Interaction term. If USER.COVARS specified, interaction model will be adjuated for covariates, 
#				  Interaction variable should not be in list of covariates, Default = NULL. Ignored if USER.CUSTOM not specified.
#	RESULTS.SORT		- Requests sorting of results of UNIVARIATE, MULTIVARIABLE and USER.CUSTOM models (TRUE/FALSE), Default = FALSE


#########################################################################################################
#Output Options
#########################################################################################################

options(scipen = 4) # bias against scientific notation
options(digits = 4) # show fewer decimal places

########################## Data Prep & Check (QC) ##########################


# Read in file distance matrix file (n x n)
#    dm_file<-read.table("C:/Documents and Settings/Rany/My Documents/Documents/Dissertation Project/D-work/P1 HV//Dmatrix/Gross-Curves/dist.hv.adaptneyman.txt",header=T)
#    dm_file<-read.table("D:/My Documents/Documents/Dissertation Project/D-work/P1 HV//Dmatrix/Gross-Curves/dist.hv.adaptneyman.txt",header=T)

    #convert to matrix
#   dmatrix<-as.matrix(dm_file)


# Read in file predictor file (n x vars)
# no missing values
#    pred_file<-read.table("C:/Documents and Settings/Rany/My Documents/Documents/Dissertation Project/D-work/P1 HV/Temp/xpred.txt",header=T, na.strings = ".")
#    pred_file<-read.table("D:/My Documents/Documents/Dissertation Project/D-work/P1 HV/Temp/xpred.txt",header=T, na.strings = ".")

# with missing values
#    pred_file<-read.table("C:/Documents and Settings/Rany/My Documents/Documents/Dissertation Project/D-work/P1 HV/Temp/wxpredict.txt",header=T, na.strings = ".")
#    pred_file<-read.table("D:/My Documents/Documents/Dissertation Project/D-work/P1 HV/Temp/wxpredict.txt",header=T, na.strings = ".")


# All HV predictors
#    pred_file<-read.table("C:/Documents and Settings/Rany/My Documents/Documents/Dissertation Project/D-work/P1 HV/Predictors/hv.xpredictors.txt",header=T, na.strings = ".")
#    pred_file<-read.table("D:/My Documents/Documents/Dissertation Project/D-work/P1 HV/Predictors/hv.predictors.txt",header=T, na.strings = ".")
#    pred_file<-read.table("C:/Documents and Settings/Rany/My Documents/Documents/Dissertation Project/D-work/P1 HV/Temp/zpred.txt",header=T, na.strings = ".")

    #convert to matrix
#    predict<-as.data.frame(pred_file)



#########################################################################################################
# TEST/DEVELOPMENT CODE (Beta)
#########################################################################################################

# Main Function: Performs MDMR Analysis, with all embedded internal fxns

MDMR<-function(FORMULA=NULL, DMATRIX=dmatrix, DATA=data, PERMS=100, 
		UNIVARIATE=FALSE, FORWARD.REG=FALSE, MULTIVARIABLE=FALSE, USER.CUSTOM=FALSE,
		USER.PREDICT=NULL, USER.COVARS=NULL, USER.INTERACTION=NULL, RESULT.SORT=FALSE)
	{

	# GOWER FXN: convert Distance Matrix into Gower's Centered Matrix, input must be a square distance matrix (Internal Sub Fxn)
	GOWER<-function(xG=DMATRIX)
		{
		dx  <- ncol(xG)                    			# Needed for calcs below
		xG2 <- (xG*xG)/(-2)                 		# Creating A Matrix to be used in constructing Gower's Centered Matrix
		aux <- (diag(dx) - matrix((1/dx), dx, dx))    	# Creating Auxiliary centered Matrix
		gmatrix <- ((aux%*%xG2)%*%aux)        		# Creating Gower's Centered Matrix
		return(gmatrix)
		}


	# Dataframe Sort FXN: sorts data frame by variable 'xzx' to insure proper order of data for later analysis (Internal Sub Fxn)
    	DF.SORT<-function (data, vars = "var_name")        
        	{
            	if (length(vars) == 0 || is.null(vars))    
            	return(data)
            	data[do.call("order", data[, vars, drop = FALSE]), , drop = FALSE]
            	}


	# Preparing variables for analysis: contrast coding of factors and set variables into list
	#	
	DUMMY.CODING<-function(xfrmla, data=DATA)
     		{
		if(class(xfrmla) =="formula")
        		{
 			Xformula<- xfrmla
			}
		else
			{
			# Generate formula notation from list of variable names
			Xformula<- as.formula(paste("y ~ ", paste(xfrmla, collapse= "+")))
			}

		#subsets starting dataset to variables used in analysis
		xz<-length(Xformula)-5
		sdata<-model.frame(Xformula[xz], data=data, na.action=na.pass, drop.unused.levels=TRUE)

		#generates contrasts for factors varaibles
  		options(contrasts = options()$contrasts)
		mdata<- model.matrix(Xformula[xz], sdata)
		mxdata<-as.data.frame(mdata[,-1])

		# gives structure of contrasts variables
		xlist<- attr(mdata, "assign")[-1]

		#retrieves names form constructed contrast variables	
		xname<-attr(attr(sdata, "terms"), "term.labels")


		DUMMY.FORMATTING<-function(nx49)		
			{
			temp.var<-as.data.frame(mxdata[,xlist %in% nx49])
      			nx<-(1:ncol(temp.var))         	# Determining n for dummy coding
                                    
			names(temp.var)<-c(as.vector(factor(nx[1:length(nx)], label=paste(xname[nx49],".", sep=""))))
			return(temp.var)
			}

		Cvars<-lapply(1:length(xname), DUMMY.FORMATTING)
		names(Cvars)<-c(xname)
		return(Cvars)
		}



	# "PVE_SORT" Function to recursively determine var w/ max PVE adjusting for all others (Internal Sub Fxn)
    	PVE.SORT<-function(dataX=pred.vars)
      	{
        	dxd<-dataX      	# Defining analysis dataset
        	npass <- c()    	# Defining Empty Results file, for sorted PVE values
        	j <- 0          	# Defining start value
        	jj<-length(dxd)   # "jj" used to define max # runs

		interZ<-interX    # Defining intercept vector w/ value = "1" (used below)

		while (j < jj)    # While loop
            	{
            	if(j<1) intercept<-cbind(interZ) else intercept<-cbind(interZ,as.data.frame(npass))
            	rrr<-MDMR1(m1DATA=dxd,INTRCPT=intercept,FXN="PVE")
            	maxP<-which.max(rrr)
            	npass<-c(npass,dxd[maxP])
            	dxd<-dxd[-maxP]
            	j <- length(npass)
            	}
        	return(npass)
		}



	# "mdmr1" Function to calculate FSTAT, PVE, and P-VALUE w/ Permutation Testing
	MDMR1<-function(m1DATA=ddd, INTRCPT=intercept, FXN="PVE", xPERMS=100)
		{
		MDMR2<-function(m2DATA)
			{

			########################## Predictor Variable File Preparation ##########################

        		varP<-as.data.frame(m2DATA)    
        		decomposition_matrix<-cbind(INTRCPT, varP)  	# Combining intercept + predictor variable into single file    
        		p3<-ncol(decomposition_matrix)            	# Needed for F-Stat & Proportion Variance Explained calcs below

			# Assessing predictor file for missing variables (determining missing subjects)
			mmm<-(1:nrow(decomposition_matrix))[grep("FALSE",(complete.cases(decomposition_matrix)))]

			# Removing subjects with missing values from both predictor and dmatrix file
        		if(length(mmm)==0) decomposition_matrix<-decomposition_matrix 	else decomposition_matrix<-as.matrix(decomposition_matrix)[-mmm,]
        		if(length(mmm)==0) gower_matrix_F<-gower_matrix 			else gower_matrix_F<-as.matrix(gower_matrix)[-mmm, -mmm]

        		p2<-nrow(gower_matrix_F)        # Needed for F-Stat & PVE Calcs below

	  		# Full Model
			########################## STEP2: Calculating Single Regression ##########################
        		########################## QR Decomposition to calculate statistics ######################
    
			qr_docompX<-qr(decomposition_matrix, tol=1e-7, LAPACK=FALSE)    	# Performing QR Decomposition
        		qhat_matrix<-qr.Q(qr_docompX)                          		# Extracting Q of QR decomposition for analysis
        		hat_matrix<-qhat_matrix%*%(t(qhat_matrix))                		# Calculate hat matrix using the decomposition
			id_minus_hat<-(diag(ncol(hat_matrix))-hat_matrix)			# Calculate (ID matrix - hat_matrix) to be used in calc of F stat

        		SS_F<-sum(diag(((hat_matrix%*%gower_matrix_F)%*%hat_matrix)))   # Calculating SS(Trace)
        
        		if(FXN=="PVE")
            		{
            		########################## Calculate Proportion of Variance Explained ##########################
            		PVE<-SS_F/(sum(diag(gower_matrix_F)))
            		return(PVE)
            		}
        		if(FXN=="MDMR.FXN")
            		{    
            		PVE<-SS_F/(sum(diag(gower_matrix_F)))

				######################### Calculate (psuedo) F-Statistic ##########################
            		fstat_a<-(SS_F/(p3-1))												
            		fstat_b<-(sum(diag((id_minus_hat%*%gower_matrix_F)%*%id_minus_hat)))/(p2-p3)
            		fstat<-fstat_a/fstat_b
            		Nobs<-nrow(gower_matrix_F)

            		######################## F-Statistic Permutation Testing ##########################
                		# Approximating p-value through permutation analysis
                		# Permutation analysis carried out by simultaneously permuting rows
                		# and columns of Gower's Center Matrix

            		# Fxn (internal) to Perform F-Stat permutation testing
            		FXN1<-function(x=gower_matrix_F, nx=xPERMS)
                			{
                			# FXN to Extract data from lower tri of distance matrix into vector and shuffle
                			FXN2<-function(x)
                    			{    
                    			mm<-x                
                    			if (ncol(mm) != nrow(mm))
                    			stop("Matrix not square.")
                    			xdata<- mm[col(mm) < row(mm)]
                    			diagx<-diag(x)        # Extract diagonal of matrix

			                	nn<-nrow(x)            
                    			full <- matrix(0, nn, nn)    # Shuffle and converting vector back into a symtric matrix (shuffled distance matrix)
			                	full[lower.tri(full)] <- sample(xdata, replace=FALSE)
	                    		full2 <- t(full)
                    			diag(full2) <- sample(diagx, replace=FALSE)    # Shuffling Diagonal of matrix    
                    			P_gower_matrix<-full + full2    # Generating Full Permuted gower_matrix

                    			P_fstat_a<-((sum(diag((hat_matrix%*%P_gower_matrix)%*%hat_matrix)))/(p3-1))						
			  			P_fstat_b<-((sum(diag((id_minus_hat%*%P_gower_matrix)%*%id_minus_hat)))/(p2-p3))				
			  			P_fstat<-P_fstat_a/P_fstat_b 											
                    			return(P_fstat)
                    			}
                			abc<-replicate(nx, FXN2(x))
                			return(length(which(abc>fstat))/length(abc))
                			}
 		       	pval<-FXN1(gower_matrix_F)
            		Nobs<-nrow(gower_matrix_F)
				results<-c(Nobs, xPERMS, SS_F, fstat, pval, PVE)												
            		return(results)
            		}
        		if(FXN=="MDMR.MULTI")
            		{    

		            # Reduced model
            		decomposition_matrix_R<-cbind(INTRCPT)      	# REDUCED MODEL only uses intercept (or Intercept + Covariates, w/ predictor of interest removed)
            		p3_R<-ncol(decomposition_matrix_R)           	# Needed for F-Stat & Proportion Variance Explained calcs below
            
            		# Removing subjects with missing values from both predictor and dmatrix file based on full model missing
            		if(length(mmm)==0) decomposition_matrix_R<-decomposition_matrix_R 	else decomposition_matrix_R<-as.matrix(decomposition_matrix_R)[-mmm,]
            		if(length(mmm)==0) gower_matrix_R<-gower_matrix 				else gower_matrix_R<-as.matrix(gower_matrix)[-mmm, -mmm]

            		p2_R<-nrow(gower_matrix_R)        # Needed for F-Stat & PVE Calcs below

            		# Reduced Model
		            ########################## STEP2: Calculating Single Regression ##########################
		            ########################## QR Decomposition to calculate statistics ######################
    
            		qr_docompX_R<-qr(decomposition_matrix_R, tol=1e-7, LAPACK=FALSE)    	# Performing QR Decomposition
		            qhat_matrix_R<-qr.Q(qr_docompX_R)                     			# Extracting Q of QR decomposition for analysis
		            hat_matrix_R<-qhat_matrix_R%*%(t(qhat_matrix_R))            		# Calculate hat matrix using the decomposition
		            id_minus_hat_R<-(diag(ncol(hat_matrix_R))-hat_matrix_R)         		# Calculate (ID matrix - hat_matrix) to be used in calc of F stat

		            SS_R<-sum(diag(((hat_matrix_R%*%gower_matrix_R)%*%hat_matrix_R)))		# Calculating SS(Trace)

		            ########################## Calculate Proportion of Variance Explained: Reduced Model, Cumulative, Added by New Var ##########################
		            PVE_R<-SS_R/sum(diag(gower_matrix_R))

		            CPVE<-SS_F/sum(diag(gower_matrix_F))  #Cumulative PVE (full model)

		            PVE<-CPVE-PVE_R              		# PVE of added variable (Full - Reduced Models)

		            ######################### Calculate (psuedo) F-Statistic ##########################
		            fstat_a<-SS_F-SS_R
		            fstat_b<-(sum(diag((id_minus_hat%*%gower_matrix_F)%*%id_minus_hat)))/(p2-p3)
		            fstat<-(fstat_a/(p3-p3_R))/fstat_b
		            Nobs<-nrow(gower_matrix_F)

		            ######################## F-Statistic Permutation Testing ##########################
	             	# Approximating p-value through permutation analysis
                		# Permutation analysis carried out by simultaneously permuting rows
                		# and columns of Gower's Center Matrix

		            # Fxn to Perform F-Stat permutation testing
            		FXN1M<-function(X1M=m2DATA)
 					{
                			# FXN to Extract data from distance matrix into vector and shuffle
                			FXN2M<-function(X2M=X1M)
						{    
						hxh<-as.data.frame(X2M) 
						nxn<-nrow(hxh)

						#shuffling predictor of interest
						P_varP<-as.data.frame(hxh[sample(1:nxn, replace=FALSE),])

						# generating decomposition matrix with permuted predictor
						P_decomposition_matrix<-cbind(INTRCPT, P_varP)  	# Combining intercept + predictor variable into single file    

	        				# Assessing predictor file for missing variables (determining missing subjects)
						P_mmm<-(1:nrow(P_decomposition_matrix))[grep("FALSE",(complete.cases(P_decomposition_matrix)))]

						# Removing subjects with missing values from both predictor and dmatrix file
						if(length(P_mmm)==0) P_decomposition_matrix<-P_decomposition_matrix 	else P_decomposition_matrix<-as.matrix(P_decomposition_matrix)[-P_mmm,]
						if(length(mmm)==0) P_gower_matrix_F<-gower_matrix 	else P_gower_matrix_F<-as.matrix(gower_matrix)[-P_mmm, -P_mmm]

						# Full Model
						########################## STEP2: Calculating Single Regression ##########################
						########################## QR Decomposition to calculate statistics ######################
    
						P_qr_docompX<-qr(P_decomposition_matrix, tol=1e-7, LAPACK=FALSE)   # Performing QR Decomposition
						P_qhat_matrix<-qr.Q(P_qr_docompX)                          		 # Extracting Q of QR decomposition for analysis
						P_hat_matrix<-P_qhat_matrix%*%(t(P_qhat_matrix))                	 # Calculate hat matrix using the decomposition
			        		P_id_minus_hat<-(diag(ncol(P_hat_matrix))-P_hat_matrix)            # Calculate (ID matrix - hat_matrix) to be used in calc of F stat

						P_SS_F<-sum(diag(((P_hat_matrix%*%P_gower_matrix_F)%*%P_hat_matrix)))   # Calculating SS(Trace)

				            ######################### Calculate (psuedo) F-Statistic ##########################
				            P_fstat_a<-P_SS_F-SS_R
				            P_fstat_b<-(sum(diag((P_id_minus_hat%*%P_gower_matrix_F)%*%P_id_minus_hat)))/(p2-p3)
				            P_fstat<-((P_fstat_a/(p3-p3_R))/P_fstat_b)
			     			}
                			abc<-replicate(xPERMS, FXN2M(X2M=X1M))
	                		perm_P<-length(which(abc>fstat))/length(abc)
                			return(perm_P)
                			}
            		pval<-FXN1M(X1M=m2DATA)
            		Nobs<-nrow(gower_matrix_F)
            		results<-c(Nobs, xPERMS, fstat_a, fstat, pval, PVE, CPVE)
            		return(results)              
				}
			}
    		if(FXN=="PVE")
        		{
        		results<-sapply(m1DATA,MDMR2)
        		return(results)
        		}
    		if(FXN=="MDMR.FXN")
        		{
        		results<-MDMR2(m1DATA)
        		return(results)
        		}
    		if(FXN=="MDMR.MULTI")
        		{
        		results<-MDMR2(m1DATA)
        		return(results)
        		}
    		}


		# MDMR TESTS: Unadjusted Variable Analysis (MARGINAL) Results
		MARGINAL.MDMR<-function(dataX=pred.vars, zPERMS=100)
			{
		    	interZ<-interX            # Defining intercept vector w/ value = "1" (used below)
			FXN.MARGINAL<-function(x, z)
				{
        			intercept<-cbind(interZ)  
				zzz<-MDMR1(m1DATA=x, INTRCPT=intercept, FXN="MDMR.FXN", xPERMS=zPERMS)
        			return(zzz)
        			}
    			aa<-length(dataX)-1
    			RESULTX<-mapply(FXN.MARGINAL, dataX, 0:aa)
    			RESULT.LABEL<-c("NOBS", "NPERMS", "SS(TRACE)", "FSTAT", "PVAL", "PVE")
    			RESULTZ<-as.data.frame(t(RESULTX))
    			RESULTZ<-RESULTZ[-7]
    			names(RESULTZ)<-RESULT.LABEL
			if(RESULT.SORT==TRUE)				# Option to Sort Reported P-Values
	        		{
    				RESULTS<-DF.SORT(RESULTZ, "FSTAT")  # Reporting Sorted Results by P-VALUE
				}
			else
				{
				RESULTS<-RESULTZ				# Reporting Unsorted Results
				}
			return(RESULTS)			      
			}


		# MDMR TESTS: Conditional Forward Regression - Variable enter model based on PVE (highest to lowest)
		FORWARD.COND.MDMR<-function(dataX=pred.vars, zPERMS=100)
    			{
    			npass<-PVE.SORT(dataX=pred.vars)
    			interZ<-interX            # Defining intercept vector w/ value = "1" (used below)
    			aa<-length(dataX)
    			RESULTZ<-matrix(-99.99,aa, 7)
    			for (i in 1:aa)
        			{
        			z<-i-1
        			varS<-npass[1:z]
        			x<-npass[i]
				if(z<1) intercept<-cbind(interZ)  else intercept<-cbind(interZ,as.data.frame(varS))
        			RESULTZ[i,]<-MDMR1(m1DATA=x, INTRCPT=intercept, FXN="MDMR.MULTI", xPERMS=zPERMS)
        			vv<-RESULTZ[i,3]
        			if(vv<.0)
            			{
            			RESULTZ[i,3]=-99.99
            			print("WARNING: Variables That Did Not Improve Model Fit were Dropped From FORWARD.REG Analysis")
            			break
            			}
        			}
			RESULTZ<-as.data.frame(RESULTZ)
    			RESULTZ<-subset(RESULTZ, V3 !=-99.99)
    			nx<-nrow(RESULTZ)
    			rownames(RESULTZ)<-names(npass[1:nx])
    			names(RESULTZ)<-c("NOBS", "NPERMS", "SS(TRACE)", "FSTAT", "PVAL", "PVE", "C.PVE")
    			return(RESULTZ)
    			}


		# MDMR TESTS: Multivariable each Variable adjusted for all other variables in model
		MULTIVARIABLE.MDMR<-function(dataX=pred.vars, zPERMS=100)
 			{
    			interZ<-interX            # Defining intercept vector w/ value = "1" (used below)
    			fxn987<-function(x, z)
        			{
        			varS<-dataX[-z]
        			intercept<-cbind(interZ,as.data.frame(varS))
        			zzz<-MDMR1(m1DATA=x, INTRCPT=intercept, FXN="MDMR.MULTI", xPERMS=zPERMS)
        			return(zzz)
        			}
    			aa<-length(dataX)
			RESULTX<-mapply(fxn987, dataX, 1:aa)
    			RESULTZ<-as.data.frame(t(RESULTX))
    			RESULTZ<-RESULTZ[-7]
    			names(RESULTZ)<-c("NOBS", "NPERMS", "SS(TRACE)", "FSTAT", "PVAL", "PVE")
			if(RESULT.SORT==TRUE)				# Option to Sort Reported P-Values
	        		{
    				RESULTS<-DF.SORT(RESULTZ, "FSTAT")  # Reporting Sorted Results by P-VALUE 
				}
			else
				{
				RESULTS<-RESULTZ				# Reporting Unsorted Results
				}
			return(RESULTS)			      
			}


		# MDMR TESTS: Custom Predictor Model - Variable adjusted for all other variables in model
		CUSTOM.MDMR<-function(dataX=pred.vars, zPERMS=100)
			{
			CUSTOM1<-function(x)
				{
				intercept<-interZ
				zzz<-MDMR1(m1DATA=x, INTRCPT=interZ, FXN="MDMR.FXN", xPERMS=zPERMS)
	      		return(zzz)
				}

			CUSTOM2<-function(x)
				{
				intercept<-interZ
				zzz<-MDMR1(m1DATA=x, INTRCPT=interZ, FXN="MDMR.MULTI", xPERMS=zPERMS)
				return(zzz)
				}

			CUSTOM3<-function(xz)			# no covariates
				{
				xix<-paste(xz,USER.INTERACTION,sep=":")
	
				if(is.null(USER.COVARS)=="TRUE")			# No Covariates
					{
					abc<-as.formula(paste("y ~ ", paste(xz, USER.INTERACTION, xix, sep="+")))
					}

				if(is.null(USER.COVARS)=="FALSE")			# 1b.  Running Custom Analysis for Model with Covariates
					{

					# checking if interaction term specified in covariates list (removed if found)
					testX<-match(USER.INTERACTION, USER.COVARS)	

					if((NA %in% testX)=="TRUE") # Interaction Variable NOT in covariates list
						{
						abc<-as.formula(paste("y ~ ", paste(xz, USER.INTERACTION, USER.COVARS, xix, sep="+")))
						}
					if((NA %in% FALSE)=="TRUE") # Interaction Variable in covariates list
						{
						abc<-as.formula(paste("y ~ ", paste(xz, USER.INTERACTION, USER.COVARS, xix, sep="+")))
						}
					}

				i.vars<-DUMMY.CODING(xfrmla=abc, data=DATA1)
				     
				zxz<-match(xix, names(i.vars))
				
				interact.term<-i.vars[zxz]					# defining variable for interaction term (full model)

				intercept<-cbind(interX,as.data.frame(i.vars[-zxz]))	# defining variables for intercept (reduced model)

				ziz<-MDMR1(m1DATA=interact.term, INTRCPT=intercept, FXN="MDMR.MULTI", xPERMS=zPERMS)
				return(ziz)
				}


			if(is.null(USER.PREDICT)=="TRUE")	# If no variables specified in USER.PREDICT all variables in in predictor file used
				{
				dataZ<-dataX
				}
			if(is.null(USER.PREDICT)=="FALSE")	# Selecting set of variables to use from predictor file
				{
				dataZ<-dataX[USER.PREDICT]
				}

			if(is.null(USER.COVARS)=="TRUE")	# No covariates specified, intercept parameter defined as vector of "1's"
				{
				interZ<-interX
				}

			if(is.null(USER.COVARS)=="FALSE")	# Covariates specified, selecting variables from predictor file + vector of "1's" (intercept)
				{
				covars<-as.data.frame(dataX[USER.COVARS])
				interZ<-cbind(interX,covars)
				}


			if(is.null(USER.INTERACTION)=="TRUE")			# 1.  Running Custom Analysis for Model w/o Interaction Term
				{
				if(is.null(USER.COVARS)=="TRUE")			# 1a. Running Custom Analysis for Model w/o Covariates
					{
					RESULTX<-sapply(dataZ, CUSTOM1)
					RESULT.LABEL<-c("NOBS", "NPERMS", "SS(TRACE)", "FSTAT", "PVAL", "PVE", "COVARIATES")
					RESULTZ<-as.data.frame(t(RESULTX))
					RESULTZ<-RESULTZ[-7]
					RESULTZ<-cbind(RESULTZ,"None")
					names(RESULTZ)<-RESULT.LABEL
					return(RESULTZ)	
					}
				if(is.null(USER.COVARS)=="FALSE")			# 1b.  Running Custom Analysis for Model with Covariates
					{
					RESULTX<-sapply(dataZ, CUSTOM2)
					RESULT.LABEL<-c("NOBS", "NPERMS", "SS(TRACE)", "FSTAT", "PVAL", "PVE", "COVARIATES")
					RESULTZ<-as.data.frame(t(RESULTX))
					RESULTZ<-RESULTZ[-7]
					RESULTZ<-cbind(RESULTZ,"See Model")
					names(RESULTZ)<-RESULT.LABEL
					}

				}
			if(is.null(USER.INTERACTION)=="FALSE")			# 2.  Running Custom Analysis for Model with Interaction Term
				{
				znames<-names(dataZ)
				RESULTX<-sapply(znames, CUSTOM3)
				RESULT.LABEL<-c("MODEL", "NOBS", "NPERMS", "SS(TRACE)", "FSTAT", "PVAL", "PVE", "INTERACTION_VAR", "COVARIATES")
				RESULTZ<-as.data.frame(t(RESULTX))
				RESULTZ<-RESULTZ[-7]
				RESULTZ<-cbind("INTERACTION", RESULTZ, USER.INTERACTION, "See Model")
				names(RESULTZ)<-RESULT.LABEL
				return(RESULTZ)	
				}			

				if(RESULT.SORT==TRUE)				# Option to Sort Reported P-Values
		        		{
    					RESULTS<-DF.SORT(RESULTZ, "FSTAT")  # Reporting Sorted Results by P-VALUE
					}
				else
					{
					RESULTS<-RESULTZ				# Reporting Unsorted Results 
					}
				return(RESULTS)
				
			}

		# UPPERCASE - Custom User specfiied Predictor Variables, specifies variables as upper case
		if(is.null(USER.PREDICT)=="FALSE") 	USER.PREDICT<-toupper(USER.PREDICT) else USER.PREDICT<-USER.PREDICT

		# UPPERCASE - Custom User specfiied Covariate Variables, specifies variables as upper case
		if(is.null(USER.COVARS)=="FALSE") USER.COVARS<-toupper(USER.COVARS) else USER.COVARS<-USER.COVARS

		# UPPERCASE - Custom User specfiied Interaction Variable, specifies variables as upper case
		if(is.null(USER.INTERACTION)=="FALSE") USER.INTERACTION<-toupper(USER.INTERACTION) else USER.INTERACTION<-USER.INTERACTION

		DATA1<-DATA
		# UPPERCASE - All Predictor Variable Names, specifies variables as upper case
		names(DATA1)<-toupper(names(DATA1)) 

		nnn<-(nrow(as.data.frame(DATA)))	# Count of subjects in predictors (used for intercept term), -1 for cat definition
		interY<-matrix(1, nnn)           	# Defining intercept vector w/ value = "1" 
		
		interX<-interY	

		# Calc Gower centered Matrix: Needed to calculate f-stat, p-value and Proportion Variance Explained (PVE)
		gower_matrix<-GOWER(DMATRIX)

		ZDATA<-names(DATA1) 
		pred.vars<-DUMMY.CODING(xfrmla=ZDATA, data=DATA1)

		RESULTS<-list(NULL,NULL,NULL,NULL,NULL)

		f.label<-paste("User Formula: ", FORMULA, sep="")

		names(RESULTS)<-c("USER.FORMULA", "UNIVARIATE", "FORWARD.REG", "MULTIVARIABLE", "USER.MODEL")

		if(class(FORMULA) == "formula")
			{
			f.vars<-DUMMY.CODING(xfrmla=FORMULA, data=DATA)
			if(length(names(f.vars)) == 1)
				{
				results.formula<-MARGINAL.MDMR(f.vars, zPERMS=PERMS)
				}
			if(length(names(f.vars)) > 1)
				{
				results.formula<-MULTIVARIABLE.MDMR(f.vars, zPERMS=PERMS)
				}
			RESULTS[1]<-list(results.formula)
			}

		if(UNIVARIATE==TRUE)
	    		{
			results.marginal<-MARGINAL.MDMR(pred.vars, zPERMS=PERMS)
			RESULTS[2]<-list(results.marginal)
			}
    
		if(FORWARD.REG==TRUE)
			{
			results.forward.conditional<-FORWARD.COND.MDMR(pred.vars, zPERMS=PERMS)
			RESULTS[3]<-list(results.forward.conditional)
			}

		if(MULTIVARIABLE==TRUE)
			{
			results.multiadjusted<-MULTIVARIABLE.MDMR(pred.vars, zPERMS=PERMS)
			RESULTS[4]<-list(results.multiadjusted)
			}

		if(USER.CUSTOM==TRUE)
			{
			if(is.null(USER.PREDICT)=="FALSE")
				{
				test2<-match(USER.PREDICT, names(DATA1))	
				if((NA %in% test2)=="TRUE") stop("WARNING: Variable(s) listed in 'USER.PREDICT' not found in PREDICT file, CHECK VARIABLE NAMES")
				}
			if(is.null(USER.COVARS)=="FALSE")
				{
				test3<-match(USER.COVARS, names(DATA1))
				if((NA %in% test3)=="TRUE") stop("WARNING: Variable(s) listed in 'USER.COVARS' not found in PREDICT file, CHECK VARIABLE NAMES")
				}
			if(is.null(USER.INTERACTION)=="FALSE")
				{
				test4<-match(USER.INTERACTION, names(DATA1))
				if((NA %in% test4)=="TRUE") stop("WARNING: Variable listed in 'USER.INTERACTION' not found in PREDICT file, CHECK VARIABLE NAMES")
				}
			results.custom<-CUSTOM.MDMR(pred.vars, zPERMS=PERMS)
	 	   	RESULTS[5]<-list(results.custom)
			}
		return(RESULTS)
		}

# v 1.Z.0

# names(predict)

# BELOW IS THE CODE FOR RUNNING THE ANALYSIS
#beta version
# system.time(aaa<-MDMR(FORMULA=DMATRIX~Age + TAR13 + Sex + TAR13:Age , DMATRIX=dmatrix, DATA=predict, PERMS=1, 
#	UNIVARIATE=TRUE, FORWARD.REG=TRUE, MULTIVARIABLE=TRUE, USER.CUSTOM=TRUE,
#	USER.PREDICT= c("tar13"), USER.COVARS= c("Sex"), USER.INTERACTION= c("Age") ))
# aaa





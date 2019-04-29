#######################
## Clean out objects ##
#######################
rm(list=ls(all=TRUE))

#################
## Source code ##
#################

source("qw_main_functions.R")

####################################################
## Function to generate a variable x.star that is ##
## correlated with x and has correlation rho      ##
####################################################

## x : original variable that x.star will be correlated with
## rho : correlation between x.star and x

x.correlated <- function(x,rho){
  ## standardize x
  x.std <- ( x - mean(x) ) / sd(x)
  
  ## coefficient needed to make x.star and x correlated
  a <- sqrt ( rho^2 / ( 1 - rho^2 ) )
  
  ## generate x.star
  x.star <- a * x.std +  rnorm(length(x),sd=1)
  
  return(x.star)
}

###############################
## Function to simulate data ##
###############################

## coef  : true coefficient values
## m     : numer of microbes
## n.i   : number of observations per diet group
## g     : number of diet groups
## x.star.correlated : indicator that we generate an extra variable that is correlated with
##   another x, but does not act on y
## rho : correlation between X[,1] and X[,1].star

data.generated <- function(coef,g,n.i,m,rho=0.8,x.star.correlated=FALSE){
  ## Total sample size
  N <- g * n.i	
  
  diet    <-  rep(c(-1,1),each=n.i)
  d.on.X  <-  c(runif(m/4*3,0.25,0.5),rep(0,m/4))
  X       <-  matrix(runif(N*m),ncol=m) 
  for (i in 1:m){
    X[,i] = X[,i] + matrix(diet*d.on.X[i],ncol=1)
  }
  
  ## Generate X[,1]^* that is correlated with X[,1] ? 
  if(x.star.correlated ==TRUE){
    X.star <- x.correlated(X[,1],rho)
  }
  
  y <- coef[1] * diet + coef[2] * X[,1] + coef[3] * X[,2] + coef[4] * X[,3] + coef[5]* X[,m] + rnorm(N,sd=1/sqrt(2))
  
  ## raw design matrix X
  if(x.star.correlated==TRUE){
    X  <-  cbind(diet,X,X.star)	
    colnames(X) <- c("Diet",paste("X_",seq(1,m),sep=""),"X_1.star")
  } else {
    X  <-  cbind(diet,X)	
    colnames(X) <- c("Diet",paste("X_",seq(1,m),sep=""))	
  }
  
  microbes <- t(X)
  phenotypes <- as.data.frame(t(y))
  rownames(phenotypes) <- "phenotype"
  
  list(microbes=microbes,phenotypes=phenotypes)
}

#####################################
## Function to do simulation study ##
#####################################

## nsimu : number of simulations
## coef  : true coefficient values
## alpha : cut-off for q-values (thresholding)
## alpha.range : different alpha values for thresholding
## delta.range : different delta values to be used in penalized loss function
## vfold : cross-validation (vfold=10) (we call this K-fold cross-validation in paper)
## m     : numer of microbes
## n.i   : number of observations per diet group
## g     : number of diet groups
## seed  : seed used for simulations
## ttest : T/F for using ttest in computing p-values for correlations
##         and partial correlations
## method : either bootstrap or smoother for computing q-values
## diet.wt : weight on diet variable
## x.star.correlated : indicator that we generate an extra variable that is correlated with
##   another x, but does not act on y
## lasso.mincp: if TRUE, we run Lasso with criterion of Mn(delta,p)
## lasso.delta.cv.mult : If TRUE, we run LASSO with delta-cv criterion multiple times (ncv times)
## rho : correlation between x_1 and x_1^*
## ncv: number of repeated 10-fold cross-validations

simu.study<- function(rho=0.8,coef=c(1.5,1,-1,-1,1),nsimu=1000,alpha=0.15,delta=2,vfold=10,
                       alpha.range=c(.01,.025,.05,.10,.15,.20),
                       delta.range=c(.25,.5,.75,1,1.5,2),
                       alpha.bh.range =c(0.01,0.05,0.1,0.15,0.2,0.5),
                       ttest=TRUE,method=c("bootstrap","smoother")[1],
                       m=40,n.i=20,g=2,seed=1110,diet.wt=100,thresh.q=FALSE,
                       pi0.true=FALSE,pi0.val=0.9,x.star.correlated=FALSE, lasso.mincp=TRUE,
			     percents.range=c(50,60,70,80,90,100),lasso.delta.cv.mult=FALSE, ncv=100){
  
######################################
## Form matrices for storing output ##
######################################
  
  ## To get dimensions and rownames for output, we run a simulated data set
  data <- data.generated(coef,g,n.i,m,rho,x.star.correlated)
  microbes <- data$microbes
  
  out.nrow <- nrow(microbes)
  out.rownames <- rownames(microbes)	

#############################################
## Output for Benjamini-Hochberg Correction ##
#############################################
  
  ## out.benhoch : stores results from applying Benjamini-Hochberg correction
  out.benhoch <- as.data.frame(matrix(0, nrow = out.nrow, ncol = length(alpha.bh.range),
                                 dimnames = list(out.rownames,paste("benhoch.alpha.",alpha.bh.range,sep="") )))
  
  ## out.benhoch.pval.adjust: store adjusted p-values from BH correction
  out.benhoch.pval.adjust <-  as.data.frame(matrix(0, nrow = out.nrow, ncol = nsimu,
                                                   dimnames = list(out.rownames, seq(1,nsimu))))

  ## out.benhoch.cor: stores results from applying Benjamini-Hochberg correction to p-values that ignore diet
  out.benhoch.cor <- as.data.frame(matrix(0, nrow = out.nrow, ncol = length(alpha),
                                          dimnames = list(out.rownames,paste("benhoch.cor.alpha.",alpha,sep=""))))
  
  
#####################################
## Output for thresholding methods ##
#####################################
  
  ## out.cor  : stores thresholding results for testing if a microbe has an effect on phenotype, but NOT 
  ##            accounting for diet
  ## That is, we test : H_0 : \beta_{x_j}=0
  
  out.cor <- as.data.frame(matrix(0, nrow = out.nrow, ncol = length(alpha),
                                  dimnames = list(out.rownames,paste("cor.alpha.",alpha,sep=""))))
  
  
  ## out.parcor  : stores thresholding results for testing if a microbe has an effect on phenotype, but AFTER
  ##            accounting for diet
  ## That is, we test : H_0 : \beta_{x_j|z}=0
  
  out.parcor <- as.data.frame(matrix(0, nrow = out.nrow, ncol = length(alpha.range),
                                     dimnames = list(out.rownames,paste("parcor.alpha.",alpha.range,sep="") )))
  
  ## out.qvalue : stores q-values from out.parcor results
  out.qvalue <- as.data.frame(matrix(0, nrow = out.nrow, ncol = nsimu,
                                     dimnames = list(out.rownames, seq(1,nsimu))))
  
  ## out.pvalue : stores p-values from out.parcor results
  out.pvalue <- as.data.frame(matrix(0, nrow = out.nrow, ncol = nsimu,
                                     dimnames = list(out.rownames, seq(1,nsimu))))
  
##########################################################
## Output for weighted Lasso with Mn(delta,p) criterion ##
##########################################################
  
  ## Weight functions
  g1 <- function(x){
    return(x)
  }
  
  g2 <- function(x){
    return(sqrt(x))
  }
  
  g3 <- function(x){
    return(1/abs(x))
  }
  
  g4 <- function(x){
    return(x^2)
  }
  
  ## out.w1 : stores results from weighted lasso when weights are set to 1, and
  ##          weight function g1
  
  out.w1 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                 dimnames = list(out.rownames,paste("w1.delta.",delta,sep=""))))
  
  ## out.w2 : stores results from weighted lasso when weights are set to q-values BEFORE taking into account diet,
  ##           and weight function g2
  
  out.w2 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                 dimnames = list(out.rownames,paste("w2.delta.",delta,sep=""))))
  
  ## out.w3 : stores results from weighted lasso when weights are set to q-values BEFORE taking into account diet,
  ##           and weight function g1
  
  out.w3 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                 dimnames = list(out.rownames,paste("w3.delta.",delta,sep=""))))
  
  ## out.w4 : stores results from weighted lasso when weights are set to q-values AFTER taking into account diet,
  ##           and weight function g2
  
  out.w4 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta.range),
                                 dimnames = list(out.rownames,paste("w4.delta.",delta.range,sep=""))))

  ## out.w5 : stores results from weighted lasso when weights are set to q-values AFTER taking into account diet,
  ##           and weight function g1
  
  out.w5 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta.range),
                                 dimnames = list(out.rownames,paste("w5.delta.",delta.range,sep=""))))
  
  ## out.w6 : stores results from weighted lasso when weights absolute value of partial correlations,
  ##           and weight function g3
  
  out.w6 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta.range),
                                 dimnames = list(out.rownames,paste("w6.delta.",delta.range,sep=""))))
  
  ## out.ndw1 : stores results from weighted lasso when weights are set to 1, and
  ##          weight function g1, diet NOT forced in the model
  
  out.ndw1 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                   dimnames = list(out.rownames,paste("ndw1.delta.",delta,sep=""))))
  
  ## out.ndw2 : stores results from weighted lasso when weights are set to q-values BEFORE taking into account diet,
  ##           and weight function g2, diet NOT forced in the model
  
  
  out.ndw2 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                   dimnames = list(out.rownames,paste("ndw2.delta.",delta,sep=""))))
  
  ## out.w3 : stores results from weighted lasso when weights are set to q-values BEFORE taking into account diet,
  ##           and weight function g1, diet NOT forced in the model
  
  out.ndw3 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                   dimnames = list(out.rownames,paste("ndw3.delta.",delta,sep=""))))
  
  ## out.ndw4 : stores results from weighted lasso when weights are set to q-values AFTER taking into account diet,
  ##           and weight function g2, diet NOT forced in the model
  
  
  out.ndw4 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                   dimnames = list(out.rownames,paste("ndw4.delta.",delta,sep=""))))

  ## out.ndw5 : stores results from weighted lasso when weights are set to q-values AFTER taking into account diet,
  ##           and weight function g1, diet NOT forced in the model
  
  
  out.ndw5 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta),
                                   dimnames = list(out.rownames,paste("ndw5.delta.",delta,sep=""))))

#######################################################################################################
## Output for weighted Lasso with M_n(delta,p) criterion but delta chosen via cross-validation x 100 ##
#######################################################################################################
  

  ## mult.cv.delta.out.w5 : stores results from weighted lasso when weights are set to q-values AFTER taking into account diet,
  ##           and weight function g1
  
  mult.cv.delta.out.w5 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=nsimu,
                                         dimnames = list(out.rownames,paste("w5.mult.nsimu.",seq(1,nsimu),sep=""))))

  mult.delta.w5 <- as.data.frame(matrix(0, nrow = 1, ncol = ncv,
                                  dimnames = list("delta", seq(1,ncv))))

  mult.cv.delta.out.w5.summary <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(percents.range),
                                          dimnames = list(out.rownames,paste("w5.mult.cv.",percents.range,sep=""))))


  ## mult.cv.delta.out.w6 : stores results from weighted lasso when weights absolute value of partial correlations,
  ##           and weight function g3
  
  mult.cv.delta.out.w6 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=nsimu,
                                         dimnames = list(out.rownames,paste("w6.mult.nsimu.",seq(1,nsimu),sep=""))))

  mult.delta.w6 <- as.data.frame(matrix(0, nrow = 1, ncol = ncv,
                                  dimnames = list("delta", seq(1,ncv))))


  mult.cv.delta.out.w6.summary <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(percents.range),
                                          dimnames = list(out.rownames,paste("w6.mult.cv.",percents.range,sep=""))))


#######################
## Start Simulations ##
#######################
  
  ## set seed
  set.seed(seed)
  
  for(j in 1:nsimu) {
#################
## Generate data #
#################	
    data <- data.generated(coef,g,n.i,m,rho,x.star.correlated)
    microbes <- data$microbes
    phenotypes <- data$phenotypes
    
    
###########################################
    ## Get correlation and partial correlation #
############################################
    cor.out <- correlations(microbes,phenotypes,partial=FALSE,ttest=ttest)
    parcor.out <- correlations(microbes,phenotypes,partial=TRUE,ttest=ttest)
    
################
    ## Get q-values #
################

    ## Results for testing if a microbe has an effect on phenotype, but NOT 
    ##            accounting for diet
    ## That is, we test : H_0 : \beta_{x_j}=0
    
    microbe.cor.out.qvalues <- q.computations(cor.out, method=method,
                                              plots=FALSE,file="cor",
                                              pi0.true=pi0.true,pi0.val=pi0.val)
    
    microbe.cor.out <- q.interest(microbe.cor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
    out.cor   <- out.cor  + c(0,t(microbe.cor.out$interest))

    ## Results from Benjamini-Hochberg adjusted p-values when p-values do not account for diet
    benhoch.cor.results <- ben.hoch.interest(cor.out$pvalues,alpha=alpha)
    out.benhoch.cor <- out.benhoch.cor + c(0,t(benhoch.cor.results$interest))
    
    
    ## Results for testing if a microbe has an effect on phenotype, but AFTER
    ##            accounting for diet
    ## That is, we test : H_0 : \beta_{x_j|z}=0
    
    microbe.parcor.out.qvalues <- q.computations(parcor.out,method=method,
                                                 plots=FALSE,file="parcor",
                                                 pi0.true=pi0.true,pi0.val=pi0.val)
    
    out.qvalue[,j] <- c(0,t(microbe.parcor.out.qvalues$qval.mat))
    out.pvalue[,j] <- c(0,t(parcor.out$pvalues))
    
    
    for(a in 1:length(alpha.range)){
      q.out <- q.interest(microbe.parcor.out.qvalues$qval.mat,alpha=alpha.range[a],criteria="less")
      out.parcor[,a] <- out.parcor[,a] + c(1,t(q.out$interest))
    }
    
#####################################################################
    ## Results for Benjamini-Hochberg Method ##
#####################################################################
    
    for(b in 1:length(alpha.bh.range)){
      benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha.bh.range[b])
      out.benhoch[,b] <- out.benhoch[,b] + c(1,t(benhoch.results$interest))
    }
    out.benhoch.pval.adjust[,j] <- c(0,t(benhoch.results$pval.adjust))
    
    
    
###############################################
    ## Weighted Lasso with Mn(delta,p) criterion ##
###############################################	
    
    if(lasso.mincp==TRUE){	
###################################################
      ## Lasso calculations: With Diet forced in model ##
###################################################
      
      ##		g1 <- function(x){
      ##			return(x)
      ##		}
      
      ##		g2 <- function(x){
      ##			return(sqrt(x))
      ##		}
      
      ##		g3 <- function(x){
      ##			return(1/abs(x))
      ##		}
      
      ##		g4 <- function(x){
      ##			return(x^2)
      ##		}
      
      include.diet <- "TRUE"
      
      ## Weights set to 1
      weights <- matrix(1,nrow=nrow(microbes)-1,ncol=nrow(phenotypes))
      lasso.w1 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,file="weight1_",
                                     include.diet=include.diet,diet.wt=diet.wt)
      out.w1 <- out.w1 + lasso.w1$interest
      
      ## Weights set to q-values BEFORE taking into account diet
      weights <- microbe.cor.out.qvalues$qval.mat

      lasso.w2 <- lasso.computations(weights,microbes,phenotypes,g2,plots=FALSE,file="weight2_",
                                     include.diet=include.diet,diet.wt=diet.wt,thresh.q=thresh.q)
      out.w2 <- out.w2 + lasso.w2$interest
      
      
      lasso.w3 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,file="weight3_",
                                     include.diet=include.diet,diet.wt=diet.wt,thresh.q=thresh.q)
      out.w3 <- out.w3 + lasso.w3$interest
      
      
      ## Weights set to q-values after taking into account diet
      weights <- microbe.parcor.out.qvalues$qval.mat
      
      for(d in 1:length(delta.range)){
        lasso.w4 <- lasso.computations(weights,microbes,phenotypes,g2,plots=FALSE,file="weight4_",
                                       include.diet=include.diet,diet.wt=diet.wt,thresh.q=thresh.q,
                                       delta=delta.range[d])
        out.w4[,d] <- out.w4[,d] + lasso.w4$interest
      }
      
      for(d in 1:length(delta.range)){
        lasso.w5 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,file="weight5_",
                                       include.diet=include.diet,diet.wt=diet.wt,thresh.q=thresh.q,delta=delta.range[d])
        out.w5[,d] <- out.w5[,d] + lasso.w5$interest
      }
    
      
      ## Weights set to absolute value of partial correlations
      weights <- parcor.out$estimate
      
      for(d in 1:length(delta.range)){
        lasso.w6 <- lasso.computations(weights,microbes,phenotypes,g3,plots=FALSE,file="weight6_",
                                       include.diet=include.diet,diet.wt=diet.wt,corr.g=TRUE,delta=delta.range[d])
        out.w6[,d] <- out.w6[,d] + lasso.w6$interest
      }
      
      
#####################################################
      ## Lasso calculations: With Diet NOT forced in model #
#####################################################
      
      include.diet <- "FALSE"
      
      weights <- matrix(1,nrow=nrow(microbes)-1,ncol=nrow(phenotypes))
      ndlasso.w1 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,file="ndweight1_",
                                       include.diet=include.diet)
      
      out.ndw1 <- out.ndw1 + ndlasso.w1$interest
      
      ## Weights set to q-values BEFORE taking into account diet
	weights <- microbe.cor.out.qvalues$qval.mat
      
      ndlasso.w2 <- lasso.computations(weights,microbes,phenotypes,g2,plots=FALSE,
                                       file="ndweight2_",include.diet=include.diet,thresh.q=thresh.q)
      
      out.ndw2 <- out.ndw2 + ndlasso.w2$interest
      
      ndlasso.w3 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                       file="ndweight3_",include.diet=include.diet,thresh.q=thresh.q)
      out.ndw3 <- out.ndw3 + ndlasso.w3$interest
      

      ## Weights set to q-values after taking into account diet
      weights <- microbe.parcor.out.qvalues$qval.mat

      ndlasso.w4 <- lasso.computations(weights,microbes,phenotypes,g2,plots=FALSE,file="ndweight4_",
                                       include.diet=include.diet,thresh.q=thresh.q)
      out.ndw4 <- out.ndw4 + ndlasso.w4$interest
      
      
      ndlasso.w5 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                       file="ndweight5_",include.diet=include.diet,thresh.q=thresh.q)
      out.ndw5 <- out.ndw5 + ndlasso.w5$interest
      
    }

############################################################################################
## Weighted Lasso with M_n(delta,p) criterion but delta chosen via cross-validation x 100 ##
############################################################################################
    
    if(lasso.delta.cv.mult==TRUE){		
###################################################
      ## Lasso calculations: With Diet forced in model ##
###################################################
      
      ##		g1 <- function(x){
      ##			return(x)
      ##		}
      
      ##		g2 <- function(x){
      ##			return(sqrt(x))
      ##		}
      
      ##		g3 <- function(x){
      ##			return(1/abs(x))
      ##		}
      
      ##		g4 <- function(x){
      ##			return(x^2)
      ##		}
      
      include.diet <- TRUE
    
      ## Weights set to q-values after taking into account diet
      weights <- microbe.parcor.out.qvalues$qval.mat
      
      for(v in 1:ncv){
        mult.cv.delta.lasso.w5 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,file="weight5_",
                                                     include.diet=include.diet,diet.wt=diet.wt,thresh.q=thresh.q,delta=delta,
                                                     cv.criterion="delta_cv",vfold=vfold)
        mult.cv.delta.out.w5[,j] <- mult.cv.delta.out.w5[,j] + mult.cv.delta.lasso.w5$interest
        mult.delta.w5[,v] <- mult.delta.w5[,v] + mult.cv.delta.lasso.w5$delta.out
      }
      
      
      ## Weights set to absolute value of partial correlations
      weights <- parcor.out$estimate
      
      for(v in 1:ncv){
        mult.cv.delta.lasso.w6 <- lasso.computations(weights,microbes,phenotypes,g3,plots=FALSE,file="weight6_",
                                                     include.diet=include.diet,diet.wt=diet.wt,corr.g=TRUE,delta=delta,
                                                     cv.criterion="delta_cv",vfold=vfold)
        mult.cv.delta.out.w6[,j] <- mult.cv.delta.out.w6[,j] + mult.cv.delta.lasso.w6$interest
        mult.delta.w6[,v] <- mult.delta.w6[,v] + mult.cv.delta.lasso.w6$delta.out
	}
  }

 }

	#############################################
	## Organize results from multiple CV-delta ##
	#############################################

	if(lasso.delta.cv.mult==TRUE){		
		for(v in 1:length(percents.range)){

###################################################
## Lasso calculations: With Diet forced in model ##
###################################################

		      ## Weights set to q-values after taking into account diet
			 mult.cv.delta.out.w5.summary[,v] <- org.mult.cv.delta(mult.cv.delta.out.w5,percents.range[v],ncv)
     
		      ## Weights set to absolute value of partial correlations
			 mult.cv.delta.out.w6.summary[,v] <- org.mult.cv.delta(mult.cv.delta.out.w6,percents.range[v],ncv)
		}
	}

  
	######################
	## Organize results ##
	######################

  out <- cbind(out.cor,out.benhoch.cor,out.parcor,out.benhoch,out.w1,out.w2,out.w3,
               out.w4,out.w5,out.w6,out.ndw1,out.ndw2,out.ndw3,out.ndw4,out.ndw5,
               mult.cv.delta.out.w5.summary, mult.cv.delta.out.w6.summary,
               mult.cv.delta.out.w5,mult.cv.delta.out.w6 )



	list(out=out,
		out.pvalue = out.pvalue,
           	out.benhoch.pval.adjust = out.benhoch.pval.adjust,
		out.qvalue = out.qvalue,
		mult.delta.w5 = mult.delta.w5,
		mult.delta.w6 = mult.delta.w6)

}

## function to organize output from multiple cv-delta process
org.mult.cv.delta <- function(mult.cv.delta.output,percent,ncv){
	make.criteria <- function(ncv,percent){
		function(x){
			return(x >= ncv * percent/100)
		}
	}

	criteria <- make.criteria(ncv,percent)
	
	check.criteria <- apply(mult.cv.delta.output,c(1,2),criteria)

	output <- apply(check.criteria,1,sum) 
	return(output)
}


# function to organize data into grouped averages

org.simu <- function(m,out.simu,nsimu,x.star.correlated=FALSE,coef){
  if(x.star.correlated==TRUE){
    num.row = 7
  } else {
    num.row=6
  }
  out <- matrix(0,nrow = num.row, ncol=ncol(out.simu))
  out <- as.data.frame(out)
  colnames(out) <- colnames(out.simu)
  
  if(x.star.correlated==TRUE){
    rownames(out) <- c("Diet","Group 1","Group 2","Group 3","Group 4","Group 5","FDR")	
  } else {
    rownames(out) <- c("Diet","Group 1","Group 2","Group 3","Group 4","FDR")	
  }

  if(coef[3]==0){
    out["Diet",] <- out.simu[1,]
    out["Group 1",] <- apply(out.simu[c(2,4),],2,mean)
    out["Group 2",] <- apply(out.simu[c(3,5:(m/4*3+1)),],2,mean)
    out["Group 3",] <- apply(out.simu[(m/4*3+2):m,],2,mean)
    out["Group 4",] <- apply(out.simu[m+1,],2,mean)
    
    if(x.star.correlated==TRUE){
      out["Group 5",] <- apply(out.simu[m+2,],2,mean)		# for X_1.star
    }
    
    out <- out * 100/nsimu
    
    R <- #1 * out["Diet",] + 
      2 * out["Group 1",] + (1+(m/4*3+1)-5+1) * out["Group 2",] + ( m - (m/4*3+2)+1 ) * out["Group 3",] + 1 * out["Group 4",] 
    
    
    F <- (1+(m/4*3+1)-5+1) * out["Group 2",] + ( m - (m/4*3+2)+1 ) * out["Group 3",] 

  } else {
    out["Diet",] <- out.simu[1,]
    out["Group 1",] <- apply(out.simu[c(2,3,4),],2,mean)
    out["Group 2",] <- apply(out.simu[5:(m/4*3+1),],2,mean)
    out["Group 3",] <- apply(out.simu[(m/4*3+2):m,],2,mean)
    out["Group 4",] <- apply(out.simu[m+1,],2,mean)
    
    if(x.star.correlated==TRUE){
      out["Group 5",] <- apply(out.simu[m+2,],2,mean)		# for X_1.star
    }
    
    out <- out * 100/nsimu
    
    R <- #1 * out["Diet",] + 
      3 * out["Group 1",] + ((m/4*3+1)-5+1) * out["Group 2",] + ( m - (m/4*3+2)+1 ) * out["Group 3",] + 1 * out["Group 4",] 
    
    
    F <- ((m/4*3+1)-5+1) * out["Group 2",] + ( m - (m/4*3+2)+1 ) * out["Group 3",] 
  }
    
    
  if(x.star.correlated==TRUE){
    R <- R + 1 * out["Group 5",]
    F <- F + 1 * out["Group 5",]
  }
  
  out["FDR",] <- F/R
  FDR <- out["FDR",]
  
  list(out=out,FDR=FDR)
  
}

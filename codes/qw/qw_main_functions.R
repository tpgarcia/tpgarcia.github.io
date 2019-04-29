##############################################################
## R code for "Structured Variable Selection with q-values" ##
##############################################################

###############
## Libraries ##
###############

library(xtable)		# to create LaTeX tables
library(lars)		# for LASSO approach
library(plotrix)		# for computing standard errors of mean in simulation study


#########################################################
## Function for making data frames for storing results ##
#########################################################

store.micropheno <- function(microbes,phenotypes){
	out <- as.data.frame(matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes)))
	rownames(out) <- rownames(microbes)
	colnames(out) <- rownames(phenotypes)
	return(out)
}

store.micro <- function(microbes){
	out <- as.data.frame(rep(0,nrow(microbes)))
	rownames(out) <- rownames(microbes)
	return(out)
}


##########################################################
## Functions to get partial correlation and its p-value #
##########################################################

# When pearson correlation is 0 (because std. deviation is 0), I am setting
# p-value  to 1.

corr.pvalue <- function(x,y,method="pearson",alternative="two.sided",ttest="FALSE"){
	x <- as.numeric(x)
	y <- as.numeric(y)

	out <- cor.test(x,y,alternative=alternative,method=method,na.action=na.omit)
	estimate <- out$estimate
	
	if(ttest=="FALSE"){
		p.value <- out$p.value
	} else {
		y1 <- y
		x1 <- x
		summary.out <- summary(lm(y1~  x1))
		p.value <- summary.out$coefficients["x1","Pr(>|t|)"]
	}
	list(p.value=p.value,estimate=estimate)
}



parcorr.pvalue <- function(x,y,z,method="pearson",alternative="two.sided",ttest="FALSE"){
	x <- as.numeric(x)
	y <- as.numeric(y)
	z <- as.numeric(z)

	xres <- residuals(lm(x~factor(z)))
	yres <- residuals(lm(y~factor(z)))
	out <- corr.pvalue(xres,yres,method,alternative,ttest="FALSE")
	estimate <- out$estimate

	if(ttest=="FALSE"){
		p.value <- out$p.value
	} else {
		y1 <- y
		x1 <- x
		summary.out <- summary(lm(y1 ~ factor(z) +  x1))
		p.value <- summary.out$coefficients["x1","Pr(>|t|)"]
	}

	list(p.value=p.value,estimate=estimate)
}



correlations <- function(microbes,phenotypes,partial="FALSE",ttest="FALSE"){

	data.phenotypes <- phenotypes
	data.microbes <- microbes[-which(rownames(microbes)=="Diet"),]
	diet <- microbes["Diet",]

	# Setting up matrices to store Pearson correlations and p-values
	correlation <- store.micropheno(data.microbes,data.phenotypes)
	pvalues <- store.micropheno(data.microbes,data.phenotypes)

	# Computing pearson correlations and p-values
	for (i in 1: nrow(data.microbes)) {
		for(j in 1:nrow(data.phenotypes)) {
			if(partial=="TRUE"){
				tmp <- parcorr.pvalue(data.microbes[i,],data.phenotypes[j,],diet,ttest=ttest)
			} else {
				tmp <- corr.pvalue(data.microbes[i,],data.phenotypes[j,],ttest=ttest)
			}
			correlation[i,j] <- tmp$estimate
			pvalues[i,j] <- tmp$p.value
		}
	}
	list(estimate = correlation, pvalues = pvalues)
}



####################
# Qvalues function #
####################

# modification of qplot so as to get individual plots

qplot2 <- function (qobj, rng = c(0, 0.1), smooth.df = 3, smooth.log.pi0 = FALSE, whichplot=1) 
{
    q2 <- qobj$qval[order(qobj$pval)]
    if (min(q2) > rng[2]) {
        rng <- c(min(q2), quantile(q2, 0.1))
    }
    p2 <- qobj$pval[order(qobj$pval)]
    #par(mfrow = c(2, 2))
    lambda <- qobj$lambda
    if (length(lambda) == 1) {
        lambda <- seq(0, max(0.9, lambda), 0.05)
    }
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
        pi0[i] <- mean(p2 > lambda[i])/(1 - lambda[i])
    }
    if (smooth.log.pi0) 
        pi0 <- log(pi0)
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    if (smooth.log.pi0) {
        pi0 <- exp(pi0)
        spi0$y <- exp(spi0$y)
    }
    pi00 <- round(qobj$pi0, 3)
	if(whichplot==1){
	#TG change, 02/07/2012 ("gamma" is lambda in manuscript)
	par(mar=c(5, 4, 4, 2)+1)
    plot(lambda, pi0, xlab = expression(gamma), ylab = expression(hat(pi)(gamma)),pch=19,cex.lab=1.5,cex.axis=1.5)
    mtext(substitute(hat(pi) == that, list(that = pi00)),cex=2)
    lines(spi0,lwd=2)
	} else if (whichplot==2){
    plot(p2[q2 >= rng[1] & q2 <= rng[2]], q2[q2 >= rng[1] & q2 <= 
        rng[2]], type = "l", xlab = "p-value", ylab = "q-value")
	} else if (whichplot==3){
    plot(q2[q2 >= rng[1] & q2 <= rng[2]], (1 + sum(q2 < rng[1])):sum(q2 <= 
        rng[2]), type = "l", xlab = "q-value cut-off", ylab = "Number of significant microbes",cex.lab=1.5,cex.axis=1.5)
	} else if (whichplot==4){
    plot((1 + sum(q2 < rng[1])):sum(q2 <= rng[2]), q2[q2 >= rng[1] & 
        q2 <= rng[2]] * (1 + sum(q2 < rng[1])):sum(q2 <= rng[2]), 
        type = "l", xlab = "significant tests", ylab = "expected false positives")
	}
    #par(mfrow = c(1, 1))
}

####################################################################
## qvalue.adj. is as used by JD Storey with some adjustments made ##
####################################################################

qvalue.adj<-function (p = NULL, lambda = seq(0, 0.9, 0.05), pi0.method = "smoother", 
    fdr.level = NULL, robust = FALSE, gui = FALSE, smooth.df = 3, 
    smooth.log.pi0 = FALSE,pi0.true=FALSE,pi0.val=0.9) 
{
    if (is.null(p)) {
        qvalue.gui()
        return("Launching point-and-click...")
    }
    if (gui & !interactive()) 
        gui = FALSE
    if (min(p) < 0 || max(p) > 1) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: p-values not in valid range.", 
                "\n"))), parent.frame())
        else print("ERROR: p-values not in valid range.")
        return(0)
    }
    if (length(lambda) > 1 && length(lambda) < 4) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.", 
                "\n"))), parent.frame())
        else print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
        return(0)
    }
    if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 
        1)) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                "\n"))), parent.frame())
        else print("ERROR: Lambda must be within [0, 1).")
        return(0)
    }
    m <- length(p)
    if (length(lambda) == 1) {
        if (lambda < 0 || lambda >= 1) {
            if (gui) 
                eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                  "\n"))), parent.frame())
            else print("ERROR: Lambda must be within [0, 1).")
            return(0)
        }
        pi0 <- mean(p >= lambda)/(1 - lambda)
        pi0 <- min(pi0, 1)
    }
    else {
      # TG added pi0.true option
      if(pi0.true==TRUE){
        pi0 <- pi0.val
      } else { 
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
          pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }
        
     # TG change: Remove any pi0 which have entry 0, and adjust lambda values
        
        if(sum(pi0==0)>0){
          ind.zero <- which(pi0==0)
          pi0 <- pi0[-ind.zero]
          lambda <- lambda[-ind.zero]
        }
        
        if (pi0.method == "smoother") {
          if (smooth.log.pi0) 
            pi0 <- log(pi0)
          spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
          pi0 <- predict(spi0, x = max(lambda))$y
          if (smooth.log.pi0) 
            pi0 <- exp(pi0)
          pi0 <- min(pi0, 1)
        }
        else if (pi0.method == "bootstrap") {
          minpi0 <- min(pi0)
          mse <- rep(0, length(lambda))
          pi0.boot <- rep(0, length(lambda))
          for (i in 1:1000) {
            p.boot <- sample(p, size = m, replace = TRUE)
            for (j in 1:length(lambda)) {
              pi0.boot[j] <- mean(p.boot > lambda[j])/(1 - 
                                                       lambda[j])
            }
            mse <- mse + (pi0.boot - minpi0)^2
          }
          pi0 <- min(pi0[mse == min(mse)])
          pi0 <- min(pi0, 1)
        }
        else {
          print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
          return(0)
        }
      }
    }

    if (pi0 <= 0) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.", 
                "\n"))), parent.frame())
        else print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
        return(0)
    }
    if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
        1)) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].", 
                "\n"))), parent.frame())
        else print("ERROR: 'fdr.level' must be within (0, 1].")
        return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
        idx <- sort.list(x)
        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)
        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl
        return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    if (robust) {
        qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
    }
    if (!is.null(fdr.level)) {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, fdr.level = fdr.level, significant = (qvalue <= 
                fdr.level), lambda = lambda)
    }
    else {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, lambda = lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}


q.computations <- function(out, method=c("smoother","bootstrap")[2],
				plots=TRUE,file="name",
				pi0.true=FALSE,pi0.val=0.9){

	qval.mat <- matrix(0,nrow=nrow(out$pvalues),ncol=ncol(out$pvalues))
	qval.mat <- as.data.frame(qval.mat)
	rownames(qval.mat) <- rownames(out$pvalues)
	colnames(qval.mat) <- colnames(out$pvalues)


	for(i in 1:ncol(qval.mat)){
		pvalues <- out$pvalues[,i]
		estimate <- out$estimate[,i]
		
		qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.95,by=0.01),robust=FALSE,pi0.true=pi0.true,pi0.val=pi0.val)
		qval <- qobj$qvalues	
		pi0 <- qobj$pi0
		qval.mat[,i] <- qval

		cnames <- colnames(out$pvalues)

		# Plots
		if(plots==TRUE){
		
			# Density histogram of p-values
			postscript(paste(file,"_histpval_",cnames[i],".eps",sep=""))
			hist(pvalues,main="",xlab="Density of observed p-values",ylab="",freq=FALSE,yaxt="n",ylim=c(0,5),cex.lab=2,cex.axis=2)
			abline(h=1,lty=1)
			abline(h=qobj$pi0,lty=2)
			axis(2,at = c(0,qobj$pi0,1,2,3,4),labels=c(0,round(qobj$pi0,2),1,2,3,4),las=2,cex.lab=2,cex.axis=2)
			box(lty=1,col="black") 	# for a box around plot
			dev.off()

			# Density histogram of estimate
			postscript(paste(file,"_histest_",cnames[i],".eps",sep=""))
			hist(estimate,main="",xlab="Density of observed statistic",ylab="",yaxt="n",freq=FALSE,ylim=c(0,5),cex.lab=2,cex.axis=2)
			axis(2,at = c(0,1,2,3,4),labels=c(0,1,2,3,4),las=2,cex.lab=2,cex.axis=2)
			dev.off()	
		
			# Plot of \hat \pi vs. \lambda
			postscript(paste(file,"_pihat_vs_lambda_",cnames[i],".eps",sep=""))
			qplot2(qobj,rng=c(0,1),whichplot=1)
			dev.off()

			# Plot of significant tests vs. q cut-off
			postscript(paste(file,"_sigtest_vs_qval_",cnames[i],".eps",sep=""))
			qplot2(qobj,rng=c(0,1),whichplot=3)
			dev.off()
		}
	}	
	list(qval.mat=qval.mat,pi0=pi0)
}





q.interest <- function(qval.mat,alpha=0.15, criteria = c("more","less")[2]){

	interest <- matrix(0,nrow=nrow(qval.mat),ncol=ncol(qval.mat))
	interest <- as.data.frame(interest)
	rownames(interest) <- rownames(qval.mat)
	colnames(interest) <- colnames(qval.mat)

	for(i in 1:ncol(interest)){
		qval <- qval.mat[,i]

		if(criteria == "less"){
			ind <- which(qval <= alpha)
		} else {
			ind <- which(qval > alpha)
		}
		if(length(ind)>0){
			interest[ind,i] <- 1
		}
	}
	list(interest=interest)
}


#############################
# Benjamini-Hochberg Method #
#############################


ben.hoch.interest <- function(pvalues,alpha=0.05){
  
  interest <- matrix(0,nrow=nrow(pvalues),ncol=ncol(pvalues))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(pvalues)
  colnames(interest) <- colnames(pvalues)
  
  pval.adjust <- matrix(0,nrow=nrow(pvalues),ncol=ncol(pvalues))
  pval.adjust <- as.data.frame(pval.adjust)
  rownames(pval.adjust) <- rownames(pval.adjust)
  colnames(pval.adjust) <- colnames(pval.adjust)
  
  for(i in 1:ncol(interest)){
    pval <- pvalues[,i]
    
    ## adjust p-values using Benjamini-Hochberg method
    pval.adjust[,i] <- p.adjust(pval,method="BH")
    
    ind <- which(pval.adjust[,i] <= alpha)
    
    if(length(ind)>0){
      interest[ind,i] <- 1
    }
  }
  list(interest=interest,pval.adjust=pval.adjust)
}



#############################
## weighted LASSO Approach ##
#############################


# function to standardize a vector x
make.std <- function(x){
	N <- length(x)
	( x-mean(x) ) / ( sd(as.vector(x)) * sqrt( N / (N-1) ) ) 
}

make.center <- function(x){
	return(x-mean(x))  
}

lasso <- function(weights,phenotype,microbe,g,file="file",plots="FALSE",include.diet="TRUE",
                  diet.wt=1000,thresh.q=FALSE,corr.g=FALSE,delta=2,
                  cv.criterion="delta_cv",vfold=10){
  y <- t(phenotype)
  X <- t(microbe)
  
  N <- length(y)
  
  ## Standardize variables
  y1 <- make.std(y)	
  X1 <- apply(X,2,make.std)
  
  ## To maintain stability, remove excessively small weights
  if(thresh.q==TRUE){
    for(i in 1:length(weights)){
      weights[i] <- max(0.0001,abs(weights)[i])
    }
  }
  
  
  ## Weights used
  if(include.diet==TRUE){
    if(corr.g==FALSE){
      wts <- c(min(abs(weights))/diet.wt,weights)
    } else {
      wts <- c(diet.wt/min(abs(weights)),weights)
    }
  } else {
    if(corr.g==FALSE){
      wts <- c(1,weights)
    } else {
      wts <- c(0.5,weights)
    }
  }
  X1 <- X1 %*% diag(1/g(wts))
  
  
  ## Fix use.Gram for Lasso
  if(ncol(X1)>500){
    use.Gram <- FALSE
  } else {
    use.Gram <- TRUE
  }
  
  ## Run Lasso
  wLasso.out <-  lars(X1, y1, type = c("lasso"), 
                      trace = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = use.Gram)
  
  ## entry in Lasso
  entry.variables <- as.numeric(wLasso.out$entry)
  
  
  if(cv.criterion=="delta_cv"){

    cv.out <- cv.delta(y1,X1,K=vfold)
    predict.out <- cv.out$predict.out
    delta.out <- cv.out$delta
    
  } else {
    ## use Cp-like criterion to find best descriptive model
    p = dim(X1)[2]
    s = length(wLasso.out$df)
    p.pos = NULL
   
    RSS = NULL
    for (i in 1:s){
      RSS[i] = sum((y1-predict(wLasso.out, X1, s=i, type = c("fit"))$fit)**2)
      p.pre = predict(wLasso.out, X1, s=i, type = c("coefficients"))$coefficients
      p.pos = c(p.pos,length(p.pre[abs(p.pre)>0]))
    }
    
    ## Get estimated MSE
    MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
    MSE <- MSE^2
    
    p.min = which.min(RSS/MSE+delta*p.pos)

    ## final best descriptive model
    predict.out <- predict(wLasso.out, X1, s=p.min, type = c("coefficients"))
    delta.out <- delta

  }
  
  ind <- which(abs(predict.out$coefficients)>1e-10)
  sig.variables <- rep(0,nrow(microbe))
  sig.variables[ind] <- 1
  
  if(plots==TRUE){
    postscript(paste(file,"_lasso1.eps",sep=""))
    plot(wLasso.out,cex.axis=1.5,cex.lab=1.5)
    dev.off()
    
    postscript(paste(file,"_lasso2.eps",sep=""))
    par(mar=c(5, 4, 4, 2)+1)
    plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),list(that=delta)), xlab="Steps")
    
    abline(v=p.min,lty=2)
    dev.off()
  }
  list(sig.variables=sig.variables,
       entry.variables=entry.variables,delta.out=delta.out)
}


lasso.computations <- function(weights,microbes,phenotypes,g,plots=TRUE,file="name",include.diet="TRUE",
					diet.wt=100,thresh.q=FALSE,corr.g=FALSE,delta=2,
					cv.criterion=FALSE,vfold=10){

	interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
	interest <- as.data.frame(interest)
	rownames(interest) <- rownames(microbes)
	colnames(interest) <- rownames(phenotypes)

	entry.var <- array(0,dim=c(nrow(phenotypes),nrow(microbes)))

	delta.out <- matrix(0,nrow=1,ncol=nrow(phenotypes))

	for(i in 1:ncol(interest)){
		lasso.out <- lasso(weights[,i],phenotypes[i,],microbes,g,file=paste(file,rownames(phenotypes)[i],sep=""),
					plots=plots,include.diet=include.diet,
					diet.wt=diet.wt,thresh.q=thresh.q,corr.g=corr.g,delta=delta,
					cv.criterion=cv.criterion,vfold=vfold)
		interest[,i] <- lasso.out$sig.variables
			
		entry.var[i,] <- lasso.out$entry.variables

		delta.out[,i] <- lasso.out$delta.out
	
	}
	list(interest=interest,entry.var=entry.var,delta.out=delta.out)
}


###############################################
## Cross-Validation function to select delta ##
###############################################

cv.delta <- function(y1,X1,K=10){
	# sequence of delta values
	delta.cv <- seq(0.75,2,by=0.1)

	# Randomly partition the data
	all.folds <- cv.folds(length(y1), K)

	# Matrix to store residuals
	residmat <- matrix(0, length(delta.cv), K)

	for(j in seq(K)){
		# data set to omit
		omit <- all.folds[[j]]

		# Run Lasso with after removing omit data
		wLasso.out <- lasso.procedure(y1[-omit],X1[-omit,,drop=FALSE])$wLasso.out

		for(d in 1:length(delta.cv)){

			# Find best-fitting model for specified delta
			beta.omit <- lasso.delta.choice(wLasso.out,y1[-omit],X1[-omit,,drop=FALSE],delta=delta.cv[d])

			# Find final fit with data omitted
			fit <- X1[omit,,drop=FALSE] %*% beta.omit$predict.out$coefficients

			# Store residual
			residmat[d,j] <-  apply((y1[omit] - fit)^2, 2, sum)
		}		
	}

	cv <- apply(residmat,1,mean)

	# Check which delta's lead to min(cv)
	delta.ind <- which(cv==min(cv))
	delta.opt <- delta.cv[delta.ind]
      delta <- mean(delta.opt)		## takes average of delta values

      wLasso.out <- lasso.procedure(y1,X1)$wLasso.out
	predict.out <- lasso.delta.choice(wLasso.out,y1,X1,delta=delta)$predict.out	
	list(predict.out=predict.out,delta=delta)
}


lasso.procedure <- function(y1,X1){
	# Setup
	N <- length(y1)

	## adjust use.Gram
	if(ncol(X1)>500){
		use.Gram <- FALSE
	} else {
		use.Gram <- TRUE
	}


	# Run Lasso 
	wLasso.out <-  lars(X1, y1, type = c("lasso"), 
                trace = FALSE, normalize = FALSE, intercept = FALSE,use.Gram=use.Gram)
	
	list(wLasso.out=wLasso.out)
}


lasso.delta.choice <- function(wLasso.out,y1,X1,delta){
	# Setup
	N <- length(y1)

	p = dim(X1)[2]

	s = length(wLasso.out$df)
	p.pos = NULL

	RSS = NULL
	for (i in 1:s){
 			RSS[i] = sum((y1-predict(wLasso.out, X1, s=i, type = c("fit"))$fit)**2)
 			p.pre = predict(wLasso.out, X1, s=i, type = c("coefficients"))$coefficients
 			p.pos = c(p.pos,length(p.pre[abs(p.pre)>0]))
	}

	MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
	MSE <- MSE^2

	p.min = which.min(RSS/MSE+delta*p.pos)


	## final best descriptive model
	predict.out <- predict(wLasso.out, X1, s=p.min, type = c("coefficients"))
	
	list(wLasso.out=wLasso.out,predict.out=predict.out)
}

## Jan. 11, 2011: This program requires no restrictions on alpha's and beta's
##			This program also implements estimating k.
## Implementation of negative-binomial integer-valued
## GARCH model by Fukang Zhu


#############
# Libraries #
#############

library(stats)
library(MASS)	# used for simulation, mvrnorm

##########
# Model: #
##########

# X_t | F_{t-1} : NB(r,p_t)
# \lambda_t = \alpha_0 + \sum_{i=1}^p \alpha_i * X_{t-i} + \sum_{j=1}^q \beta_j * lambda_{t-j}
# k>0, alpha_0 >0, \alpha_i \geq 0, i=1,..,p, \beta_j \geq 0
# theta= (k,alpha_0,...,alpha_p,beta_1,...,beta_p)

####################
# Useful functions #
####################
get.Xlambda <- function(X,covariates,p,q){
	X <- as.matrix(X)
	covariates <- as.matrix(covariates)
	n <- nrow(X)
	Xlambda<- matrix(0,nrow=n,ncol=p+1)				# Xlambda is needed to construct lambda_t
	Xlambda[,1] <-1
	for (i in 1:p){
		zero 		    <- rep(0,p-i+1)
		Xlambda[,p-i+2] <- c(zero,X[1:(n-(p-i+1))]) 
	}
	
	Xlambda <- cbind(Xlambda,covariates)

	if(q==0){
		Xlambda <- Xlambda[(p+1):nrow(Xlambda),]			# Xlambda starts from p+1
	}
	X.use   <- X[(p+1):n,]							# X time series from p+1,...
	list(Xlambda=Xlambda,X.use=as.matrix(X.use))
}


# this is WRONG!!!
#get.lambda <- function(Xlambda,p,q,alpha,beta,se=0){		# lambda when q>0 
#	n 	    <- nrow(Xlambda)
#	lambda    <- rep(0,n)					
#	if((p+1)<=q){
#		for (t in (p+1):q){
#			zeroq     <- rep(0,q-t+1)
#			lambda[t] <- sum(alpha * Xlambda[t,]) + sum(beta *c(lambda[(t-1):1],zeroq) ) 
#		}
#		for (t in (q+1):n){
#			lambda[t] <- sum(alpha * Xlambda[t,]) + sum(beta *lambda[(t-1):(t-q)]) 
#		}
#	}else{
#		for (t in (p+1):n){
#			lambda[t] <- sum(alpha * Xlambda[t,]) + sum(beta *lambda[(t-1):(t-q)]) 
#		}
#	}
#	if(se==1){ 
#		lambda.use <- lambda
#	}else{
#		lambda.use <- lambda[(p+1):n]
#	}
#	return(lambda.use)
#}


# made changes!
get.lambda <- function(Xlambda,p,q,alpha,beta,se=0){		# lambda when q>0 
	n 	    <- nrow(Xlambda)
	lambda.calc    <- rep(0,n+q)					

	lambda.calc[1+q] <- sum(alpha * Xlambda[1,])

	for(t in 2:n){
		lambda.calc[t+q] <- sum(alpha *Xlambda[t,]) + sum ( beta * lambda.calc[(t+q-1):(t+q-q)] )
	}

	lambda <- lambda.calc[(q+1):length(lambda.calc)]
	
	if(se==1){ 
		lambda.use <- lambda
	}else{
		lambda.use <- lambda[(p+1):length(lambda)]
	}
	return(lambda.use)
}




##################
# Initial Values #
##################

# input:
# X 	: time-series data
# p,q : parameters from NBINGARCH(p,q)
# r	: parameter from NB model

make.f.init <- function(X,covariates,p,q){
	out      <- get.Xlambda(X,covariates,p,q)
	Xlambda <- out$Xlambda
	X.use    <- out$X.use	
	covariates <- as.matrix(covariates)	
	m.cov    <- ncol(covariates)	

	function(theta){
		if(q==0){
			r     <- theta[1]
			alpha <- theta[2:length(theta)]
			lambda <- Xlambda %*% alpha				# Construct lambda_t
		} else{
			r     <- theta[1]
			alpha <- theta[2:(p+2+m.cov)]
			beta  <- theta[(p+3+m.cov):length(theta)]	
			lambda    <- get.lambda(Xlambda,p,q,alpha,beta)	# Construct lambda_t
		}

		max.lambda <- rep(100,length(lambda))
		control.lambda <- apply(cbind(max.lambda,lambda),1,min)		# to control input to gamma(), not make it too large!

				
		out <- sum ( (X.use - exp(control.lambda)) ^ 2 ) 
		return(out)
	}
}


########################
# Function to Maximize #
########################

# input:
# X 	: time-series data
# p,q : parameters from NBINGARCH(p,q)

make.f <- function(X,covariates,p,q){
	out <- get.Xlambda(X,covariates,p,q)
	Xlambda <- out$Xlambda
	X.use   <- out$X.use
	n      <- length(X)
	covariates <- as.matrix(covariates)
	m.cov   <- ncol(covariates)

	function(theta){								# function that will be optimized if q=0
		if(q==0){
			r     <- theta[1]
			alpha <- theta[2:length(theta)]
			lambda <- Xlambda %*% alpha				# Construct lambda_t
		} else{
			r     <- theta[1]
			alpha <- theta[2:(p+2+m.cov)]
			beta  <- theta[(p+3+m.cov):length(theta)]	
			lambda    <- get.lambda(Xlambda,p,q,alpha,beta)	# Construct lambda_t
		}
		max.gamma <- rep(100,length(X.use))
		control.gamma <- apply(cbind(max.gamma,X.use+exp(r)),1,min)		# to control input to gamma(), not make it too large!
		
		out <- sum( -log(factorial(X.use))+ log(gamma(control.gamma)) + X.use * lambda - (exp(r)+X.use) * log(exp(r)+exp(lambda)) )- (n-p)*(log(gamma(exp(r))) - r* exp(r))

		out <- -out
		return(out)
	}
}

make.f.simu <- function(X,covariates,p,q){
	out <- get.Xlambda(X,covariates,p,q)
	Xlambda <- out$Xlambda
	X.use   <- out$X.use
	n      <- length(X)
	covariates <- as.matrix(covariates)
	m.cov   <- ncol(covariates)

	function(theta){								# function that will be optimized if q=0
		if(q==0){
			r     <- theta[1]
			alpha <- theta[2:length(theta)]
			lambda <- Xlambda %*% alpha				# Construct lambda_t
		} else{
			r     <- theta[1]
			alpha <- theta[2:(p+2+m.cov)]
			beta  <- theta[(p+3+m.cov):length(theta)]	
			lambda    <- get.lambda(Xlambda,p,q,alpha,beta)	# Construct lambda_t
		}
		max.gamma <- rep(100,length(X.use))
		control.gamma <- apply(cbind(max.gamma,X.use+exp(r)),1,min)		# to control input to gamma(), not make it too large!
		
		#out <- sum( -log(factorial(X.use))+ log(gamma(control.gamma)) + X.use * lambda - (exp(r)+X.use) * log(exp(r)+exp(lambda)) )- (n-p)*(log(gamma(exp(r))) - r* exp(r))
		out <- sum( log(gamma(control.gamma)) + X.use * lambda - (exp(r)+X.use) * log(exp(r)+exp(lambda)) )- (n-p)*(log(gamma(exp(r))) - r* exp(r))

		out <- -out
		return(out)
	}
}



##############################
# Asymptotic Standard Errors #
##############################

# Input:
# theta 	: parameter estimates
# X		: time series data
# p,q 	: parameters from NBINGARCH(p,q)

se.mle <- function(theta,X,covariates,p,q){
	X <-as.matrix(X)
	n <- nrow(X)					# length of time series

	r <- theta[1]
	covariates <- as.matrix(covariates)
	m.cov <- ncol(covariates)
	
	if(q==0){
		Xlambda 	<- matrix(0,nrow=n,ncol=p+1)		# Xlambda is needed to construct lambda_t, not cut off from 1 to p
		Xlambda[,1] <-1
		for (i in 1:p){
			zero 		    <- rep(0,p-i+1)
			Xlambda[,p-i+2] <- c(zero,X[1:(n-(p-i+1))]) 
		}
		Xlambda <- cbind(Xlambda,covariates)

		alpha <- theta[2:length(theta)]			# Construct lambda_t
		lambda <- Xlambda %*% alpha	

		###########################
		# partial ell / partial r #
		###########################

		temp.r <- exp(r)* ( digamma(X+exp(r))-digamma(exp(r)) + 1+r - log(exp(lambda) + exp(r))-(exp(r)+X)/ ( exp(lambda)+exp(r) ) )  

		###############################
		# partial ell^2 / partial r^2 #
		###############################
		
		temp.rr <- exp(r) * ( digamma(X + exp(r) )-digamma(exp(r)) )+
			     exp(r)^2 * ( trigamma(X + exp(r)) - trigamma(exp(r)))+
			     exp(r) * (2+r) - 
			     exp(r) * log ( exp(lambda) + exp(r) ) -
			     exp(r) * (2*exp(r)  + X ) /(exp(lambda) + exp(r) ) - 
	       	     exp(r)^2 * ( exp(lambda)-X ) / ( exp(lambda) + exp(r) )^2


		S.hat <- diag(0,nrow=length(theta))			# We build \hat S, and \hat D
		D.hat <- diag(0,nrow=length(theta))

		for(i in (p+1):n){
			temp1 <- ( X[i] - ( exp(r)+ X[i] ) * exp(lambda[i]) / (exp(r) + exp(lambda[i]) ) )*Xlambda[i,]
			temp2 <- c(temp.r[i],temp1)
			S.hat <- S.hat + temp2 %*% t(temp2)

			temp3 <- Xlambda[i,] %*% t(Xlambda[i,])	
			temp.theta.r <- -exp(lambda[i]) * exp(r) * (exp(lambda[i])-X[i])/ (exp(lambda[i])+exp(r))^2
			temp.theta.2 <- (exp(r)+X[i]) * exp(lambda[i])*exp(r)/ ( exp(lambda[i])+exp(r) )^2

			A <- diag(0,nrow=length(theta))

			A[1,] <- c(temp.rr[i], temp.theta.r *Xlambda[i,])
			A[,1] <- t(A[1,])
			A[2:length(theta),2:length(theta)] <- -temp.theta.2 * temp3
			D.hat <- D.hat + A
			}
	} else{
			alpha <- theta[2:(p+2+m.cov)]
			beta  <- theta[(p+3+m.cov):length(theta)]	
			
			######################
			# Construct lambda_t #
			######################
			Xlambda <- get.Xlambda(X,covariates,p,q)$Xlambda	
			lambda    <- get.lambda(Xlambda,p,q,alpha,beta,se=1)	# Construct lambda_t

			###########################
			# partial ell / partial r #
			###########################

			temp.r <- exp(r)* ( digamma(X+exp(r))-digamma(exp(r)) + 1+r - log(exp(lambda) + exp(r))-(exp(r)+X)/ ( exp(lambda)+exp(r) ) )  

			###############################
			# partial ell^2 / partial r^2 #
			###############################
		
			temp.rr <- exp(r) * ( digamma(X + exp(r) )-digamma(exp(r)) )+
			     exp(r)^2 * ( trigamma(X + exp(r)) - trigamma(exp(r)))+
			     exp(r) * (2+r) - 
			     exp(r) * log ( exp(lambda) + exp(r) ) -
			     exp(r) * (2*exp(r)  + X ) /(exp(lambda) + exp(r) ) - 
	       	     exp(r)^2 * ( exp(lambda)-X ) / ( exp(lambda) + exp(r) )^2


			######################
			# Construct L.lambda #
			######################
			L.lambda 	<- matrix(0,nrow=n,ncol=q)		# L.lambda is needed to construct partial derivatives
			for (i in 1:q){
				zero 		    <- rep(0,q-i+1)
				L.lambda[,q-i+1] <- c(zero,lambda[1:(n-(q-i+1))]) 
			}

			##############################################
			# Construct partial lambda_t / partial alpha #
			##############################################
			#par.alpha        <- matrix(0,nrow=n,ncol=p+m.cov+1)			# constructing partial lambda_t / partial alpha
			#par.alpha[1,]    <- Xlambda[1,] 
			#for (j in 1:(1+p+m.cov)){
			#	if((p+1)<=q){
			#		for(t in (p+1):q){
			#			zeroq          <- rep(0,q-t+1)
			#			par.alpha[t,j] <- Xlambda[t,j] + sum (beta * c(par.alpha[(t-1):1,j],zeroq) ) 
			#		}
			#		for (t in (q+1):n){
			#			par.alpha[t,j] <- Xlambda[t,j] + sum( beta * par.alpha[(t-1):(t-q),j] )
			#		}
			#	}else{
			#		for (t in (p+1):n){
			#			par.alpha[t,j] <- Xlambda[t,j] + sum( beta * par.alpha[(t-1):(t-q),j] )
			#		}
			#	}
			#}

			par.alpha.tmp <- matrix(0,nrow=n+q,ncol=p+m.cov+1)

			par.alpha.tmp[1+q,1] <- 1
			for (t in 2:n){
				par.alpha.tmp[t+q,1] <- 1+ sum( beta * par.alpha.tmp[(t+q-1):t,1])
			}

			for( j in 2:(p+m.cov+1)){
				for( t in (j+1):(n+1)){
					par.alpha.tmp[t+q-1,j] <- X[t-j] + sum(beta * par.alpha.tmp[(t+q-2):(t-1),j] )
				}
			}
			
			par.alpha <- par.alpha.tmp[(1+q):(n+q),]



			#############################################
			# Construct partial lambda_t / partial beta #
			#############################################
			#par.beta	     <- matrix(0,nrow=n,ncol=q)			# constructing partial lambda_t / partial beta	
			#for (j in 1:q){
			#	if((p+1)<=q){
			#		for(t in (p+1):q){
			#			zeroq          <- rep(0,q-t+1)
			#			par.beta[t,j]  <- L.lambda[t,j] + sum(beta * c(par.beta[(t-1):1,j],zeroq) )		
			#		}
			#		for (t in (q+1):n){
			#			par.beta[t,j]  <- L.lambda[t,j] + sum( beta * par.beta[(t-1):(t-q),j] )
			#		}
			#	}else{
			#		for (t in (p+1):n){
			#			par.beta[t,j]  <- L.lambda[t,j] + sum( beta * par.beta[(t-1):(t-q),j] )
			#		}
			#	}
			#}

			par.beta.tmp <- matrix(0,nrow=n+q,ncol=q)
			for(j in 1:q){
				for(t in (j+1):n){
					par.beta.tmp[t+q,j] <- lambda[t-j] + sum(beta*par.beta.tmp[(t+q-1):t,j] )
				}
			}

			par.beta <- as.matrix(par.beta.tmp[(1+q):(n+q),1:q])

			

			###################
			# Construct S.hat #
			###################
			S.hat <- diag(0,nrow=length(theta))			# We build \hat S

			for(i in (p+1):n){
				temp0 <- ( X[i] - ( exp(r)+ X[i] ) * exp(lambda[i]) / (exp(r) + exp(lambda[i]) ) )*
						c(par.alpha[i,],par.beta[i,])
				temp1 <- c(temp.r[i],temp0)
				S.hat <- S.hat+ temp1 %*% t(temp1) 
			}

			###################
			# Construct D.hat #
			###################
			D.hat <- diag(0,nrow=length(theta))

			#par2.ab <- array(0,dim=c(n,p+m.cov+1,q))
			#par2.bb <- array(0,dim=c(n,q,q))
		
			#for (i in 1:(p+m.cov+1)){
			#	for (j in 1:q){
			#		for (t in max(j+1,p+1):n){
			#			par2.ab[t,i,j] <- par.alpha[t-j,i] + beta[j] * par2.ab[t-j,i,j]
			#		}
			#	}
			#}
		
			#for (j in 1:q){
			#	for(l in 1:q){
			#		for(t in max(l+1,p+m.cov+1):n){
			#			par2.bb[t,j,l] <- par.beta[t-l,j] + beta[l] * par2.bb[t-l,j,l]
			#		}
			#	}
			#}

			#for (j in 1:q){
			#	for (l in 1:q){
			#		for(t in max(j+1,p+m.cov+1):n){
			#			par2.bb[t,j,l] <- par2.bb[t,j,l]  + par.beta[t-j,l]
			#		}
			#	}
			#}			
			
			par2.ab.tmp <- array(0,dim=c(n+q,p+m.cov+1,q))
			for(j in 1:(p+m.cov+1)){
				for(i in 1:q){
					for(t in (i+1):n){
						par2.ab.tmp[t+q,j,i] <- par.alpha[t-i,j] + sum(beta*par2.ab.tmp[(t+q-1):t,j,i])
					}
				}
			}
			par2.ab <- par2.ab.tmp[(1+q):(n+q),,]	

			par2.bb <- array(0,dim=c(n,q,q))
			
			for(l in 1:q){
				for(j in 1:q){
					for(t in (j+1):n){
						par2.bb[t,j,l] <- par.beta[t-j,l] + par2.bb[t,j,l]
					}
				}
			}


			for(j in 1:q){
				for(l in 1:q){
					for(t in (l+1):n){
						par2.bb[t,j,l] <- par.beta[t-l,j] + par2.bb[t,j,l]
					}
				}
			}

			for(l in 1:q){
				for(j in 1:q){
					for(s in 1:q){
						for(t in (s+1):n){
							par2.bb[t,j,l] <- par2.bb[t,j,l] + beta[s]*par2.bb[t-s,j,l]
						}
					}
				}
			}

			

			for(s in (p+1):n){
				par.theta <- c(par.alpha[s,],par.beta[s,])
				temp.theta.1 <- X[s] - (exp(r) + X[s]) * exp(lambda[s])/ (exp(lambda[s])+exp(r) )
				temp.theta.2 <- (exp(r)+X[s]) * exp(lambda[s])*exp(r)/ ( exp(lambda[s])+exp(r) )^2
				

				temp.theta.r <- -exp(lambda[s]) * exp(r) * (exp(lambda[s])-X[s])/ (exp(lambda[s])+exp(r))^2

				
				A <- diag(0,nrow=length(theta))
				A[1,] <- c(temp.rr[s],temp.theta.r*par.theta)
				A[,1] <- t(A[1,])
				A[2:(2+p+m.cov),2:(2+p+m.cov)] <- - temp.theta.2 * par.alpha[s,] %*% t(par.alpha[s,])
				if(q==1){
					A[2:(2+p+m.cov),(p+m.cov+3):(p+m.cov+q+2)] <- temp.theta.1 * par2.ab[s,] - temp.theta.2 * par.alpha[s,] %*% t(par.beta[s,])
				}else{
					A[2:(2+p+m.cov),(p+m.cov+3):(p+m.cov+q+2)] <- temp.theta.1 * par2.ab[s,,] - temp.theta.2 * par.alpha[s,] %*% t(par.beta[s,])
				}
				A[(p+m.cov+3):(p+m.cov+q+2),2:(2+p+m.cov)] <- t(A[2:(2+p+m.cov),(p+m.cov+3):(p+m.cov+q+2)])
				A[(p+m.cov+3):(p+m.cov+q+2),(p+m.cov+3):(p+m.cov+q+2)] <- temp.theta.1 * par2.bb[s,,] - temp.theta.2 * par.beta[s,] %*% t(par.beta[s,])

				D.hat <- D.hat + A
			}
	}
	S.hat <- S.hat/(n-p)
	D.hat <- -D.hat/(n-p)
	cov <- solve ( D.hat %*% solve(S.hat) %*% t(D.hat) ) 		
	cov <- cov/(n-p)
	return(cov)	
}


#####################
# Pearson Residuals #
#####################

pears.resid <- function(X,covariates,theta,p,q){
	out     <- get.Xlambda(X,covariates,p,q)
	Xlambda <- out$Xlambda
	X.use   <- out$X.use
	covariates <- as.matrix(covariates)
	m.cov   <- ncol(covariates)	

	if(q==0){
		r     <- theta[1]
		alpha <- theta[2:length(theta)]
		lambda <- Xlambda %*% alpha					# Construct lambda_t
	} else{
		r     <- theta[1]
		alpha <- theta[2:(p+2+m.cov)]
		beta  <- theta[(p+m.cov+3):length(theta)]	
		lambda<- get.lambda(Xlambda,p,q,alpha,beta)		# Construct lambda_t
	}
	resid <- ( X.use - exp(lambda) ) / sqrt( exp(lambda) * (1 + exp(lambda) /exp(r) ) )  
}


###############################
# Un-Constrained Optimization #
###############################
# This function will calculate MLE and plot residuals

# iprint : if TRUE, will plot residuals, ACF, and PACF in one plot; if FALSE, will produce 3 separate pdf files
# file   : number for file name

optim.resid <- function(X,covariates,p,q=0,theta0,iprint=TRUE,file=1,digits=3){
	X    <- as.matrix(X)
	covariates <- as.matrix(covariates)
	num  <- p+q+2+ncol(covariates)
	if(length(theta0)!= num) return("not valid theta0")

	# Get initial value via CLS
	f.init     <- make.f.init(X,covariates,p,q)
	init.optim <- optim(theta0,f.init)
	theta.init <- init.optim$par	

	# Get MLE
	f.max     <- make.f.simu(X,covariates,p,q)
	out.optim <- optim(theta.init,f.max)		# start at CLS result
	#out.optim <- optim(theta0,f.max)		# start at initial

	theta.hat <- out.optim$par
	loglik    <- -out.optim$value
	AIC       <- 2*length(theta.hat) - 2 *loglik
	BIC	    <- -2 * loglik +length(theta.hat) * log((length(X)-p))
	c.ini     <- init.optim$convergence
	c.mle     <- out.optim$convergence	

	#est.var <- rep(0,length(theta.hat))
	#ratio <- rep(0,length(theta.hat))

	est.var   <- sqrt(diag(se.mle(theta.hat,X,covariates,p,q)))
	ratio     <- abs(theta.hat) / est.var		
	results   <- rbind(theta.init,theta.hat,est.var,ratio,c(loglik,rep(NA,num-1)),c(AIC,rep(NA,num-1)),c(BIC,rep(NA,num-1)),c(c.ini,rep(NA,num-1)),c(c.mle,rep(NA,num-1)))
	rownames(results) <- c("est.ini","est","se","ratio","loglik","AIC","BIC","c.ini","c.mle")
	colnames(results) <- c("log(k)",c(rep("alpha",p+1),rep("gamma",ncol(covariates)),rep("beta",q)))


	pears <- pears.resid(X,covariates,theta.hat,p,q)
	time  <- seq((p+1):length(X))

	if(iprint==TRUE){
		par(mfrow=c(3,1))
		plot(time,pears,lty=2,main="",cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Pearson Residuals")
		acf(pears,main="",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
		pacf(pears,main="",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
	}else{
		postscript(paste("pearresid",file,".ps",sep=""))
		plot(time,pears,lty=2,main="Pearson Residuals") 
		dev.off()

		postscript(paste("acfresid",file,".ps",sep=""))
		acf(pears,main="ACF Pearson Residuals")
		dev.off()
		
		postscript(paste("pacfresid",file,".ps",sep=""))
		pacf(pears,main="PACF Pearson Residuals")
		dev.off()
	}
	list(results=round(results,digits=digits),resid=pears)
}




##############
# Prediction #
##############

# X		   : truncated time series
# X.full	   : full time series
# out		   : results from optim.resid function
# n.pred	   : number of predicted values
# B		   : number of bootstrap samples
# alpha1	   : for confidence interval

# only works for q=0!!!
predict.nbin <- function(X,X.full,covariates.full,p,q=0,out,n.pred=10,alpha1=0.05,B=10000,future=FALSE){
	X <- as.matrix(X)
	X.new <- X								# copy of time series X
	X.full <- as.matrix(X.full)
	covariates <- as.matrix(covariates.full)
	
	n <- nrow(X)
	k <- exp(out$results["est",1])
	alpha <- out$results["est",2:(p+2)]					# estimated alpha's
	gamma <- out$results["est",(p+3):(p+3+ncol(covariates)-1)]
	if(q>0){
		beta  <- out$results["est", (p+3+ncol(covariates)):ncol(out$results)]
	}

	pred.ts <- matrix(NA,nrow=n.pred,ncol=5)			# matrix to store estimates, variance, confidence interval, MSE
	colnames(pred.ts) <- c("pred","var","l.bound","u.bound","MSE")

	for(i in 1:n.pred){
		t <- n+i
		X.temp   <- c(1,X.new[(t-1):(t-p)]) 		# get values X_{t-1},...,X_{t-p}}
		lambda.t <- sum(alpha * X.temp)+sum(gamma*covariates[t])			# Calculate lambda_t
		samp     <- rnbinom(B, size=k, mu=exp(lambda.t))		# Sample B many negative binomials
		sort.samp <- sort(samp)

		#pred.ts[i,"pred"]     <- mean(sort.samp)		# Taking mean as prediction, may not be integer!
		pred.ts[i,"pred"]     <- median(sort.samp)	# Taking median as prediction, will be integer!
		#pred.ts[i,"pred"]	    <- round(mean(sort.samp))	# Rounding to nearest integer   
		pred.ts[i,"var"]      <- var(sort.samp)
		pred.ts[i,"l.bound"]  <- pred.ts[i,"pred"] - 2 * sqrt(pred.ts[i,"var"])
		pred.ts[i,"u.bound"]  <- pred.ts[i,"pred"] + 2 * sqrt(pred.ts[i,"var"])
		# pred.ts[i,"l.bound"]  <- sort.samp[B*alpha1/2]
		# pred.ts[i,"u.bound"]  <- sort.samp[B*(1-alpha1/2)]
	
		X.new <- c(X.new,pred.ts[i,"pred"])
	
		if(future==FALSE){
			pred.ts[i,"MSE"] <- pred.ts[i,"var"] + (pred.ts[i,"pred"]-X.full[t])^2
		}
	}

	if(future==FALSE){
		y.lim.min <- min(X.new,pred.ts[,"l.bound"])
		y.lim.max <- max(X.new,pred.ts[,"u.bound"])

		n.total <- n+n.pred
		index   <- 1:n.total
		index1  <- (n+1):(n+n.pred)	
	
		plot(index,X.full[index],type="l",main="True vs. Predicted Number of Crashes \n NBINGARCH Model",ylim=c(y.lim.min,y.lim.max),xlab="Day",ylab="")

		lines(index1,pred.ts[,"pred"],col="red")
		lines(index1,pred.ts[,"l.bound"],col="blue",lty="dashed")
		lines(index1,pred.ts[,"u.bound"],col="blue",lty="dashed")

		MSE.final <- sum(pred.ts[,"MSE"])

	}else{
		y.lim.min <- min(X.new,pred.ts[,"l.bound"])
		y.lim.max <- max(X.new,pred.ts[,"u.bound"])

		n.total   <- n+n.pred
		index     <- (n.total-100):n.total
		index1    <- (n+1):(n+n.pred)

		plot(index,X.new[index],type="l",main="True vs. Predicted Number of Crashes\n NBINGARCH Model",ylim=c(y.lim.min,y.lim.max),xlab="Day",ylab="")
	
		lines(index1,pred.ts[,"pred"],col="red")
		lines(index1,pred.ts[,"l.bound"],col="green",lty=2)
		lines(index1,pred.ts[,"u.bound"],col="green",lty=2)

		MSE.final <- 0

	}
	list(pred.ts=pred.ts,MSE=MSE.final,X.new=X.new)
}

















###############
# Simulations #
###############


optim.resid.simu <- function(X,covariates,p,q=0,theta0){
	X    <- as.matrix(X)
	covariates <- as.matrix(covariates)
	num  <- p+q+2+ncol(covariates)
	if(length(theta0)!= num) return("not valid theta0")

	# Get initial value via CLS
	f.init     <- make.f.init(X,covariates,p,q)
	init.optim <- optim(theta0,f.init)
	theta.init <- init.optim$par	

	# Get MLE
	f.max     <- make.f.simu(X,covariates,p,q)
	#out.optim <- optim(theta.init,f.max)
	out.optim <- optim(theta0,f.max)

	theta.hat <- out.optim$par
	loglik    <- -out.optim$value
	AIC       <- 2*length(theta.hat) - 2 *loglik
	BIC	    <- -2 * loglik +length(theta.hat) * log((length(X)-p))
	c.ini     <- init.optim$convergence
	c.mle     <- out.optim$convergence	

	#est.var <- rep(0,length(theta.hat))
	#ratio <- rep(0,length(theta.hat))

	est.var   <- sqrt(diag(se.mle(theta.hat,X,covariates,p,q)))
	#ratio     <- abs(theta.hat) / est.var		
	results   <- c(theta.init,theta.hat,est.var,c.ini,c.mle)
	#rownames(results) <- c("est.ini","est","se","ratio","loglik","AIC","BIC","c.ini","c.mle")
	
	
	return(results)
}

# get covariates
get.covar <- function(gamma,n){
	m <- length(gamma)
	mu <- rep(0,m)
	sigma <- diag(1,nrow=m)
	cov.val <- mvrnorm(n,mu,sigma)
	return(cov.val)
}

get.coval<- function(gamma,covar){
	return(sum(gamma*covar))
}


# gen.data.simple NOT flexible!!
# ge.data.simple ONLY handles p=1 and p=2, q=1,q=2!!!
# gen.data.simpley ONLY handles length(gamma)=2!!!
gen.data.simple <- function(theta,p,q=0,n,k,covar){
	if(q==0){
		if(p==1){
			a0 <- theta[1]
			a1 <- theta[2]
			gamma <- theta[3:length(theta)]
			x <- rep(0,n)			
			x[1]=rnbinom(1,size=k,mu=exp(a0+get.coval(gamma,covar[1,])))
			for (i in 2:n){
				lambda=a0+a1*x[i-1]+get.coval(gamma,covar[i,])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda))
			}
		}else{
			a0 <- theta[1]
			a1 <- theta[2]
			a2 <- theta[3]
			gamma <- theta[4:length(theta)]
			x <- rep(0,n)
			x[1]<- rnbinom(1,size=k,mu=exp(a0+get.coval(gamma,covar[1,])))
			x[2]<- rnbinom(1,size=k,mu=exp(a0+a1*x[1]+get.coval(gamma,covar[2,])))
			for (i in 3:n){
				lambda=a0+a1*x[i-1]+a2*x[i-2]+get.coval(gamma,covar[i,])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda))
			}
		}
	}
	if(q==1){
		if(p==1){
			a0 <- theta[1]
			a1 <- theta[2]
			gamma <- theta[3:(length(theta)-1)]
			b1 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0+get.coval(gamma,covar[1,])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			
		
			for (i in 2:n){
				lambda[i] <- a0+a1*x[i-1] + b1*lambda[i-1]+get.coval(gamma,covar[i,])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}else{
			a0 <- theta[1]
			a1 <- theta[2]
			a2 <- theta[3]
			gamma <- theta[4:(length(theta)-1)]
			b1 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0 + get.coval(gamma,covar[1,])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			lambda[2] <- a0+a1*x[1] +b1*lambda[1] + get.coval(gamma,covar[2,])
			x[2]<- rnbinom(1,size=k,mu=exp(lambda[2]))
			
			for (i in 3:n){
				lambda[i]=a0+a1*x[i-1]+a2*x[i-2] + b1*lambda[i-1]+get.coval(gamma,covar[i,])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}
	}
	if(q==2){
		if(p==1){
			a0 <- theta[1]
			a1 <- theta[2]
			gamma <- theta[3:(length(theta)-2)]
			b1 <- theta[(length(theta)-1)]
			b2 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0 + get.coval(gamma,covar[1,])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			
		
			for (i in 2:n){
				if(i==2){
					lambda[i] <- a0+a1*x[i-1] + b1*lambda[i-1] + get.coval(gamma,covar[i,])
				}else{
					lambda[i] <- a0+a1*x[i-1] + b1*lambda[i-1]+b2*lambda[i-2]+ get.coval(gamma,covar[i,])
				}
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}else{
			a0 <- theta[1]
			a1 <- theta[2]
			a2 <- theta[3]
			gamma <- theta[4:(length(theta)-2)]
			b1 <- theta[(length(theta)-1)]
			b2 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0+get.coval(gamma,covar[1,])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			lambda[2] <- a0+a1*x[1] +b1*lambda[1]+get.coval(gamma,covar[2,])
			x[2]<- rnbinom(1,size=k,mu=exp(lambda[2]))
			
			for (i in 3:n){
				lambda[i]=a0+a1*x[i-1]+a2*x[i-2] + b1*lambda[i-1]+b2*lambda[i-2]+
						get.coval(gamma,covar[i,])	
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}
	}
	return(x)
}


gen.data.simple2 <- function(theta,p,q=0,n,k,covar){
	if(q==0){
		if(p==1){
			a0 <- theta[1]
			a1 <- theta[2]
			gamma <- theta[3:length(theta)]
			x <- rep(0,n)			
			x[1]=rnbinom(1,size=k,mu=exp(a0+get.coval(gamma,covar[1])))
			for (i in 2:n){
				lambda=a0+a1*x[i-1]+get.coval(gamma,covar[i])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda))
			}
		}else{
			a0 <- theta[1]
			a1 <- theta[2]
			a2 <- theta[3]
			gamma <- theta[4:length(theta)]
			x <- rep(0,n)
			x[1]<- rnbinom(1,size=k,mu=exp(a0+get.coval(gamma,covar[1])))
			x[2]<- rnbinom(1,size=k,mu=exp(a0+a1*x[1]+get.coval(gamma,covar[2])))
			for (i in 3:n){
				lambda=a0+a1*x[i-1]+a2*x[i-2]+get.coval(gamma,covar[i])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda))
			}
		}
	}
	if(q==1){
		if(p==1){
			a0 <- theta[1]
			a1 <- theta[2]
			gamma <- theta[3:(length(theta)-1)]
			b1 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0+get.coval(gamma,covar[1])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			
		
			for (i in 2:n){
				lambda[i] <- a0+a1*x[i-1] + b1*lambda[i-1]+get.coval(gamma,covar[i])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}else{
			a0 <- theta[1]
			a1 <- theta[2]
			a2 <- theta[3]
			gamma <- theta[4:(length(theta)-1)]
			b1 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0 + get.coval(gamma,covar[1])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			lambda[2] <- a0+a1*x[1] +b1*lambda[1] + get.coval(gamma,covar[2])
			x[2]<- rnbinom(1,size=k,mu=exp(lambda[2]))
			
			for (i in 3:n){
				lambda[i]=a0+a1*x[i-1]+a2*x[i-2] + b1*lambda[i-1]+get.coval(gamma,covar[i])
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}
	}
	if(q==2){
		if(p==1){
			a0 <- theta[1]
			a1 <- theta[2]
			gamma <- theta[3:(length(theta)-2)]
			b1 <- theta[(length(theta)-1)]
			b2 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0 + get.coval(gamma,covar[1])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			
		
			for (i in 2:n){
				if(i==2){
					lambda[i] <- a0+a1*x[i-1] + b1*lambda[i-1] + get.coval(gamma,covar[i])
				}else{
					lambda[i] <- a0+a1*x[i-1] + b1*lambda[i-1]+b2*lambda[i-2]+ get.coval(gamma,covar[i])
				}
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}else{
			a0 <- theta[1]
			a1 <- theta[2]
			a2 <- theta[3]
			gamma <- theta[4:(length(theta)-2)]
			b1 <- theta[(length(theta)-1)]
			b2 <- theta[length(theta)]
			x <- rep(0,n)
			lambda <- rep(0,n)
			lambda[1] <- a0+get.coval(gamma,covar[1])
			x[1]<- rnbinom(1,size=k,mu=exp(lambda[1]))
			lambda[2] <- a0+a1*x[1] +b1*lambda[1]+get.coval(gamma,covar[2])
			x[2]<- rnbinom(1,size=k,mu=exp(lambda[2]))
			
			for (i in 3:n){
				lambda[i]=a0+a1*x[i-1]+a2*x[i-2] + b1*lambda[i-1]+b2*lambda[i-2]+
						get.coval(gamma,covar[i])	
				#print(lambda)
				#x[i]=rnbinom(1,mu=lambda,size=k)
				x[i]=rnbinom(1,size=k,mu=exp(lambda[i]))
			}
		}
	}
	return(x)
}





simulation.study <- function(theta,p,q=0,n,simu,covar){
	k <- exp(theta[1])
	alphabeta <- theta[2:length(theta)]
	output <- matrix(0,ncol=(3*length(theta)+2),nrow=simu)
	set.seed(1)
	for(i in 1:simu){
		#data <- gen.data.simple(alphabeta,p,q,n,k,covar)
		data <- gen.data.simple2(alphabeta,p,q,n,k,covar)
		output[i,] <-optim.resid.simu(data,covar,p,q,theta)
		print(i)
	}
	return(output)
}

# this function only outputs good convergence results
simulation.study2 <- function(theta,p,q=0,n,simu,covar){
	k <- exp(theta[1])
	alphabeta <- theta[2:length(theta)]
	output <- matrix(0,ncol=(3*length(theta)+2),nrow=simu)
	set.seed(1)
	i <- 1
	while(i <= simu){
		#data <- gen.data.simple(alphabeta,p,q,n,k,covar)
		data <- gen.data.simple2(alphabeta,p,q,n,k,covar)
		output[i,] <-optim.resid.simu(data,covar,p,q,theta)
		if(output[i,ncol(output)]<1){
			i <- i+1
		}
		#print(i)
	}
	return(output)
}




simulation.output <- function(theta,p,q=0,out){
	simu.summary <- matrix(0,ncol=(length(theta)*2),nrow=4)
	colnames(simu.summary) <- c("log(k)",c(rep("alpha",p+1),rep("gamma",1),rep("beta",q)),"log(k)",c(rep("alpha",p+1),rep("gamma",1),rep("beta",q)))
	rownames(simu.summary) <- c("mean","bias","emp.se","est.se")
	simu.summary["mean",] <- apply(out[,1:(2*length(theta))],2,mean)
	simu.summary["bias",] <- apply(out[,1:(2*length(theta))],2,mean)-c(theta,theta)
	simu.summary["emp.se",] <- apply(out[,1:(2*length(theta))],2,sd)
	simu.summary["est.se",] <- c(rep(0,length(theta)),apply(out[,(2*length(theta)+1):(ncol(out)-2)],2,mean))
	#return(simu.summary)

	conv.summary <- apply(out[,(ncol(out)-1):ncol(out)],2,sum)	# convergence results; first is for initial estimate and second is for NBINGARCH MLE approach; both should be 0!

	list(simu.summary=simu.summary,conv.summary=conv.summary)

}	










###############
# Source code #
###############

source("qw_simulation_function.R")

########################
## Parameter settings ##
########################

m <- 40				# or m <- 1000 (In R code, we call `m1' as 'm'.	
coef<- c(4.5,3,-3,-3,3)		# or coef <- c(2.5,1.5,-1.5,-1.5,1.5)
x.star.correlated <- FALSE 	# or x.star.correlated <- TRUE (will include x_1^* so we have Groups 1-5)
lasso.mincp <- TRUE		# or lasso.mincp <- FALSE
lasso.delta.cv.mult <- FALSE 	# or lasso.delta.cv.mult <- TRUE to get repeated 10-fold cross-validation.
pi0.true <- FALSE


nsimu <- 1000	## number of simulations

## value of pi0.val
if(x.star.correlated==TRUE){
  pi0.val <-  (((m/4*3+1)-5+1) + ( m - (m/4*3+2)+1 ) + 1)/(m +1)
} else {
  pi0.val <-  (((m/4*3+1)-5+1) + ( m - (m/4*3+2)+1 ) )/m
}

## method for estimating pi, either smoother or bootstrap
method <- "smoother"

## thresholding alpha
alpha <- 0.15

## M_n(delta,p) criterion
delta <- 2

## Number of folds for CV criterion
vfold <- 10

## thresholding alpha range
alpha.range <- c(.01,.025,.05,.075,.10,.15,.20)

## range of delta values in M_n(delta,p) criterion
delta.range <- c(.25,.4,.5,.75,1,1.5,2)

## range of alpha values for Benjamini-Hochberg
alpha.bh.range <- c(0.01,.025,0.05,.075,0.1,0.11,0.12,0.13,0.14,0.15,0.2)

## t-test used to compute p-values
ttest <- TRUE

## number of mice in each diet group
n.i <- 20

## number of diet groups
g <- 2

## seed for random simulations
seed <- 1110

## weight on diet variable
diet.wt <- 1000

## we cut off small weight values so we avoid division by 0
thresh.q <- TRUE

## correlation between x1 and x1^*
rho <- 0.8

##
if(lasso.delta.cv.mult==TRUE){
  ## for testing!
  nsimu <- 1
}

## number of repeated cross-validations
ncv <- 1

## range of percents for multiple cross-validated data
percents.range <-c(50,60,70,80,90,100)

## number of significant digits
sig.digits <- 3


####################
# Simulation study #
####################
## Function to call simulation study
out2.FDR <-  simu.study(rho = rho, coef = coef,nsimu = nsimu,alpha = alpha,
                        delta = delta, vfold = vfold,
                        alpha.range = alpha.range, delta.range = delta.range,
                        alpha.bh.range = alpha.bh.range,
                        ttest = ttest,
                        method = method,m = m ,n.i = n.i,g = g,seed = seed,diet.wt = diet.wt,thresh.q=thresh.q,
                        pi0.val = pi0.val, pi0.true = pi0.true, 
                        x.star.correlated = x.star.correlated, lasso.mincp = lasso.mincp, 
                        percents.range = percents.range, lasso.delta.cv.mult = lasso.delta.cv.mult, ncv = ncv)

########################
## Summary of Results ##
########################

out2.FDR.org <- org.simu(m,out2.FDR$out,nsimu,x.star.correlated=x.star.correlated,coef=coef)


## Table A1 : H_{0,j} : beta_{x_j} =0, (q-values based on NOT accounting for diet)
xtable(out2.FDR.org$out[,c(paste("ndw3.delta.",delta,sep=""),
                           paste("ndw2.delta.",delta,sep=""),
                           paste("w3.delta.",delta,sep=""),
                           paste("w2.delta.",delta,sep=""),
                           paste("ndw5.delta.",delta,sep=""),
                           paste("w1.delta.",delta,sep=""),
                           paste("ndw1.delta.",delta,sep="")
                                 )],
             digits=sig.digits)

## Table A2 : H_{0,j} : beta_{x_j | z} =0, (q-values based on accounting for diet)

## Thresholding Benjamini-Hochberg adjusted p-values results
xtable(out2.FDR.org$out[,paste("benhoch.alpha.",alpha.bh.range,sep="")],digits=sig.digits)


## Thresholding q-values results
xtable(out2.FDR.org$out[, paste("parcor.alpha.",alpha.range,sep="")],digits=sig.digits)

# q-weighted Lasso results
xtable(out2.FDR.org$out[, paste("w5.delta.",delta.range,sep="")], digits=sig.digits)

## rho-weighted Lasso results (partial correlation weights)
xtable(out2.FDR.org$out[,c(paste("w6.delta.",delta.range,sep=""))],
                                         digits=sig.digits)


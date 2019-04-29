###########
# Library #
###########
library(stats)

###################
# Cases for paper #
###################

#source("../main_codes/RINGARCH_PQKCOV_UNRES.R")
source("RINGARCH_PQKCOV_UNRES.R")

# Situation 1: NBINGARCH(1,1) low mean low heterogeneity data (p=1,q=1)  

k <- log(10)		# this is really k.tilde
n <- 500		
gamma <- -0.1
theta <- c(k,0.5,0.05,gamma,0.1)
p <- 1
q <- 1
simu <- 100
covariates <- arima.sim(model=list(ar=c(1, -0.1)),n=500) #generate AR variable
out.sit1 <- simulation.study2(theta,p,q,n,simu,covariates)
simulation.output(theta,p,q,out.sit1)

# July 14, 2011
$simu.summary
          log(k)      alpha        alpha        gamma        beta     log(k)
mean   2.6068170 0.52803221  0.046167733 -0.106615921  0.06234213 2.36884088
bias   0.3042319 0.02803221 -0.003832267 -0.006615921 -0.03765787 0.06625579
emp.se 4.8425709 0.15820738  0.019994574  0.028801816  0.24654949 0.47527041
est.se 0.0000000 0.00000000  0.000000000  0.000000000  0.00000000 0.59945864
            alpha        alpha        gamma        beta
mean   0.52902191  0.048277295 -0.106697783  0.05565825
bias   0.02902191 -0.001722705 -0.006697783 -0.04434175
emp.se 0.14293456  0.018499411  0.025927033  0.22617552
est.se 0.09496969  0.019881900  0.025147938  0.15646459

$conv.summary
[1] 0 0



# Situation 2 : NBINGARCH(1,1) high mean high heterogeneity data (p=1,q=1) 

k <- log(5)			# this is really k.tilde
n <- 500		
gamma <- 0.3
theta <- c(k,4,-0.05,gamma,0.3)
p <- 1
q <- 1
simu <- 200
covariates <- arima.sim(model=list(ar=c(1, -0.1)),n=500) #generate AR variable
out.sit2 <- simulation.study2(theta,p,q,n,simu,covariates)
simulation.output(theta,p,q,out.sit2)


$simu.summary
          log(k)        alpha        alpha      gamma        beta      log(k)
mean   2.0564202 4.0004718177 -0.051371774 0.30478687 0.308353714 -0.16024493
bias   0.4469823 0.0004718177 -0.001371774 0.00478687 0.008353714 -1.76968285
emp.se 6.7151010 0.1602108883  0.006453238 0.02868007 0.033357092  0.07209991
est.se 0.0000000 0.0000000000  0.000000000 0.00000000 0.000000000  0.12027590
              alpha         alpha        gamma         beta
mean   4.000048e+00 -5.006000e-02 0.3006751791 0.3005596508
bias   4.778084e-05 -5.999988e-05 0.0006751791 0.0005596508
emp.se 5.262290e-02  1.514769e-03 0.0093730366 0.0144972606
est.se 1.141141e-01  1.295877e-03 0.0011747062 0.0281116178

$conv.summary
[1] 2 0


# Situation 3: NBINGARCH(2,1) low mean low heterogeneity data (p=2,q=1) 		

k <- log(10)		# this is really k.tilde
n <- 500		
gamma <- -0.1
theta <- c(k,0.5,0.05,0.01,gamma,0.1)
p <- 2
q <- 1
simu <- 100
covariates <- arima.sim(model=list(ar=c(1,-0.1)),n=500) #generate AR variable
out.sit3 <- simulation.study2(theta,p,q,n,simu,covariates)
simulation.output(theta,p,q,out.sit3)

$simu.summary
          log(k)      alpha       alpha        alpha        gamma         beta
mean   2.8066867 0.52498817  0.04438428  0.003725016 -0.104611095  0.095112291
bias   0.5041016 0.02498817 -0.00561572 -0.006274984 -0.004611095 -0.004887709
emp.se 3.9419949 0.15083508  0.02176061  0.028927055  0.029231927  0.272941293
est.se 0.0000000 0.00000000  0.00000000  0.000000000  0.000000000  0.000000000
          log(k)      alpha        alpha        alpha        gamma         beta
mean   2.4323928 0.51927210  0.046440679  0.006413903 -0.102855925  0.091211267
bias   0.1298077 0.01927210 -0.003559321 -0.003586097 -0.002855925 -0.008788733
emp.se 0.4168965 0.13114254  0.020353378  0.026046734  0.024642333  0.232638274
est.se 0.4575930 0.08291046  0.020265955  0.024671561  0.022698290  0.147593635

$conv.summary
[1] 28  0





# Situation 4: NBINGARCH(2,1) high mean high heterogeneity data (p=2,q=1) 

k <- log(5)			# this is really k.tilde
n <- 500		
gamma <- 0.3
theta <- c(k,4,-0.05,-0.01,gamma,0.3)
p <- 2
q <- 1
simu <- 100
covariates <- arima.sim(model=list(ar=c(-0.9, -0.1)),n=500) #generate AR variable
out.sit4 <- simulation.study2(theta,p,q,n,simu,covariates)
simulation.output(theta,p,q,out.sit4)

$simu.summary
           log(k)       alpha         alpha         alpha        gamma
mean    2.0296912  3.97421954 -0.0504687648 -0.0097121083 0.3005782694
bias    0.4202533 -0.02578046 -0.0004687648  0.0002878917 0.0005782694
emp.se 10.6744765  0.20192245  0.0058009355  0.0026383970 0.0504357149
est.se  0.0000000  0.00000000  0.0000000000  0.0000000000 0.0000000000
              beta      log(k)      alpha         alpha         alpha
mean   0.308642406  0.03539028 4.05336784 -0.0498866946 -0.0105958534
bias   0.008642406 -1.57404763 0.05336784  0.0001133054 -0.0005958534
emp.se 0.054285397  0.10075831 0.22469197  0.0018075830  0.0025310003
est.se 0.000000000  0.11142349 0.29352288  0.0016482499  0.0029692793
              gamma        beta
mean   0.3008422641  0.28674249
bias   0.0008422641 -0.01325751
emp.se 0.0199264219  0.05194654
est.se 0.0012909513  0.05978634

$conv.summary
[1] 63  0



###############
# source file #
###############

source("RINGARCH_PQKCOV_UNRES.R")

##Poisson(1,1) data
m<-1
while(m <=100){
	
	set.seed(1)
	
	N=500 				#sample size of dataset
	lamda=c(rep(0,N)) 		#initialize lambda vector
	x=c(rep(0,N)) 			#initialize the independt variabe vector
	ar=arima.sim(model=list(ar=c(1, -0.1)),n=N) 	# covariates
	#ar <- rnorm(N)
	a0=1
	a1=0.5
	gamma= -0.01
	beta1=0.1
	
	
	lamda[1]=a0+gamma*ar[1]
	x[1]=rpois(1,lambda=lamda[1])
	for (i in 2:N) {
    		lamda[i]=a0+a1*x[i-1]+gamma*ar[i]+beta1*lamda[i-1]
    		x[i]=rpois(1,lambda=lamda[i])
	}

	data <- x
	covariates <- ar

	# setup
	k <- 2
	theta0 <- c(log(k),a0,a1,a1,gamma,beta1)
	out <- optim.resid(data,covariates,p=2,q=1,theta0,iprint=TRUE)
	if(out$results["c.mle",1] >0){
		m <-m+1
	}
	print(m)
	print("conv.results")
	print(out$results["c.mle",1])
}

#pdf("poisson11.pdf")
out <- optim.resid(data,covariates,p=2,q=1,theta0,iprint=TRUE)
#dev.off()
out$results
Box.test(out$resid,lag=12,type="Ljung")	


> out$results
         log(k) alpha alpha  gamma   beta
est.ini   1.478 0.473 0.161 -0.028 -0.011
est       3.991 0.378 0.177 -0.025  0.031
se        1.494 0.487 0.026  0.195  1.220
ratio     2.671 0.778 6.719  0.126  0.025
loglik  -53.392    NA    NA     NA     NA
AIC     116.784    NA    NA     NA     NA
BIC     137.847    NA    NA     NA     NA
c.ini     0.000    NA    NA     NA     NA
c.mle     0.000    NA    NA     NA     NA

> Box.test(out$resid,lag=12,type="Ljung")

        Box-Ljung test

data:  out$resid 
X-squared = 9.6694, df = 12, p-value = 0.645

###############
# source file #
###############

source("RINGARCH_PQKCOV_UNRES.R")

##Independent NB data
set.seed(2)
k=1
a0=1
beta1=1
x1=rnorm(500,0,1)
err=rgamma(500,k,k)
lamda=exp(a0+beta1*x1)*err
y=rpois(500,lamda)

data <- y
covariates <- x1
k <- 1
p <- 1
q <- 0
gamma <- beta1
theta0 <- c(log(k),a0,a0,gamma)

# check if data are independent, ACF, and PACF, they are!
acf(data)
pacf(data)
 
# run NBINGARCH model anyway....

#pdf("independent_nb_k1.pdf")
out <- optim.resid(data,covariates,p,q,theta0,iprint="TRUE")
out <- optim.resid(data,covariates,p,q,theta0,file="_NBIN",iprint=FALSE)
#dev.off()
out$results
Box.test(out$resid,lag=12,type="Ljung")	


> Box.test(out$resid,lag=12,type="Ljung")

        Box-Ljung test

data:  out$resid 
X-squared = 10.8813, df = 12, p-value = 0.5391

# find k and var(k)
kappa <- out$results[2,1]
k <- exp(kappa)
var.kappa <- (out$results[3,1])^2
var.k <- var.kappa * exp(kappa)^2
ratio <- k / sqrt(var.k)
round(k,3)
round(sqrt(var.k),3)
round(ratio,3)

> round(k,3)
[1] 0.945
> round(sqrt(var.k),3)
[1] 0.084

















# Old code



# Another model

source("RINGARCH_PQKCOV_UNRES.R")

##NBINGARCH(1,1) low mean low heterogeneity data (p=1,q=1)
set.seed(1)

N=500 #sample size of dataset
lamda=c(rep(0,N)) #initialize the mean vector
x=c(rep(0,N)) #initialize the independt variabe vector
ar=arima.sim(model=list(ar=c(1, -0.1)),n=500) #generate AR variable
a0=0.5
a1=0.05
theta1= -0.1
beta1=0.1
k=10
lamda[1]=a0+theta1*ar[1]
x[1]=rnbinom(1,mu=exp(lamda[1]),size=k)
for (i in 2:N) {
    lamda[i]=a0+a1*x[i-1]+theta1*ar[i]+beta1*lamda[i-1]
    x[i]=rnbinom(1,mu=exp(lamda[i]),size=k)
}

data <- x
covariates <- ar
k <- log(10)		# this is really k.tilde
n <- 500		
gamma <- -0.1
theta0 <- c(k,0.5,0.05,gamma,0.1)
p <- 1
q <- 1
out <- optim.resid(data,covariates,p,q,theta0,iprint=TRUE)
out$results

# Analyze it ignoring correlation of data 
library(MASS)
nb1=glm.nb(x~ar)
summary(nb1)
plot(resid(nb1))


# Third model

##NBINGARCH(2,1)
N=500 #sample size of dataset
lamda=c(rep(0,N)) #initialize the mean vector
x=c(rep(0,N)) #initialize the independt variabe vector
ar=arima.sim(model=list(ar=c(-0.9, -0.1)),n=500) #generate AR variable
a0=0.1
a1=-0.5
a2=0.1
theta1= 0.03
beta1=0.3
k=5

lamda[1]=a0+theta1*ar[1]
lamda[2]=a0+theta1*ar[2]
x[1]=rnbinom(1,mu=exp(lamda[1]),size=k)
x[2]=rnbinom(1,mu=exp(lamda[2]),size=k)
for (i in 3:N) {
    lamda[i]=a0+a1*x[i-1]+a2*x[i-2]+theta1*ar[i]+beta1*lamda[i-1]
    x[i]=rnbinom(1,mu=exp(lamda[i]),size=k)
}

exp(lamda)
x
 

library(MASS)
nb1=glm.nb(x~ar)
summary(nb1)
plot(resid(nb1))
# Metropolis-Hastings-within-Gibbs algorithm 
# for Weibull posterior simulation

# Data
X <- c(103, 157 , 39 ,145 , 24  ,22 ,122, 126 , 66 , 97)


# Value of lambda1
lambda1 <- 1

# Other hyperparameters
m <- 2
eta0 <- 100

# Euler constant (0.57...)
euler <- -digamma(1)

# Prior log-density for beta (up to a given constant)
log.prior.beta <- function(x)
{
	res <- -(lambda1+1)*log(x) - lambda1*euler/x -lambda1*lgamma(1+1/x)
	return(res)
}

# Prior log-density for eta
log.prior.eta <- function(x)
{
	a <- m
	b <- m/eta0
	res <- a*log(b)-lgamma(a) +(a-1)*log(x) - b*x
	return(res)
}

# Weibull log-likelihood
log.lik <- function(eta,beta)
{
	n <- length(X)
	d <- length(eta)
	GX <- matrix(X,n,d) 
	Geta <- t(matrix(eta,d,n))
	Gbeta <- t(matrix(beta,d,n))
	res <- n*(log(beta)-beta*log(eta)) + (beta-1)*log(sum(X)) - apply( ((GX/Geta)^(Gbeta)),2,sum,na.rm=T)
	return(res)  
}


# Instrumental log-density for eta or beta (random walk on the real line (log(eta) or log(beta)) with constant standard deviation sd.dev)
log.instr <- function(x,y,sd.dev)
{
	return(log(dnorm(log(x),mean=log(y),sd=sd.dev)))
}

# Initialization for 3 parallel MCMC chains
eta.hist <- eta <- runif(3,50,200)
beta.hist <- beta <- runif(3,1,5)
nb.iter.total <- 5000 # number of total MCMC iterations
k <- 0

# Initializing graphical windows
layout(matrix(c(1:2),1,2))
layout.show(2)

# Gibbs simulation loop for 3 parallel MCMC chains
while (TRUE)
{
	
	# Metropolis-Hastings sampling for eta
	sd.eta <- 10
	eta.instr <- exp(rnorm(3,log(eta),sd.eta)) # log-normal sampling 
	u <- runif(3) # uniform sampling
	log.alpha.1 <- log.prior.eta(eta.instr) - log.prior.eta(eta) + log.lik(eta.instr,beta) -  log.lik(eta,beta)
	log.alpha.2 <- log.instr(eta,eta.instr,sd.eta) - log.instr(eta.instr,eta,sd.eta) # this is useless since always 0 because of symetry
	log.alpha.eta <- log.alpha.1 + log.alpha.2
	test <- log(u)  <= log.alpha.eta
	
	eta[test] <- eta.instr[test] # updating 
	eta.hist <- rbind(eta.hist,eta) # history of eta values

	# Metropolis-Hastings sampling for beta
	sd.beta <- 1
	beta.instr <- exp(rnorm(3,log(beta),sd.beta)) # log-normal sampling 
	u <- runif(3) # uniform sampling
	log.alpha.1 <- log.prior.beta(beta.instr) - log.prior.beta(beta) + log.lik(eta,beta.instr) -  log.lik(eta,beta)
	log.alpha.2 <- log.instr(beta,beta.instr,sd.beta) - log.instr(beta.instr,beta,sd.beta) # this is useless since always 0 because of symetry
	log.alpha.beta <- log.alpha.1 + log.alpha.2 
	test <- log(u)  <= log.alpha.beta

	beta[test] <- beta.instr[test] # updating 
	beta.hist <- rbind(beta.hist,beta) # history of beta values



   # Dislaying how the chains change
   matplot(eta.hist,type="l",xlab="nb. iterations",ylab="eta")
   matplot(beta.hist,type="l",xlab="nb. iterations",ylab="beta")

    k <- k+1 
	if (k>=nb.iter.total){break}
}

# Here you should remove autocorrelation 
# to compute the posterior predictive mean

posterior.mean <- median(as.double(eta.hist*gamma(1+1/beta.hist)),na.rm=T)
print(paste("posterior mean is approximately equal to ",posterior.mean))
    
  


# Acceptation-rejection for sampling the prior distribution of beta
# (Exercise on maximum entropy)

# Number of wished prior samples of beta
M <- 100000

# Value of lambda1
lambda1 <- 1

# Euler constant (0.57...)
euler <- -digamma(1)

# Target log-density to sample from (log prior)
log.prior <- function(x)
{
	res <- -(lambda1+1)*log(x) - lambda1*euler/x -lambda1*lgamma(1+1/x)
	return(res)
}

# Instrumental log-density (inverse gamma)
log.instr <- function(x)
{
	lambda1*log(lambda1/euler) - lgamma(lambda1) -(lambda1+1)*log(x) - lambda1*euler/x
}

# Log-constant K
log.K <- euler*log(3/sqrt(pi)) + lgamma(euler)


# Simulation loop (acceptation-rejection)

k <- 0
beta.ok <- beta.non.ok <- c()
while (TRUE)
{
	beta.prop <- 1/rgamma(1,lambda1,lambda1/euler) # instrumental sampling
	u         <- runif(1)						   # uniform sampling
	test      <- log(u) <= log.prior(beta.prop) - log.K - log.instr(beta.prop) # acceptatin-rejection test
	beta.ok   <- c(beta.ok,beta.prop[test])
	beta.non.ok   <- c(beta.non.ok,beta.prop[!test])
	k <- k+1
	if (k>=M){break}
}

# Calculation of the proportion of acceptance
p <- length(beta.ok)/(length(beta.ok)+length(beta.non.ok))

# Calculation of the integration constant
int.const <- (exp(log.K))*p



# Full target (prior) density for the reparametrization Y = 1/beta
target <- function(x)
{
	return((int.const)*exp(log.prior(1/x))/x^2)
}


# Plotting the histogram of accepted samples Y=1/beta and the target density
hist(1/beta.ok,freq=F,xlim=c(0,10),ylim=c(0,0.7))
curve(target,col=2,lwd=2,add=T)  
  
    
  


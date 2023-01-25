#================== Calcul bayésien par MCMC pour une vraisemblance de Gumbel ==========
#                        et un prior obtenu par maximum d'entropie
#======================================================================================== 


X <- c(107.6,
                     72.4,
                     204.5,
                      83.8,
                      142,
                       95.5,
                        316.1,
                        177.9,
                         87.3,
                         81.9,
                         109.1,
                          89.5,
                          150.7,
                          122.1,
                           98.2,
                           113.2,
                           104.4,
                            66.9,
                            136.4,
                             275.4,
                             125,
                             199.8, 
                             51.2,
                              75,
                             168.2,
                              106,
                               72.8)
   
n <- length(X) 
 
# nombre de chaines MCMC
nb.chains <- 3  
 
# nombre d'itération maximal des MCMC
M <- 1000
   
 # log-vraisemblance de Gumbel 
 # (attention, toutes les fonctions de mu et sigma acceptent des vecteurs pour chacune des composantes)
 log.vrais <-  function(mu, sigma)
 {
 	groX <- matrix(X,n,nb.chains)
 	groMu <- t(matrix(mu,nb.chains,n))
 	groSigma <- t(matrix(sigma,nb.chains,n))
 	
 	res <- (- sum(X) + mu)/sigma  -  apply(exp((groX-groMu)/groSigma),2,sum)  
        return(res)
   }     
                               
 # log-densité a priori
  log.dens.prior <- function(mu,sigma)
  {
  	lambda2.tilde <- 0.3155
  	lambda1.tilde <- 0.013155
  	lambda2  <- 1 + lambda2.tilde
  	lambda1 <- - lambda1.tilde
  	euler <- -digamma(1)  # constante d'Euler
  	res <-  (lambda2-2)*log(sigma) + lambda1*sigma*euler + lambda1*mu
  	return(res)
  }

# Lois instrumentales choisies comme des normales calibrées par l'itération précédente (marches aléatoires)

# log-densité d'une loi instrumentale pour un paramètre générique "param"
  cv <- 20/100
  log.dens.instr <- function(x,param.past)
  {
    sdev <- cv*abs(param.past)
    moy <- param.past
    res <- -log(sdev) - (x-param.past)^2/sdev^2
    return(res)
  }

  # Initialisation des chaines MCMC
  mu.past <- mu <- runif(nb.chains,0,100)
  sigma.past <- sigma <- runif(nb.chains,10,30)
  
  # Boucle MCMC
  k = 0
  while (k<M)
  {
    k <- k+1
    
    # Test et mise à jour pour mu -----------------
    mu.instr <- rnorm(nb.chains,mean=mu,sd=cv*abs(mu))   # tirage instrumental
    log.ratio <- log.vrais(mu.instr,sigma) - log.vrais(mu,sigma) +  log.dens.prior(mu.instr,sigma) - log.dens.prior(mu,sigma) + log.dens.instr(mu,mu.instr) - log.dens.instr(mu.instr,mu)   # log-rapport de Hastings-Metropolis
    test <- log(runif(nb.chains)) <=  log.ratio    # test de Hastings-Metropolis
    mu.instr[!test] <- mu[!test]
    mu <- mu.instr
    mu.past <- cbind(mu.past,mu)
    
    # Test et mise à jour pour sigma -----------------
    sigma.instr <- rnorm(nb.chains,mean=sigma,sd=cv*abs(sigma))   # tirage instrumental
    log.ratio <- log.vrais(mu,sigma.instr) - log.vrais(mu,sigma) +  log.dens.prior(mu,sigma.instr) - log.dens.prior(mu,sigma) + log.dens.instr(sigma,sigma.instr) - log.dens.instr(sigma.instr,sigma)   # log-rapport de Hastings-Metropolis
    test <- log(runif(nb.chains)) <=  log.ratio    # test de Hastings-Metropolis
    sigma.instr[!test] <-sigma[!test]
   sigma <- sigma.instr
  sigma.past <- cbind(sigma.past,sigma)
    
    
    # Representation dynamique (reprendre les anciens TPs)
    
  
  }
  
  
  # Representation finale
  iters <- matrix(rep(c(0:M),nb.chains),(M+1),nb.chains)
  layout(matrix(c(1:2),2,1))
  layout.show(2)
  matplot(iters,t(mu.past),type="l",xlab="iterations",ylab="MU")
  matplot(iters,t(sigma.past),type="l",xlab="iterations",ylab="SIGMA")
  
  # Sélection de la période de chauffe
  chauffe <- seq(floor(M/2),M)
  
  # Echantillon approximativement dans la loi a posteriori
  mu.POST <- mu.past[,chauffe]
  sigma.POST <- sigma.past[,chauffe]
  
  # Déccorélation brutale (utiliser autocorr pour faire plus fin)
  pas.decorr <- 5
  MM <- dim(mu.POST)[2]
  steps <- seq(1,MM,by=5)
  mu.POST.decorr <- as.double(mu.POST[,steps])
  sigma.POST.decorr <- as.double(sigma.POST[,steps])
  
  # Nuage de valeurs a posteriori
  plot(mu.POST.decorr,sigma.POST.decorr)
  
 
  
  
  
  
  
  
    
  


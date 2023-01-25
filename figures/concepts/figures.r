#------------- Représentations utiles pour le chapitre 1 -------------

option = 3


if (option==1)
{
#------------ histogramme / densité / fonction de répartition (version francaise)
jpeg("chapitre1-fig1-1.jpeg")
layout(matrix(c(1:2),1,2))
f1 <- function(x){pnorm(x)}
curve(f1,from=(-3),to=3,xlab="valeurs de X",ylab="probabilité",axes=FALSE,lwd=2,col=4)
axis(1);axis(2)
dev.off()
jpeg("chapitre1-fig1-2.jpeg")
M <- 5000
X <- rnorm(M,0,1)
X <- cbind(X,rnorm(M,2,1))
U <- runif(M)
Y <- c(X[U<0.3,1],X[U>=0.3,2])
xo <- c(-3,5)
yo <- c(0,0.4)
hist(Y,freq=F,30,col=gray(0.8),xlim=xo,ylim=yo,axes=FALSE,main="",ylab="",xlab="")
par(new=T)
plot(density(Y),xlab="valeurs de X",ylab="fréquences",xlim=xo,ylim=yo,main="",axes=FALSE,lwd=2,col=2)
axis(1);axis(2,labels=F)
dev.off()


#------------ histogramme / densité / fonction de répartition (version anglaise)
jpeg("chapitre1-fig1-EN-1.jpeg")
layout(matrix(c(1:2),1,2))
f1 <- function(x){pnorm(x)}
curve(f1,from=(-3),to=3,xlab="X",ylab="probability",axes=FALSE,lwd=2,col=4)
axis(1);axis(2)
dev.off()
jpeg("chapitre1-fig1-EN-2.jpeg")
M <- 5000
X <- rnorm(M,0,1)
X <- cbind(X,rnorm(M,2,1))
U <- runif(M)
Y <- c(X[U<0.3,1],X[U>=0.3,2])
xo <- c(-3,5)
yo <- c(0,0.4)
hist(Y,freq=F,30,col=gray(0.8),xlim=xo,ylim=yo,axes=FALSE,main="",ylab="",xlab="")
par(new=T)
plot(density(Y),xlab="X",ylab="Frequency",xlim=xo,ylim=yo,main="",axes=FALSE,lwd=2,col=2)
axis(1);axis(2,labels=F)
dev.off()


}

if (option==2)
{
  #========================== FRANCAIS
# Représentation de quelques densités  
jpeg("chapitre1-fig2-1.jpeg")
 # layout(matrix(c(1:2),1,2))
#layout.show(2)

# Loi binomiale    
f1 <- function(x){dbinom(x,100,0.1)}

# Loi de Poisson
f2 <- function(x){dpois(x,15)}
xx <- c(0:30)
plot(xx,f1(xx),"l",col=2,lty=1,lwd=2,xlab="",ylab="",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
par(new=T)
plot(xx,f1(xx),col=2,lty=1,lwd=2,xlab="",ylab="",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
par(new=T)
plot(xx,f2(xx),"l",col=4,lty=2,pch=2,add=T,lwd=2,xlab="valeurs de X",ylab="densité",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
axis(1);axis(2,labels=F)
par(new=T)
plot(xx,f2(xx),col=4,lty=2,pch=2,add=T,lwd=2,xlab="valeurs de X",ylab="densité",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
axis(1);axis(2,labels=F)
dev.off()

jpeg("chapitre1-fig2-2.jpeg")
# Loi normale
f3 <- function(x){dnorm(x,0,1)}

# Loi du chi2
f4 <- function(x){dchisq(x,1)}

curve(f3,from=(-3),to=5,col=1,lty=1,lwd=2,ylim=c(0,1.2),xlab="valeurs de X",ylab="densité",main="",axes=F)
curve(f4,add=T,col=2,lty=2,lwd=2,axes=F)
axis(1);axis(2,labels=F)
dev.off()
  

#========================== ANGLAIS

# Représentation de quelques densités  
jpeg("chapitre1-fig2-1-EN.jpeg")
# layout(matrix(c(1:2),1,2))
#layout.show(2)

# Loi binomiale    
f1 <- function(x){dbinom(x,100,0.1)}

# Loi de Poisson
f2 <- function(x){dpois(x,15)}
xx <- c(0:30)
plot(xx,f1(xx),"l",col=2,lty=1,lwd=2,xlab="",ylab="",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
par(new=T)
plot(xx,f1(xx),col=2,lty=1,lwd=2,xlab="",ylab="",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
par(new=T)
plot(xx,f2(xx),"l",col=4,lty=2,pch=2,add=T,lwd=2,xlab="X",ylab="density",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
par(new=T)
plot(xx,f2(xx),col=4,lty=2,pch=2,add=T,lwd=2,xlab="X",ylab="density",main="",axes=F,xlim=c(0,30),ylim=c(0,0.17))
axis(1);axis(2,labels=F)
dev.off()

jpeg("chapitre1-fig2-2-EN.jpeg")
# Loi normale
f3 <- function(x){dnorm(x,0,1)}

# Loi du chi2
f4 <- function(x){dchisq(x,1)}

curve(f3,from=(-3),to=5,col=1,lty=1,lwd=2,ylim=c(0,1.2),xlab="X",ylab="density",main="",axes=F)
curve(f4,add=T,col=2,lty=2,lwd=2,axes=F)
axis(1);axis(2,labels=F)
dev.off()





}  
  
if (option==3)
{

 X <- X.last <- 100
  for (i in c(1:100))
 {
    X.old <- X  
   X.last <- X.last + rnorm(1,X.last/500,1) 
   X <- c(X.old,X.last)
  }
 # Représentation d'une chaîne de Markov non-stationnaire
 jpeg("MC-non-stationnaire.jpeg")
 plot(c(0:100),X,xlab="itérations",ylab="valeurs de X",main="",ylim=c(90,150))
dev.off()  
jpeg("MC-non-stationary.jpeg")
plot(c(0:100),X,xlab="iterations",ylab="X",main="",ylim=c(90,150))
dev.off()  
}  
  
  


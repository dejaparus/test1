curve(dnorm(x,2,0.75),-4,4,col=2,lwd=2,lty=2,xlab="theta",main="",axes=FALSE,ylab="")
curve(dnorm(x),add=T,col=1,lwd=2,xlab="theta",main="",axes=FALSE,ylab="")
curve(dcauchy(x,),-4,4,add=T,lwd=2,lty=2,col=4,xlab="theta",main="",axes=FALSE,ylab="")
axis(1):axis(2,labels=FALSE)
legend(-3,0.5,c("posterior","rho 1","rho 2"),col=c(1,2,4),lty=c(1,2,2),lwd=c(2,2,2),bty="n")


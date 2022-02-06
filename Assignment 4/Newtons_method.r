## Newtons_method.R
## Stefano Allesina and Sarah Cobey
## Code for "Theoretical Ecology"


Funct<-function(x,R0){
  # return f(x)
	return(1-exp(-R0*x)-x)
}
DeriFunct<-function(x,R0){
  # return f'(x)
	return(-R0*x*exp(-R0*x)-1)
}

## Implementation to plot fraction recovered
RInf<-function(R0,Funct,DeriFunct){
	xn<-1000000 # a very large number
	xnp1<-0.5 # a guesstimate of the root
	Eps<-10^{-3} # tolerance
	while(abs(xn-xnp1)>Eps){
		xn<-xnp1
		xnp1<-xn-Funct(xn,R0)/DeriFunct(xn,R0)
	}
	return(xnp1)
}

R0s<-seq(0,5,length.out=500)
RecFracs<-rep(0,500)
for (i in 1:500){
    RecFracs[i]<-RInf(R0s[i],Funct,DeriFunct)
}
HerdImmunity <- 1-1/R0s
plot(RecFracs~R0s,t="l",xlab="R0",ylab="Fraction Recovered")
lines(HerdImmunity~R0s, col="blue")
abline(h=1,col="red")
legend(x = "bottomright",          # Position
                legend = c("R(inf)", "herd immunity"),  # Legend texts
                col = c("black", "blue"),           # Line colors
                lwd = 2) 

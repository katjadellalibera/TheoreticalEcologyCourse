## Newtons_method.R
## Stefano Allesina and Sarah Cobey
## Code for "Theoretical Ecology"

## Prototype example of the method
xn<-1000000 # a very large number
xnp1<-1 # a guesstimate of the root
Eps<-10^(-3) # tolerance
while(abs(xn-xnp1)>Eps){
	xn<-xnp1
	xnp1<-xn-Funct(xn,par1,par2)/DeriFunct(xn,par1,par2)
}

## Template for you to change
Funct<-function(x,R0){
	return(???????)
}
DeriFunct<-function(x,R0){
	return(???????)
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
plot(RecFracs~R0s,t="l",xlab="R0",ylab="Fraction Recovered")
abline(h=1,col="red")
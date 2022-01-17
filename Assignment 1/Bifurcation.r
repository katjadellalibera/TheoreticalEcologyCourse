## Bifurcation.R 
## Stefano Allesina sallesina@uchicago.edu allesinalab.uchicago.edu
## Code for "Theoretical Ecology"


## This function returns the values of the min and max 
peaks <- function(x) {
  if (min(x)==max(x)) return(min(x)) ## Does not oscillate 
  l <- length(x) 
  xm1 <- c(x[-1], x[l])
  xp1 <- c(x[1], x[-l]) 
  z<-x[x > xm1 & x > xp1 | x < xm1 & x < xp1] 
  if (length(z)==0) return(min(x)) ## It has not converged yet 
  return (z)
} 

## This function creates a simulation of the logistic map 
LogisticMap<-function(N0,r,TimeSteps){
  Results<-rep(0,TimeSteps) 
  Results[1]<-N0 
  for (j in 2:TimeSteps){
    Results[j]<-r*Results[j-1]*(1-Results[j-1])
  } 
  return(Results)
} 

## Plot the Diagram 
plot(0,0, xlim=c(0,4), ylim=c(-0.05,1.05),type="n", xlab="r", ylab="X") 
for (r in seq(0.001,4,0.005)) { # These are the initial and final values for r
  out <- LogisticMap(0.5,r,2500) # Initial conditions 
  l <- length(out) %/% 10 # use only the last 250 steps 
  out <- out[(9*l):(10*l)] 
  p <- peaks(out) 
  l <- length(out) 
  points(rep(r, length(p)), p, pch=".")
}

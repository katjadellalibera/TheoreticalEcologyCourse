## ExpGrowth.R 
## Stefano Allesina sallesina@uchicago.edu allesinalab.uchicago.edu
## Code for "Theoretical Ecology"

require(deSolve) 
## package for integrating numerically ODEs see 
##http://cran.r-project.org/web/packages/deSolve/

## This function takes a step in time
ExponentialGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  dX<-r*X 
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunExponentialGrowth<-function(MaxTime=10,GrowthRate=0.1,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = ExponentialGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

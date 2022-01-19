## ExpGrowth.R 
## Stefano Allesina sallesina@uchicago.edu allesinalab.uchicago.edu
## Code for "Theoretical Ecology"

require(deSolve) 
## package for integrating numerically ODEs see 
##http://cran.r-project.org/web/packages/deSolve/

### Exponential Growth ### 

## This function takes a step in time
ExponentialGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  dX<-r*X 
  return(list(dX)) ## for some reason, you have to return a list
}

## This function runs the model and produces the trajectory
RunExponentialGrowth<-function(MaxTime=10,GrowthRate=0.2,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = ExponentialGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
out <- RunExponentialGrowth()
# plot the solution 
ggplot(as.data.frame(out)) +
  aes(x=time, y=X) +
  geom_line() +
  ggtitle("Exponential growth") +
  ylab("Density of X") +
  xlab("Time") + 
  theme_bw()


### Logistic Growth ###

## This function takes a step in time
LogisticGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  k <- parameters[2] ## the second parameter is the carrying capacity
  dX<-r*X*(1-X/k) 
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunLogisticGrowth<-function(MaxTime = 10, GrowthRate = 0.5,
                            CarryingCapacity = 5, InitialX = 1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate,k=CarryingCapacity)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = LogisticGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
logisticOut <- RunLogisticGrowth()
# plot the solution 
ggplot(as.data.frame(logisticOut)) +
  aes(x=time, y=X) +
  geom_line() +
  ggtitle("Logistic growth") +
  ylab("Density of X") +
  xlab("Time") + 
  theme_bw()


### Allee Growth ###

## This function takes a step in time
AlleeGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  u <- parameters[1] ## the first parameter is the upper limit of the growth rate
  v <- parameters[2] ## the second parameter is critical point
  c <- parameters[3] ## the third parameter is the density dependent death rate
  dX<-((u*X/(v+X)-c*X)*X)
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunAlleeGrowth<-function(MaxTime = 10, UpperGrowthRate = 2,
                          CriticalPoint = 5, DeathRate = 0.1, InitialX = 1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(u=UpperGrowthRate,v=CriticalPoint, c= DeathRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = AlleeGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
AlleeOut <- RunAlleeGrowth()
# plot the solution 
ggplot(as.data.frame(AlleeOut)) +
  aes(x=time, y=X) +
  geom_line() +
  ggtitle("Allee growth") +
  ylab("Density of X") +
  xlab("Time") + 
  theme_bw()


### Allee Growth ###

## This function takes a step in time
AlleeGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  u <- parameters[1] ## the first parameter is the upper limit of the growth rate
  v <- parameters[2] ## the second parameter is critical point
  c <- parameters[3] ## the third parameter is the density dependent death rate
  dX<-((u*X/(v+X)-c*X)*X)
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunAlleeGrowth<-function(MaxTime = 10, UpperGrowthRate = 2,
                         CriticalPoint = 5, DeathRate = 0.1, InitialX = 1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(u=UpperGrowthRate,v=CriticalPoint, c= DeathRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = AlleeGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
AlleeOut <- RunAlleeGrowth()
# plot the solution 
ggplot(as.data.frame(AlleeOut)) +
  aes(x=time, y=X) +
  geom_line() +
  ggtitle("Allee growth") +
  ylab("Density of X") +
  xlab("Time") + 
  theme_bw()


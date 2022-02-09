library(deSolve)
library(ggplot2)


#constants
parameters <- c(R0 = 18, 
                gamma = 0.2)
N <- 100000

#initial values
state <- c(S = 0.9999,
           R = 0,
           I = 0.0001)

Epidemic <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    dS <- -R0*gamma*I
    dI <- R0*gamma*S*I-gamma*I
    dR <- gamma*I
    
    list (c(dS, dI, dR))
  })
}

times <- seq(0, 30, by = 1)

epidemicOutRadau <- ode(state, times, Epidemic, parameters, method = "radau",
            atol = 1e-4, rtol = 1e-4)
epidemicOutRadau <- data.frame(epidemicOutRadau)
epidemicOutRadau[c("S","I","R")] <- epidemicOutRadau[c("S","I","R")]*N

epidemicOutrk4 <- ode(state, times, Epidemic, parameters, method = "rk4")
epidemicOutrk4 <- data.frame(epidemicOutrk4)
epidemicOutrk4[c("S","I","R")] <- epidemicOutrk4[c("S","I","R")]*N

library(tidyr)
dfRadau <- gather(epidemicOutRadau, key = Category, value = population, 
             c("S", "I", "R"))
dfrk4 <- gather(epidemicOutrk4, key = Category, value = population, 
                c("S", "I", "R"))

ggplot(dfRadau, aes(x=time, y = population, group = Category, colour = Category)) + 
  geom_line() +
  labs(x = "Day",
       y = "Population")

ggplot(dfrk4, aes(x=time, y = population, group = Category, colour = Category)) + 
  geom_line() +
  labs(x = "Day",
       y = "Population")





  
  
  

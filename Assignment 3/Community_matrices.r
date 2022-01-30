## Community_matrices.R
## Stefano Allesina and Sarah Cobey 
## Code for Theoretical Ecology

library(ggplot2)
library(psych)
S <- 250
sigma <- 0.3

MayMatrix<-function(S,C,sigma){
  ## This matrix determines the connections
  A<-matrix(runif(S*S),S,S)
  ## Contains the values for the connections
  B<-matrix(rnorm(S*S,0.0,sigma),S,S)
  A<-(A<= C)*1 # A matrix contains 1 when A[i,j] <= C
  M<-A*B
  diag(M)<- -1
  return(M)
}

PPMatrix<-function(S,C,sigma){
  ## Determine the signs for the connections
  MyS<-sign(rnorm(S*(S-1)/2))
  A<-matrix(0,S,S)
  A[upper.tri(A,diag=F)]<-MyS
  D<-matrix(runif(S*S),S,S)
  D<-(D <= C)*1
  A<-D*A
  A<- A-t(A)
  ## Contains the  values for the connections
  B<-matrix(abs(rnorm(S*S,0.0,sigma)),S,S)
  M<-A*B
  diag(M)<- -1
  return(M)
}

stableQ<- function(m){
  if(tr(m)<0 & det(m)>0){
    return(1)
  }else{
    return(0)
  }
}


C <-seq(0.01,0.5,by=0.01)

mays <- c()
for (cs in C){
  mays <- append(mays,sum(replicate(1000, stableQ(MayMatrix(S,cs,sigma)), simplify = TRUE )))
}
     
pred <- c()
for (cs in C){
  pred <- append(pred,sum(replicate(1000, stableQ(PPMatrix(S,cs,sigma)), simplify = TRUE )))
}

data <- data.frame(C,mays,pred)

ggplot(data, aes(x = C)) +
  geom_line(aes(y = mays), color = "black") +
  geom_line(aes( y = pred), color = "red")

## Community_matrices.R
## Stefano Allesina and Sarah Cobey 
## Code for Theoretical Ecology

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

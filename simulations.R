## __________________________________________________________
##
## class.ind
##
## ----------------------------------------------------------

class.ind<-function (cl) {
  n  <- length(cl);
  cl <- as.factor(cl);
  x  <- matrix(0, n, length(levels(cl)));
  x[(1:n) + n * (unclass(cl) - 1)] <- 1;
  dimnames(x) <- list(names(cl), levels(cl));
  return(x);
}

## __________________________________________________________
##
## graph.affiliation
## warning : could be faster.
##
## INPUT  n: number of vertex
##        alphaVect: vecteur of class proportion
##        lambda: proba of edge given  same classe
##        epsilon: proba of edge given two different classes
## OUTPUT x: adjacency matrix
##        cluster: class vector
## ----------------------------------------------------------

graph.affiliation<-function(n=100,alphaVect=c(1/2,1/2),lambda=0.7,epsilon=0.05) {
   x <- matrix(0,n,n);
   Q <- length(alphaVect);
   rmultinom(1, size=n, prob = alphaVect)->nq;
   Z <- class.ind(rep(1:Q,nq));
   Z <- Z[sample(1:n,n),];
   for (i in 1:n) {
     for (j in i:n) {
     # if i and j in same class
       if (which.max(Z[i,])  == which.max(Z[j,])) p<-lambda else  p<-epsilon
          if ((rbinom(1,1,p))&(i != j)) {
            x[i,j]<-1;
            x[j,i]<-1;
          }
      }
   }
  return(list(x=x,cluster=apply(Z,1,which.max)) )
}


## __________________________________________________________
##
## CreateMu
##
## INPUT: num.classes: Number of classes
##	  Mu1: Mean for the nodes belonging to the same cluster
##	  Mu2: Mean for the nodes belonging to different cluster
##
## OUTPUT: Mu: Matrix of the means for the variable Y
##
## __________________________________________________________

CreateMu <- function(num.classes, Mu1, Mu2) {
  size       <- num.classes;
  Mu         <- matrix( Mu2, nrow=size, ncol=size);
  diag(Mu)   <- Mu1;
  return(Mu);
}

## __________________________________________________________
##
## SimDataYcondZ
##
## INPUT: node.classes: Vector of node class labels
##	 MU: Matrix of the means
##	 Sigma: Standard Error
##
## OUTPUT: Y: Similarity matrix
##
## Simulate a similarity matrix conditionally to Z
## __________________________________________________________

SimDataYcondZ <- function(node.classes, Mu, Sigma, SelfLoop = FALSE) {
  num.nodes <- length (node.classes);
  Y <- matrix(0, num.nodes, num.nodes);
  for (i in 1:num.nodes) {
    for (j in 1:i) {
      Y[i,j] <- rnorm(1,Mu[ node.classes[i], node.classes[j]], Sigma)
      Y[j,i] <- Y[i,j]
    }
    if (SelfLoop == FALSE){
      Y[i,i] <- 0;
    }
  }
  return (Y);
}

## __________________________________________________________
##
## SimDataYcondXZ
##
## INPUT: node.classes: Vector of node class labels
##	  X: Adjency matrix
##	  MuX0: Matrix of the means (case X=0)
##	  MuX1: Matrix of the means (case X=1)
##	  Sigma: Variance
##
## OUTPUT: Y: Similarity matrix
##
## Simulate a similarity matrix conditionally to X and Z
## __________________________________________________________

SimDataYcondXZ <- function(node.classes, X, MuX0, MuX1, Sigma, SelfLoop = FALSE) {
   num.nodes <- length (node.classes);
   Y <- matrix(0, num.nodes, num.nodes);

   for (i in 1:num.nodes) {
     for (j in 1:i) {
       if (X[i,j] ==0 ){
         Y[i,j] <- rnorm(1,MuX0[ node.classes[i], node.classes[j]], Sigma);
         Y[j,i] <- Y[i,j] ;
       } else {
         Y[i,j] <- rnorm(1,MuX1[node.classes[i], node.classes[j]], Sigma);
         Y[j,i] <- Y[i,j]  ;
       }
     }
    if (SelfLoop == FALSE){
      Y[i,i] = 0;
    }
  }
  return (Y);
}

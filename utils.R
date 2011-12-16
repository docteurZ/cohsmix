## __________________________________________________________
##
## SortMatrix
##
## INPUT :	node.classes	-> Vector of the classes
##	   	Matrix		-> Matrix we want to sort according 
##			   	   to the classes
##
## OUTPUT : 	The sorted matrix
##
## Sort a matrix according to the node classes
## __________________________________________________________

SortMatrix<-function(node.classes, Matrix){
  rank <- order(node.classes);
  Matrix[rank,rank] -> Matrix;
  return(Matrix);
}

##------------------------------------------------------------
##
## TauInit 
##
## INPUT : g :cmgraph
##
## OUTPUT : Class Matrix
##
## Create a class matrix : use class.ind
##
TauInit<-function(g) {
  return(class.ind(g$nodeclasses));
}

## __________________________________________________________
##
## getCluster
## 
## __________________________________________________________

getCluster <-function(Tau){
  mincut <- log(.Machine$double.xmin);
  Tau[is.nan(Tau) == TRUE] <- exp(mincut);
  cluster <- apply(Tau,1, which.max);
  return(cluster);
}

## __________________________________________________________
##
## randError
##
## INPUT :
##   x, y : partion assignment vectors
##
## OUTPUT:
##   scalar : adjusted rand index
##
## Compute the adjusted rand index between two partitions
## __________________________________________________________

randError<-function(x, y) {

 # first, get crosstabs
 ctab <- table(x,y);

 # now calculate 4 intermediary sums
 cellsum <- sum(ctab * (ctab-1) / 2);
 totsum  <- sum(ctab) * (sum(ctab)-1) / 2;

 # use matrix multiplication to get row and column marginal sums
 rows   <- ctab %*% rep(1,ncol(ctab));
 rowsum <- sum(rows*(rows-1)/2);
 cols   <- rep(1,nrow(ctab)) %*% ctab;
 colsum <- sum(cols*(cols-1)/2);

 # now put them together
 adj.rand <-  (cellsum - (rowsum*colsum/totsum)) /
              (.5*(rowsum +colsum) - (rowsum*colsum/totsum));

 if ( adj.rand < 0) {
 	adj.rand <- 0;
 }
  return (adj.rand);
}


## __________________________________________________________
##
## ClusteringCoef
##
##
## INPUT :		X            -> Adjency matrix
##			node.classes -> Vector of the classes
##
## Calculation of the clustering coefficient
## __________________________________________________________

ClusteringCoef <- function(X, node.classes){
  prop  <- 0;
  clust <- 0;
  
  for (i in 1:length(node.classes)){
    for (j in 1:length(node.classes)){
      for (k in 1:length(node.classes)){
        if ((X[i,j]==1) && (X[i,k]==1)){
           prop = prop +1;
	}
	if ((X[i,j]==1) && (X[i,k]==1) && (X[j,k]==1)){
          clust = clust+1;
        }
      }
    }
  }
  result = clust/prop;
  return(result);
}


## __________________________________________________________
##
## Modularity
##
##
## INPUT :		graph -> Adjency matrix
##			clus  -> Vector of the classes
##
## Calculation of the modularity
## __________________________________________________________

Modularity <- function(graph, clus) {
  m <- sum(graph);
  maxId <- dim(graph)[1];	
  tmp <- 0;
  out <- apply(graph, 1 , sum)  ;
  for (i in 1:maxId) {
    for (j in 1:maxId) {
     if (clus[i] == clus[j]) {
      tmp <- tmp + (graph[i,j] - (out[i] * out[j])/(2*m));
     }
    }
  }
  return(tmp / (2*m));
}

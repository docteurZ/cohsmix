## __________________________________________________________
##
## RepresentationXGroup
##
## INPUT: x: Adjency matrix
##        cluste: Vector of the classes
##
## Representation of the nodes and the corresponding groups after reorganization
## __________________________________________________________

RepresentationXGroup <- function(x, cluster){
  x <- SortMatrix(cluster, x); # reorganize the adjacency matrix  
  dim(x)[1]->n;
  m<-which(x==1,arr.ind=TRUE);
  plot(1, 1, xlim = c(0, n + 1), ylim = c(n + 1, 0), type = "n", axes= FALSE,xlab="classes",ylab="classes",main="Reorganized Adjacency matrix")
  rect(m[,1]-0.5,m[,2]-0.5,m[,1]+0.5,m[,2]+0.5,col=1);
  table(cluster)->limits; # find the class limits
  cumsum(limits)[1:(length(limits)-1)]+0.5->limits;
  abline(v=c(0.5,limits,n+0.5),h=c(0.5,limits,n+0.5),col="red");
}


## __________________________________________________________
##
## RepresentationXY
##
## INPUT: X: Adjency matrix
##        Y: Similarity matrix
##	  node.classes: Vector of the classes
##	  Sigma: Variance	
##
## Represention of the variables Yij and the affiliation matrix
## __________________________________________________________

RepresentationXY <- function(X, Y, node.classes, DrawGroup=TRUE) {
  OrderX <- SortMatrix(node.classes,X); #sorted X#
  OrderY <- SortMatrix(node.classes,Y); #sorted Y#
  image(OrderY);
  size=length(node.classes);
  Xlie <- which(OrderX==1,arr.ind=TRUE);
  points((Xlie-1)/(size-1), pch=20);		#Dilatation#
  temp <- rep(0, (max(node.classes)[1]+1));
  if (DrawGroup==TRUE){
    for (i in 1:max(node.classes)){
     mq <- table(node.classes);
     size <- length(node.classes);
     temp[i+1] <-  mq[i];
     sum <- sum(temp);
     axe <- sum/size;
     abline(v=axe,h=axe,col="blue");
   }
 }
}

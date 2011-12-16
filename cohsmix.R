source("utils.R");
source("simulations.R");
source("representations.R");
source("estimation_model_1.R");
source("estimation_model_2.R")
source("estimation_model_3.R")
source("icl.R");

##---------------##
## generateGraph ##
##---------------##
generateGraph<-function(nodes=100, classes=2, alphaVect=rep(1/classes,classes), intraproba=0.1, crossproba=0.02) {
  g <- graph.affiliation(n=nodes,alphaVect=alphaVect,lambda=intraproba,epsilon=crossproba);
  cmgraph <- list(nodes=nodes,
                  classes=classes,
                  adjacency=g$x,
                  nodeclasses=g$cluster,
                  alpha=alphaVect,
                  intraproba=intraproba,
                   crossproba=crossproba);
  attr(cmgraph,'class') <- c('cmgraph');
  return(cmgraph);
}

## get labels ##
labels.cmgraph<-function(object,...) {
  c("Nodes","Classes","Adjacency Matrix","Node Classification","Class Probability Distribution","Intra Class Edge Probability","Cross Class Edge Probability");
}

## summary ##
summary.cmgraph<-function(object,...)  {
   cat(c("Nodes                         : ",object$nodes,"\n",
         "Edges                         : ",length(which(object$adjacency!=0)),"\n",
         "Classes                       : ",object$classes,"\n",
         "Class Probability Distribution: ",object$alpha,"\n"));
}

## plot adjacency graph ##
plot.cmgraph<-function(x,...) {
   RepresentationXGroup(x$adjacency,x$nodeclasses);
}

##--------------------##
## generateCovariates ##
##--------------------##
## conditionaly to Z  ##
generateCovariatesCondZ<-function(g,sameclustermean=0,otherclustermean=2,sigma=1) {
   mu=CreateMu(g$classes,sameclustermean,otherclustermean);
   res=SimDataYcondZ(g$nodeclasses,mu,sigma);
   cmcovars=list(graph=g,sameclustermean=sameclustermean,otherclustermean=otherclustermean,sigma=sigma,mu=mu,y=res);
   attr(cmcovars,'class')<-c('cmcovarz','cmcovar');
   return(cmcovars);
}

## conditionaly to Z and X ##
generateCovariatesCondXZ<-function(g,sameclustermean=c(0,3),otherclustermean=c(2,5),sigma=1) {
  mux0=CreateMu(g$classes,sameclustermean[1],otherclustermean[1]);
  mux1=CreateMu(g$classes,sameclustermean[2],otherclustermean[2]);
  res=SimDataYcondXZ(g$nodeclasses,g$adjacency,mux0,mux1,sigma);
  cmcovars=list(graph=g,sameclustermean=sameclustermean,otherclustermean=otherclustermean,sigma=sigma,mu=c(mux0,mux1),y=res);
  attr(cmcovars,'class')<-c('cmcovarxz','cmcovar');
  return(cmcovars);
}

## summary ##
summary.cmcovar<-function(x,...) {
  cat("Classes           : ",x$graph$classes,"\n",
      "Intra cluster mean: ",x$sameclustermean,"\n",
      "Cross cluster mean: ",x$otherclustermean,"\n",
      "Variance          : ",x$sigma,"\n",
      "Covariates       :\n",x$y,"\n")
}

## plot covariates with graph ##
plot.cmcovar<-function(x,...) {
   RepresentationXY(x$graph$adjacency,x$y, x$graph$nodeclasses);
}

##-------------##
## Estimation  ##
##-------------##
## estimateCondZ ##
estimateCondZ<-function(graph,covars,maxiterations,initialclasses,selfloops) {
  res=EMalgorithmZ(initialclasses,covars$y,graph$adjacency,maxiterations,FALSE,selfloops);
  cmestimation=list(mean=res$MuEstimated,variance=res$VarianceEstimated,
                    pi=res$PIEstimated,alpha=res$AlphaEstimated,
                    tau=res$TauEstimated,jexpected=res$EJ,graph=graph);
  attr(cmestimation,'class')<-c('cmestimationz');
  return(cmestimation);
}

privateestimate<-function(covars,graph,maxiterations,initialclasses,selfloops,...) UseMethod("privateestimate")

privateestimate.cmcovarz<-function(covars,graph,maxiterations,initialclasses,selfloops,...) {
   res=estimateCondZ(graph,covars,maxiterations,initialclasses,selfloops);
   attr(res,'class')<-c(attr(res,'class'),'cmestimation');
   return(res);
}

## estimateCondXZ ##
estimateCondXZ<-function(graph,covars,maxiterations,initialclasses,selfloops) {
  res=EMalgorithmXZ(initialclasses,covars$y,graph$adjacency,maxiterations,selfloops);
  cmestimation=list(mean=c(res$MuEstimated1,res$MuEstimated2),variance=res$VarianceEstimated,
                    pi=res$PIEstimated,alpha=res$AlphaEstimated,
                    tau=res$TauEstimated,jexpected=res$EJ,graph=graph);
  attr(cmestimation,'class')<-c('cmestimationxz');
  cmestimation
}

privateestimate.cmcovarxz<-function(covars,graph,maxiterations,initialclasses,selfloops,...) {
   res=estimateCondXZ(graph,covars,maxiterations,initialclasses,selfloops);
   attr(res,'class')<-c(attr(res,'class'),'cmestimation');
   return(res);
}

estimate<-function(graph,covars,...) UseMethod("estimate");

estimate.cmgraph<-function(graph,covars,maxiterations=20,initialclasses=TauInit(graph),selfloops=FALSE,method=NULL,...) {
  if (length(method)  == 0) {
    res=privateestimate(covars,graph,maxiterations,initialclasses,selfloops,...)
  } else {
    res=method(graph,covars,maxiterations,initialclasses,selfloops)
    attr(res,'class')<-c(attr(res,'class'),'cmestimation')
  }
  return(res);
}

plot.cmestimation<-function(x,...) {
  print("plot");
  par(mfrow = c(1,2))
  plot(x$jexpected)
  title("Expected value of J: Convergence criterion")
  clusters=getCluster(x$tau);
  gplot(x$graph$adjacency,vertex.col=clusters+2)
  title("Network with estimated classes")
}

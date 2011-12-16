## __________________________________________________________
##
## EMalgorithmNodes
##
## INPUT: Tau	-> Initial classification matrix
##	  Y		-> Matrix of the covariates
##	  X		-> Adjacency matrix
##        NbIteration	-> Number of iteration
##        SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
##
## OUTPUT: Mu		-> Matrix of the estimation of the mean 
##	   VarCovEstimated	-> Estimation variance-covariance matrix
##	   PIEstimated	-> Estimation of the connectivity matrix
##	   AlphaEstimated	-> Estimation of the probability
##			           for a node i to belong to class q
##         TauEstimated	-> Estimation of the variational paramater
##	   EJ		-> Value of the expected value of J for 
##			   each iteration
##		
## => Plot of EJ curve and network with the estimated classe
## 
## __________________________________________________________

EMalgorithmNodes <-function(Tau,Y, X, NbIteration, Plot=TRUE, SelfLoop=FALSE){
############# Notations ##################
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(X)[1];
   nbCov <- dim(Y)[2];
   AlphaEstimated <- vector(length=nbGroup);
   facteur <- matrix(1/nbNodes,nrow=1,ncol=nbNodes);
   EJ <- 0;
   TauEstimated <- Tau;
   EJnew<-0;
   Mu<-matrix(0,nrow=nbCov,ncol=nbGroup);

##########################################
#	ALGO	                         #
##########################################
   cat("Iteration: ");
   for (i in 1:NbIteration) {
     cat(i," ");
     TauPrim <- t(TauEstimated);
     if (SelfLoop == FALSE){
       div = t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated)-(TauPrim%*%TauEstimated);
     } else {
       div = t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated);
     }		
####### Estimation of the mean Mu for the Q groups#######
    for (q in 1:nbGroup) {
      Mu[,q]<-  colSums(Y*TauEstimated[,q])/sum(TauEstimated[,q])
    }

####### Estimation of PI #######
   PItemp <- (TauPrim%*%X%*%TauEstimated);
   PIEstimated <- PItemp/(div);
   PIEstimated[is.nan(PIEstimated) ==TRUE] <- exp(mincut);	
   PIEstimated[PIEstimated == 'Inf'] <- exp(mincut);
   PIEstimated[PIEstimated == 1] <- (1-exp(mincut));
   PIEstimated[PIEstimated < exp(mincut)] <- exp(mincut);
		
####### Estimation of Alpha #######
   AlphaEstimated <- facteur %*% TauEstimated;
####### Estimation of the variance-covariance matrix #######
   sigma2<-0;
   for (q in 1:nbGroup) {
     sigma2 <- sigma2 + (1/sum(TauEstimated[,q])) * sum(TauEstimated[,q]  * apply((t(Y) - Mu[,q])^2,2, sum));
   }
   sigma2<-sigma2/nbCov;
   VarCovEstimated<-diag(rep(sigma2,nbCov));

####### Expected value of J #######
   EJnew<-ExpectedJNodes(TauEstimated,Y,X,Mu,PIEstimated,AlphaEstimated,VarCovEstimated,SelfLoop=SelfLoop);
   #print("NEw");
   #print(EJnew);
   EJ <- c(EJ,EJnew);
   EJold<-EJ[i];
   #print("old");
   #print(EJold)
   if (EJold!=0){
    critere <- abs((EJnew-EJold)/EJold);
    print(critere);
    if (critere<=0.0000001){
      break();
    }
   }
		
####### Estimation of Tau #######
   Tautemp <- EstimTauNodes(TauEstimated,Y,X,Mu,PIEstimated,AlphaEstimated,VarCovEstimated,SelfLoop=SelfLoop);
   TauEstimated <- Tautemp$TauEstimated;	
  }
  cat("\n");
  
  if (Plot==TRUE){
   #################### Representation #######################
   #x11()
   #par(mfrow = c(1,2))
   #cat("Plotting the convergence criterion...","\n")
   #plot(EJ[2:NbIteration])
   #title("Expected value of J: Convergence criterion")
   #cat("Plotting the estimated network structure...","\n")
   #map=MAP(TauEstimated)
   #gplot(Adjacente,vertex.col=map$node.classes+2)
   #title("Network with estimated classes")
   ############################################################
   }
   return(list(Mu=Mu,VarCovEstimated=VarCovEstimated,
               PIEstimated=PIEstimated,
               AlphaEstimated=AlphaEstimated,
               TauEstimated=TauEstimated,
               EJ=EJ[2:length(EJ)]));
}
## __________________________________________________________
##
## EstimTauNodes
##
## INPUT: Tau		-> Initial classification matrix
##	  Y		-> Matrix of the covariates
##	  X		-> Adjacency matrix
##	  Mu		-> Matrix of the estimations of the 
##			   means
##	  PI		-> Connectivity matrix
##	  Alpha		-> Vector  of the probability
##	 		   for a node i to belong to class q
##	  VarCov		-> Variance-Covariance matrix
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
##
## OUTPUT : 	TauEstimated	-> Estimation of the variational paramater
##				   according to the inputs 
## __________________________________________________________

EstimTauNodes <-function(Tau, Y, X, Mu, PI, Alpha, VarCov, SelfLoop = FALSE){
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(Y)[1];
   TauPrim <- t(Tau);
   TauEstimated <- Tau;
   LogTauEstimated <- matrix(0,nrow=nbNodes,ncol=nbGroup);
   PIPrim <- t(PI);
   H <- 1-PI;
   HPrim <- t(H);
   VarCovInv <- ginv(VarCov);
   detVar <- det(VarCov);
   maxcut <- log(.Machine$double.xmax) - log(nbGroup);
   mincut <- log(.Machine$double.xmin);
   for (q in 1:nbGroup){
     HPrim[q,] <- pmin(HPrim[q, ], exp(maxcut));
     HPrim[q,] <- pmax(HPrim[q, ], exp(mincut));
     PIPrim[q,] <- pmin(PIPrim[q, ], exp(maxcut));
     PIPrim[q,] <- pmax(PIPrim[q, ], exp(mincut));
     Alpha[q] <- pmin(Alpha[q], exp(maxcut));
     Alpha[q] <- pmax(Alpha[q], exp(mincut));	
   }

   HPrim[is.nan(HPrim) == TRUE] <- exp(mincut);
   HPrim[HPrim < exp(mincut) ] <- exp(mincut);
   H[is.nan(H) == TRUE] <- exp(mincut);
   H[H < exp(mincut) ] <- exp(mincut);	

   LogAlpha <- log(Alpha);
   LogPIPrim <- log(PIPrim);
   LogHPrim <-  log(HPrim);
	
   BernMatrix = (((X%*%Tau)%*%(LogPIPrim)) + (((1-X)%*%Tau)%*%(LogHPrim)));
   TauLogHPrim <- Tau%*%LogHPrim;
   for (i in 1:nbNodes){
    for (q in 1:nbGroup){
       Bern <- BernMatrix[i,q];
       Norm <-  log(1/(2*pi^(nbNodes/2)*(detVar^(1/2))))-1/2 *(Y[i,] - Mu[,q])%*%VarCovInv%*%t(t((Y[i,] - Mu[,q])));
       EgaliteIJ <-  - TauLogHPrim[i,q];
       if (SelfLoop == FALSE){ 
          LogTauEstimated[i,q] <-  LogAlpha[q] + Bern + Norm + EgaliteIJ;
       } else {
          LogTauEstimated[i,q] <-  LogAlpha[q] + Bern + Norm;
       }
       if (is.nan(LogTauEstimated[i,q])==TRUE){
          LogTauEstimated[i,q]<-mincut;
       }
       if (LogTauEstimated[i,q]=='Inf'){
         LogTauEstimated[i,q]<-mincut;
       }
    }
    LogTauEstimated[i, ] <- pmin(LogTauEstimated[i, ], maxcut);
    LogTauEstimated[i, ] <- pmax(LogTauEstimated[i, ], mincut);
    TauEstimated[i, ] <- exp(LogTauEstimated[i, ]);
    Normalize <- 1/sum(TauEstimated[i,]);
    TauEstimated[i,] <- TauEstimated[i,] * Normalize;
    TauEstimated[i, ][TauEstimated[i, ] < .Machine$double.xmin] <- .Machine$double.xmin;	
  }
  return( list(TauEstimated=TauEstimated));
}

## __________________________________________________________
##
## ExpectedJNodes
##
## INPUT: Tau		-> Initial classification matrix
##	  Y		-> Matrix of the covariates
##	  X		-> Adjacency matrix
##	  Mu		-> Matrix of the estimations of the 
##			   means
##	  PI		-> Connectivity matrix
##	  Alpha		-> Vector  of the probability
##			   for a node i to belong to class q
##	  VarCov		-> Variance-Covariance matrix
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
##
## OUTPUT : 	Expected	-> Expected value of J
##
## __________________________________________________________

ExpectedJNodes <-function(Tau, Y, X, Mu, PI, Alpha, VarCov, SelfLoop=FALSE){
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(Y)[1];
   
   entropieTemp <- 0;
   PIPrim <- t(PI);
   H <- 1-PI;
   HPrim <- t(H);
   maxcut <- log(.Machine$double.xmax) - log(nbGroup);
   mincut <- log(.Machine$double.xmin);
   VarCovInv <- ginv(VarCov);
   detVar <- det(VarCov);
   Norm<-0;
	
################ Entropie #############################
    for(q in 1:nbGroup){
      HPrim[q,] <- pmin(HPrim[q, ], exp(maxcut));
      HPrim[q,] <- pmax(HPrim[q, ], exp(mincut));
      PIPrim[q,] <- pmin(PIPrim[q, ], exp(maxcut));
      PIPrim[q,] <- pmax(PIPrim[q, ], exp(mincut));
      Alpha[q] <- pmin(Alpha[q], exp(maxcut));
      Alpha[q] <- pmax(Alpha[q], exp(mincut));
      PI[q,] <- pmin(PI[q,], exp(maxcut));
      PI[q,] <- pmax(PI[q,], exp(mincut));
      H[q,] <- pmin(H[q,], exp(maxcut));
      H[q,] <- pmax(H[q,], exp(mincut));
    }
    
    HPrim[is.nan(HPrim) == TRUE] <- exp(mincut);
    HPrim[HPrim < exp(mincut) ] <- exp(mincut);
    H[is.nan(H) == TRUE] <- exp(mincut);
    H[H < exp(mincut) ] <- exp(mincut);	
    for(i in 1:nbNodes){
      Tau[i, ] <- pmin(Tau[i, ], exp(maxcut));
      Tau[i, ] <- pmax(Tau[i, ], exp(mincut));
      for (q in 1:nbGroup){
         entropieTemp <- c(entropieTemp,Tau[i,q]*log(Tau[i,q]));
         tmpNorm <- Tau[i,q]*(log(1/(2*pi^(nbNodes/2)*detVar^(1/2)))-1/2 *((Y[i,] - Mu[,q])%*%VarCovInv)%*%t(t((Y[i,] - Mu[,q]))));
         if (tmpNorm=='Inf'){
           tmpNorm<-mincut;
         }
         Norm <- c(Norm, tmpNorm);
      }
    }
    entropie <- sum(entropieTemp);
    AlphaPrim <- t(Alpha);
    TauPrim <- t(Tau);
    Bern1 <- TauPrim%*%X%*%Tau;
    Bern2 <-  TauPrim%*%(1-X)%*%Tau;
    #print(Norm)
    #print(AlphaPrim)
    #print(TauPrim)
    if (SelfLoop == FALSE){
       Expected <- sum(Tau%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H))-sum((TauPrim%*%Tau)*log(H)) + sum(Norm) - entropie;
       print("alpha");
       print(sum(Tau%*%AlphaPrim));
       print("pi");
       print(sum(Bern1*log(PI) + Bern2*log(H)) -sum((TauPrim%*%Tau)*log(H)));
       print("norm");
       print( sum(Norm));
       print("entropie")
       print(entropie);
    } else {
      Expected <- sum(Tau%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H)) +sum(Norm)- entropie;
    }	
  return(Expected);
}

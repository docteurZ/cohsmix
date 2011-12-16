## __________________________________________________________
##
## EMalgorithmZ
##
## INPUT: Tau		-> Initial classification matrix
##	  Y		-> Symilarity matrix
##	  X		-> Adjacency matrix
##	  NbIteration	-> Number of iteration
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
##
## OUTPUT: MuEstimated	-> Matrix of the estimation of the 
##			   mean 
##	  VarianceEstimated->Estimation of the variance
##	  PIEstimated	-> Estimation of the connectivity 
##			   matrix
##	 AlphaEstimated	-> Estimation of the probability
##			   for a node i to belong to class q
##	 TauEstimated	-> Estimation of the variational paramater
##	 EJ		-> Value of the expected value of J for 
##			   each iteration
##		
## => Plot of EJ curve and network with the estimated classes
##
##
## Algorithm EM when Y and X are independent conditionally to Z 
## __________________________________________________________0

EMalgorithmZ <-function(Tau, Y, X, NbIteration, Plot=TRUE, SelfLoop=FALSE){

############# Notations ##################
  critere <- .Machine$double.xmax;
  nbGroup <- dim(Tau)[2];
  nbNodes <- dim(Tau)[1];

########### Parameters ##
  MuEstimated <- matrix(0,nrow=nbGroup,ncol=nbGroup);
  PIEstimated <- matrix(0,nrow=nbGroup,ncol=nbGroup);
  AlphaEstimated <- vector(length=nbGroup);
  facteur <- matrix(1/nbNodes,nrow=1,ncol=nbNodes);
  TauEstimated <- Tau;
  
  EJ <- 0;
  EJnew <- 0;

  maxcut <- log(.Machine$double.xmax) - log(nbGroup);
  mincut <- log(.Machine$double.xmin);


####### Algorithm ##################
  cat("Iteration: ")
  for(i in 1:NbIteration){
    cat(i," ")
    EJold <- EJ;
    TauPrim <- t(TauEstimated)
    if (SelfLoop == FALSE){
      div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated)-(TauPrim%*%TauEstimated);
    } else {
      div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated);
    }

####### Estimation of the mean #######
   Mutemp <- TauPrim%*%Y%*%TauEstimated;
   MuEstimated <- Mutemp/(div);

   if (length(which(MuEstimated==Inf))>0 || length(which(MuEstimated==-Inf))>0){
     warning("Empty class problem when the number of classes is: ",nbGroup);
     break();
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

####### Estimation of the variance #######
   Ysquare <- Y*Y;
   Musquare <- MuEstimated * MuEstimated;

   Variancetemp <- TauPrim %*% Ysquare %*% TauEstimated -2 * MuEstimated*Mutemp + Musquare*div;
   VarianceEstimated <- sum(Variancetemp)/sum(div);

####### Expected value of J #######
   EJ = c(EJ,ExpectedJ(TauEstimated,Y,X,MuEstimated,PIEstimated,AlphaEstimated,VarianceEstimated,SelfLoop=SelfLoop));	
   EJold=EJ[i];
   if (EJold != 0){
    critere <- abs((EJ[i+1]-EJold)/EJold);
    print(critere);
    if (is.nan(critere)==TRUE){
      break();
    }
    if (critere<=.Machine$double.xmin){
      break();
    }
   } else {
    print("init");
   }


####### Estimation of Tau #######
   Tautemp <- EstimTau(TauEstimated,Y,X,MuEstimated,PIEstimated,AlphaEstimated,VarianceEstimated,SelfLoop=SelfLoop);
   TauEstimated <- Tautemp$TauEstimated;
  }
  cat("\n");

###################### Representation #######################
  if (Plot==TRUE){
    x11();
    par(mfrow = c(1,2));
    cat("Plotting the convergence criterion...","\n");
    plot(EJ[2:length(EJ)]);
    title("Expected value of J: Convergence criterion");
    cat("Plotting the estimated network structure...","\n");
    clusters=getCluster(TauEstimated);
    gplot(X,vertex.col=clusters+1);
    title("Network with estimated classes");
  }

return(list(MuEstimated=MuEstimated,
            VarianceEstimated=VarianceEstimated,
            PIEstimated=PIEstimated,
            AlphaEstimated=AlphaEstimated,
            TauEstimated = TauEstimated,
            EJ = EJ[2:length(EJ)]))
}

## __________________________________________________________
##
## EstimTau
##
## INPUT :	Tau		-> Initial classification matrix
##		Y		-> Symilarity matrix
##		X		-> Adjacency matrix
##		MuX0		-> Matrix of the estimations of the 
##				   means
##		PI		-> Connectivity matrix
##		Alpha		-> Vector  of the probability
##				   for a node i to belong to class q
##		Variance	-> Variance of the random variable Y
##		SelfLoop	-> Equal to FALSE if the self loops 
##				   are not considered
##
## OUTPUT : 	TauEstimated	-> Estimation of the variational paramater
##				   according to the inputs
##
##
## 
## __________________________________________________________

EstimTau <-function(Tau, Y, X, Mu, PI, Alpha, Variance, SelfLoop = FALSE) {
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(Tau)[1];
   TauPrim <- t(Tau);
   TauEstimated <- Tau;
   LogTauEstimated <- matrix(0,nrow=nbNodes,ncol=nbGroup);
   PIPrim <- t(PI);
   H <- 1-PI;
   HPrim <- t(H);

   Ysquare <- Y*Y;
   Musquare <- Mu * Mu;
   MusquarePrim <- t(Musquare);
   MuPrim <- t(Mu);
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
	
   BernMatrix <- ((X%*%Tau)%*%(LogPIPrim)) + (((1-X)%*%Tau)%*%(LogHPrim));
   YsquareTau <- Ysquare%*%Tau;
   MusquareTauPrim <- Musquare%*%TauPrim;
   YTauMuPrim <- (Y%*%Tau)%*%MuPrim;
   TauMusquarePrim <- Tau%*%MusquarePrim;
   TauLogHPrim <- Tau%*%LogHPrim;
	
   rowSumYsquareTau <- rowSums(YsquareTau);
   rowSumMusquareTauPrim <- rowSums(MusquareTauPrim);
   rowSumTau <- rowSums(Tau);
   sumTau <-  sum(Tau);

   for (i in 1:nbNodes){
     for (q in 1:nbGroup){
       Bern <- BernMatrix[i,q];
       Norm <-  (-1/2 * log(2*pi*Variance) * sumTau - 1/(2*Variance) * (rowSumYsquareTau[i] + rowSumMusquareTauPrim[q] - 2 * (YTauMuPrim)[i,q]));
       if (SelfLoop == TRUE){
         LogTauEstimated[i,q] =  LogAlpha[q] + Bern + Norm;	
       } else {
         EgaliteIJ =  (1/2 * log(2*pi*Variance) * rowSumTau[i] + 1/(2*Variance) * (TauMusquarePrim)[i,q] - (TauLogHPrim)[i,q])
         LogTauEstimated[i,q] =  LogAlpha[q] + Bern + (Norm + EgaliteIJ);
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
     Normalize = 1/sum(TauEstimated[i,]);
     TauEstimated[i,] = TauEstimated[i,] * Normalize;
     TauEstimated[i, ][TauEstimated[i, ] < .Machine$double.xmin] <- .Machine$double.xmin;
   }
  return( list(TauEstimated=TauEstimated));
}

## __________________________________________________________
##
## ExpectedJ
##
## INPUT: Tau		-> Initial classification matrix
##	  Y		-> Symilarity matrix
##	  X		-> Adjacency matrix
##	  MuX0		-> Matrix of the estimations of the 
##			   means
##	  PI		-> Connectivity matrix
##	  Alpha		-> Vector  of the probability
##			   for a node i to belong to class q
##	  Variance	-> Variance of the random variable Y
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
## OUTPUT : 	Expected	-> Expected value of J
##
## __________________________________________________________

ExpectedJ <-function(Tau, Y, X, Mu, PI, Alpha, Variance, SelfLoop=FALSE){
   nbGroup = dim(Tau)[2];
   nbNodes = dim(Tau)[1];
   entropieTemp = 0;
   PIPrim = t(PI);
   H = 1-PI;
   HPrim = t(H);

   mincut <- log(.Machine$double.xmin);
   maxcut <- log(.Machine$double.xmax) - log(nbGroup);

################ Entropie #############################
   for (q in 1:nbGroup){
      HPrim[q,]  <- pmin(HPrim[q, ], exp(maxcut));
      HPrim[q,]  <- pmax(HPrim[q, ], exp(mincut));
      PIPrim[q,] <- pmin(PIPrim[q, ], exp(maxcut));
      PIPrim[q,] <- pmax(PIPrim[q, ], exp(mincut));
      Alpha[q]   <- pmin(Alpha[q], exp(maxcut));
      Alpha[q]   <- pmax(Alpha[q], exp(mincut));
      PI[q,]     <- pmin(PI[q,], exp(maxcut));
      PI[q,]     <- pmax(PI[q,], exp(mincut));
      H[q,]      <- pmin(H[q,], exp(maxcut));
      H[q,]      <- pmax(H[q,], exp(mincut));
   }
	
   HPrim[is.nan(HPrim) == TRUE] <- exp(mincut);
   HPrim[HPrim < exp(mincut) ] <- exp(mincut);
   H[is.nan(H) == TRUE] <- exp(mincut);
   H[H < exp(mincut) ] <- exp(mincut);
   for (i in 1:nbNodes){
     Tau[i, ] <- pmin(Tau[i, ], exp(maxcut));
     Tau[i, ] <- pmax(Tau[i, ], exp(mincut));
     for (q in 1:nbGroup){
        entropieTemp = c(entropieTemp,Tau[i,q]*log(Tau[i,q]));
     }
   }
   entropie = sum(entropieTemp);
   AlphaPrim = t(Alpha);
   TauPrim = t(Tau);

   if (SelfLoop == FALSE){
     div = t(t(as.vector(colSums(Tau))))%*%colSums(Tau)-(TauPrim%*%Tau);
   } else {
     div = t(t(as.vector(colSums(Tau))))%*%colSums(Tau);
   }
   
   Mutemp = TauPrim%*%Y%*%Tau;
   Musquare = Mu * Mu;
   Ysquare = Y*Y;
   Bern1 = TauPrim%*%X%*%Tau;
   Bern2 =  TauPrim%*%(1-X)%*%Tau;
   Norm = -1/2 *log(2*pi*Variance)*sum(div) -sum((1/(2*Variance)) * (TauPrim %*% Ysquare %*% Tau -2*Mu*(TauPrim%*%Y%*%Tau) + Musquare*div));
   if (SelfLoop == FALSE){
     Expected = sum(Tau%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H)) - (sum((TauPrim%*%Tau)*log(H)) +Norm - entropie);
   } else {
     Expected = sum(Tau%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H)) + (Norm - entropie);
   }
   return(Expected);
}

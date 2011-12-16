## __________________________________________________________
##
## EMalgorithmXZ
##
## INPUT: Tau		-> Initial classification matrix
##	  Y		-> Symilarity matrix
##	  X		-> Adjacency matrix
##	  NbIteration	-> Number of iteration
##	  SelfLoop	-> Equal to FALSE if the self loops 
##				   are not considered
##
## OUTPUT : 	MuEstimated1	-> Matrix of the estimation of the 
##				   mean for X=1
##		MuEstimated2	-> Matrix of the estimation of the 
##				   mean for X=0
##		VarianceEstimated->Estimation of the variance
##		PIEstimated	-> Estimation of the connectivity 
##				   matrix
##		AlphaEstimated	-> Estimation of the probability
##				   for a node i to belong to class q
##		TauEstimated	-> Estimation of the variational paramater
##		EJ		-> Value of the expected value of J for 
##				   each iteration
##		
## => Plot of EJ curve and network with the estimated classes
## 
## __________________________________________________________

EMalgorithmXZ <- function(Tau, Y, X, NbIteration, Plot=TRUE, SelfLoop=FALSE){

############# Notations ##################
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(Tau)[1];
   AlphaEstimated <- vector(length=nbGroup);
   facteur <- matrix(1/nbNodes,nrow=1,ncol=nbNodes);
   EJ <- 0;
   TauEstimated <- Tau;
   ProdYX <- Y*X;
   ProdY1X <- Y*(1-X);
   EJnew <- 0;
   maxcut <- log(.Machine$double.xmax) - log(nbGroup);
   mincut <- log(.Machine$double.xmin);

   cat("Iteration: ");
   for (i in 1:NbIteration){
     #cat(i," ");
     #print("1)");
     #ptm <- proc.time();
     TauPrim = t(TauEstimated);
     if (SelfLoop == FALSE){
        div = t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated)-(TauPrim%*%TauEstimated);
        divMu1 = TauPrim%*%X%*%TauEstimated;
        divMu2 = TauPrim%*%(1-X)%*%TauEstimated -(TauPrim%*%TauEstimated);
      } else {
        div = t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated);
        divMu1 = TauPrim%*%X%*%TauEstimated
        divMu2 = TauPrim%*%(1-X)%*%TauEstimated
      }

####### Estimation of the mean 1#######
   Mutemp1 = TauPrim%*%ProdYX%*%TauEstimated;
   if (length(which(divMu1==Inf))>0){
     warning("Empty class problem when the number of classes is: ",nbGroup);
     print("moy1div");
     break();
   }
   ## Matrix nbGroup x nbGroup ##
   if (length(divMu1[is.nan(divMu1) == TRUE]) > 0) {
     warning("NAN divMu1: Empty class problem when the number of classes is: ",nbGroup);
     break();
   }
   divMu1[divMu1 < exp(mincut)] <- exp(mincut);
   Mutemp1[Mutemp1 < exp(mincut) ] <- exp(mincut);
   ## MuEstimated1 ##
   MuEstimated1=Mutemp1/(divMu1)
	
   if (length(which(MuEstimated1==Inf))>0){
     warning("Empty class problem when the number of classes is: ",nbGroup)
     break();
   }
####### Estimation of the mean 2#######
   Mutemp2 = TauPrim%*%ProdY1X%*%TauEstimated;
   if (length(divMu2[is.nan(divMu2) == TRUE]) > 0) {
     warning("NAN divMu2: Empty class problem when the number of classes is: ",nbGroup);
     break();
   }
   divMu2[divMu2 < exp(mincut)] <- exp(mincut);
   Mutemp2[Mutemp2 < exp(mincut) ] <- exp(mincut);
   ## MuEstimated2 ##
   MuEstimated2=Mutemp2/(divMu2);
   if (length(which(MuEstimated2==Inf))>0){
     warning("Empty class problem when the number of classes is: ",nbGroup);
     break();
  }

####### Estimation of PI #######
   PItemp = (TauPrim%*%X%*%TauEstimated)
   PIEstimated = PItemp/(div)
   if (length(PIEstimated[is.nan(PIEstimated) == TRUE]) > 0) {
     warning("Empty class problem when the number of classes is: ",nbGroup);
     print("PiNAN");
   }	
   PIEstimated[is.nan(PIEstimated) ==TRUE] <- exp(mincut);
   if (length(PIEstimated[PIEstimated == 'Inf']) > 0) {
     warning("Empty class problem when the number of classes is: ",nbGroup);
     print("PiInf");
   }	
   PIEstimated[PIEstimated == 'Inf'] <- exp(mincut);
   PIEstimated[PIEstimated == 1] <- (1-exp(mincut));
   PIEstimated[PIEstimated < exp(mincut)] <- exp(mincut);

####### Estimation of Alpha #######
   AlphaEstimated <- facteur %*% TauEstimated;

####### Estimation of the variance #######
   Ysquare <- Y*Y;
   Musquare0 <- MuEstimated1 * MuEstimated1;
   Musquare1 <- MuEstimated2 * MuEstimated2;
   Variancetemp = TauPrim %*% Ysquare %*% TauEstimated - 2*MuEstimated2*(TauPrim%*%Y%*%TauEstimated) + 
                 Musquare1*div +(-2*MuEstimated1*Mutemp1 + Musquare0*divMu1)-(-2*MuEstimated2*Mutemp1 + Musquare1*divMu1);

   VarianceEstimated = sum(Variancetemp)/sum(div);
   #print(proc.time() - ptm);
   #print("2)");
   #ptm <- proc.time();
####### Expected value of J #######
   EJnew=ExpectedJXZ(TauEstimated,Y,X,MuEstimated2,MuEstimated1,PIEstimated,AlphaEstimated,VarianceEstimated,SelfLoop=SelfLoop);
   EJ = c(EJ,EJnew);
   EJold=EJ[i];
   if (is.nan(EJold)==TRUE){
      warning("Empty class problem when the number of classes is: ",nbGroup);
      break();
   }
    #print(proc.time() - ptm);
    if (EJold!=0){
      critere <- abs((EJnew-EJold)/EJold);
      print(critere);
      if (is.nan(critere)==TRUE){
        break();
      }
      if (critere<=.Machine$double.xmin){
         break()
      }
    } else {
      print("init");
    }
    #print("3)");
    #ptm <- proc.time();

###### Estimation of Tau #######
   Tautemp <- EstimTauXZ(TauEstimated,Y,X,MuEstimated2,MuEstimated1,PIEstimated,AlphaEstimated,VarianceEstimated,SelfLoop=SelfLoop);
   TauEstimated <- Tautemp$TauEstimated;
   #print(proc.time() - ptm);
   }
   cat("\n")
   
   if (Plot==TRUE){
     #x11()
     par(mfrow = c(1,2));
     cat("Plotting the convergence criterion...","\n");
     plot(EJ[2:NbIteration]);
     title("Expected value of J: Convergence criterion");
     cat("Plotting the estimated network structure...","\n");
     cluster=getCluster(TauEstimated);
     gplot(X,vertex.col=cluster+2);
     title("Network with estimated classes");
   }
   return(list(MuEstimated1=MuEstimated1,
               MuEstimated2=MuEstimated2,
               VarianceEstimated=VarianceEstimated,
               PIEstimated=PIEstimated,
               AlphaEstimated=AlphaEstimated,
               TauEstimated=TauEstimated,
               EJ=EJ[2:length(EJ)]));
}



## __________________________________________________________
##
## EstimTauXZ
##
## INPUT: Tau		-> Initial classification matrix
##	  Y		-> Symilarity matrix
##	  X		-> Adjacency matrix
##	  MuX0		-> Matrix of the estimations of the 
##			   means for X=0
## 	  MuX1		-> Matrix of the estimations of the 
##			   means for X=1
##	  PI		-> Connectivity matrix
##	  Alpha		-> Vector  of the probability
##	         	   for a node i to belong to class q
##	  Variance	-> Variance of the random variable Y
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
##
## OUTPUT : 	TauEstimated	-> Estimation of the variational paramater
##				   according to the inputs
## __________________________________________________________

EstimTauXZ <- function(Tau,Y,X,	MuX0,MuX1, PI, Alpha, Variance, SelfLoop = FALSE, C = FALSE){
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(Tau)[1];
   
   TauPrim <- t(Tau);
   TauEstimated <- Tau;
   LogTauEstimated <- matrix(0,nrow=nbNodes,ncol=nbGroup);

   PIPrim <- t(PI);
   H <- 1-PI;
   HPrim <- t(H);

   Ysquare <- Y*Y;
   Musquare1 <- MuX1 * MuX1;
   Musquare0 <- MuX0 * MuX0;
   Musquare1Prim <- t(Musquare1);
   Mu1Prim <- t(MuX1);
   Musquare0Prim <- t(Musquare0);
   Mu0Prim <- t(MuX0);
   XY2 <- X*Ysquare;
   XY <- X*Y;

   maxcut <- log(.Machine$double.xmax) - log(nbGroup);
   mincut <- log(.Machine$double.xmin);
   for (q in 1:nbGroup){
     HPrim[q,]  <- pmin(HPrim[q, ], exp(maxcut));
     HPrim[q,]  <- pmax(HPrim[q, ], exp(mincut));
     PIPrim[q,] <- pmin(PIPrim[q, ], exp(maxcut));
     PIPrim[q,] <- pmax(PIPrim[q, ], exp(mincut));
     Alpha[q]   <- pmin(Alpha[q], exp(maxcut));
     Alpha[q]   <- pmax(Alpha[q], exp(mincut));
   }	
   ## stupid ? ##
   HPrim[HPrim == "-Inf"] <- exp(mincut);
   HPrim[HPrim == "Inf"] <- exp(mincut);
   HPrim[is.nan(HPrim) == TRUE] <- exp(mincut);
   HPrim[HPrim < exp(mincut) ] <- exp(mincut);
   H[H == "-Inf"] <- exp(mincut);
   H[H == "Inf"] <- exp(mincut);
   H[is.nan(H) == TRUE] <- exp(mincut);
   H[H < exp(mincut) ] <- exp(mincut);
   ##         ##	

   LogAlpha  <- log(Alpha);
   LogPIPrim <- log(PIPrim);
   LogHPrim  <-  log(HPrim);

   BernMatrix = ((X%*%Tau)%*%(LogPIPrim)) + (((1-X)%*%Tau)%*%(LogHPrim));
   sumTau <- sum(Tau);
   YsquareTau <- Ysquare%*%Tau;
   Musquare0TauPrim <- Musquare0%*%TauPrim;
   YTauMu0Prim <- (Y%*%Tau)%*%Mu0Prim;
   XY2Tau <- XY2%*%Tau;
   XYTauMu1Prim <- XY%*%(Tau%*%Mu1Prim);
   XTauMusquare1Prim <- X%*%(Tau%*%Musquare1Prim);
   XY2Tau <- XY2%*%Tau;
   XYTauMu0Prim <- XY%*%(Tau%*%Mu0Prim);
   XTauMusquare0Prim <- X%*%(Tau%*%Musquare0Prim);
   TauMusquare0Prim <- Tau%*%Musquare0Prim;
   TauLogHPrim <- Tau%*%LogHPrim;
   rowSumYsquareTau <- rowSums(YsquareTau);
   rowSumMusquare0TauPrim <- rowSums(Musquare0TauPrim);
   rowSumXY2Tau <- rowSums(XY2Tau);
   rowSumTau <- rowSums(Tau);
   for (i in 1:nbNodes){
     for (q in 1:nbGroup){
        Bern <- BernMatrix[i,q];
	Norm <-  (-1/2 * log(2*pi*Variance) * sumTau - 1/(2*Variance) * (rowSumYsquareTau[i] + rowSumMusquare0TauPrim[q] -2 * (YTauMu0Prim)[i,q])) +
 1/(2*Variance) * (-(rowSumXY2Tau[i] -2*(XYTauMu1Prim)[i,q] + (XTauMusquare1Prim)[i,q]) + (rowSumXY2Tau[i] -2*(XYTauMu0Prim)[i,q] + (XTauMusquare0Prim)[i,q]));		
        if (SelfLoop == FALSE){
           EgaliteIJ <-  (1/2 * log(2*pi*Variance) * rowSumTau[i] + 1/(2*Variance) * (TauMusquare0Prim)[i,q] - (TauLogHPrim)[i,q])
	   LogTauEstimated[i,q] <-  LogAlpha[q] + Bern + Norm + EgaliteIJ;
        } else {
          LogTauEstimated[i,q] <-  LogAlpha[q] + Bern + Norm;
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
## ExpectedJXZ
##
## INPUT :	Tau		-> Initial classification matrix
##		Y		-> Symilarity matrix
##		X		-> Adjacency matrix
##		MuX0		-> Matrix of the estimations of the 
##				   means for X=0
## 		MuX1		-> Matrix of the estimations of the 
##				   means for X=1
##		PI		-> Connectivity matrix
##		Alpha		-> Vector  of the probability
##				   for a node i to belong to class q
##		Variance	-> Variance of the random variable Y
##		SelfLoop	-> Equal to FALSE if the self loops 
##				   are not considered
##
## OUTPUT : 	Expected	-> Expected value of J
##
## __________________________________________________________

ExpectedJXZ <-function(Tau, Y, X, MuX0, MuX1, PI, Alpha, Variance, SelfLoop=FALSE){
   nbGroup <- dim(Tau)[2];
   nbNodes <- dim(Tau)[1];
   entropieTemp <- 0;
   PIPrim <- t(PI);
   H <- 1-PI;
   HPrim <- t(H);
   ProdYX <- Y*X;
   ProdY1X <- Y*(1-X);
   maxcut <- log(.Machine$double.xmax) - log(nbGroup);
   mincut <- log(.Machine$double.xmin);

################ Entropie #############################
   for (q in 1:nbGroup){
      HPrim[q,]  <- pmin(HPrim[q, ], exp(maxcut));
      HPrim[q,]  <- pmax(HPrim[q, ], exp(mincut));
      PIPrim[q,] <- pmin(PIPrim[q, ], exp(maxcut));
      PIPrim[q,] <-  pmax(PIPrim[q, ], exp(mincut));
      Alpha[q]   <- pmin(Alpha[q], exp(maxcut));
      Alpha[q]   <- pmax(Alpha[q], exp(mincut));
      PI[q,]     <- pmin(PI[q,], exp(maxcut));
      PI[q,]     <- pmax(PI[q,], exp(mincut));
      H[q,]      <- pmin(H[q,], exp(maxcut));
      H[q,]      <- pmax(H[q,], exp(mincut));
   }
   # usefull ##
   HPrim[is.nan(HPrim) == TRUE] <- exp(mincut);
   HPrim[HPrim < exp(mincut) ] <- exp(mincut);
   H[is.nan(H) == TRUE] <- exp(mincut);
   H[H < exp(mincut) ] <- exp(mincut);
   ## ##########	
   for (i in 1:nbNodes){
     Tau[i, ] <- pmin(Tau[i, ], exp(maxcut))
     Tau[i, ] <- pmax(Tau[i, ], exp(mincut))
     for (q in 1:nbGroup){
        entropieTemp = c(entropieTemp,Tau[i,q]*log(Tau[i,q]));
     }
   }
   entropie <- sum(entropieTemp);
   AlphaPrim <- t(Alpha);
   TauPrim <- t(Tau);
   Bern1 <- TauPrim%*%X%*%Tau;
   Bern2 <-  TauPrim%*%(1-X)%*%Tau;
   Ysquare <- Y*Y;
   Musquare0 <- MuX0 * MuX0;
   Musquare1 <- MuX1 * MuX1;	
   if (SelfLoop == FALSE){
     div <- t(t(as.vector(colSums(Tau))))%*%colSums(Tau)-(TauPrim%*%Tau);
   } else {
     div <- t(t(as.vector(colSums(Tau))))%*%colSums(Tau);
   }

   divMu1 <- TauPrim%*%X%*%Tau
   Mutemp1 <- TauPrim%*%ProdYX%*%Tau
   if (SelfLoop == FALSE){
      Expected <- sum(Tau%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H))-sum((TauPrim%*%Tau)*log(H))- 1/2 *log(2*pi*Variance)* sum(div) + sum((1/(2*Variance)) * (-(TauPrim %*% Ysquare %*% Tau -2*MuX0*(TauPrim%*%Y%*%Tau) + Musquare0*div) +(-2*MuX0*Mutemp1 + Musquare0*divMu1)-(-2*MuX1*Mutemp1 + Musquare1*divMu1)))- entropie;
    } else {
     Expected <- sum(Tau%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H))- 1/2 *log(2*pi*Variance)* sum(div) +sum((1/(2*Variance)) * (-(TauPrim %*% Ysquare %*% Tau -2*MuX0*(TauPrim%*%Y%*%Tau) + Musquare1*div) +(-2*MuX0*Mutemp1 + Musquare0*divMu1)-(-2*MuX0*Mutemp1 + Musquare1*divMu1)))- entropie;
   }
  return(Expected);
}

## __________________________________________________________
##
## ICLXZ
##
## INPUT: Y		-> Symilarity matrix
##	  X		-> Adjacency matrix
##	  Qmax		-> Maximal number of classes
##	  NbIteration	-> Number of iteration
##	  Plot		-> Equal to TRUE if we want a 
##	         	   plot of ICL
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
## OUTPUT: 	Value of ICL criterion for each number of classes
##		
## => Plot of ICL curve 
##
## Calculation of the ICL when value of Y depends on the edges
## __________________________________________________________

ICLXZ <-function(Y, X, Qmin, Qmax, NbIteration, loop = 10, Plot=TRUE, SelfLoop=FALSE){
  nbNodes = dim(Y)[1];
  result=0;
  resultEM=vector(mode="list",length=((Qmax-Qmin)+1));
  TempEst=0;
  maxRatioAlpha <- 0;
  maxRatioEdge <- 0;
  maxRatioCov <- 0;
  for (q in Qmin : Qmax){
    maxcut <- log(.Machine$double.xmax) - log(q);
    mincut <- log(.Machine$double.xmin);
    ICL<- -Inf;
    ratioAlpha <- 0;
    ratioEdge <- 0;
    ratioCov <- 0;
    cat("Number of classes :",q,"\n")
    for (j in 1 : loop){
      cat("Number of loops :",j,"\n")
      Tau <- t(rmultinom(size=1,rep(1/q,q),n=nbNodes));
      CalculEst <- EMalgorithmXZ(Tau,Y,X,NbIteration,Plot=FALSE,SelfLoop=SelfLoop)
      TauEstimated <- CalculEst$TauEstimated
      TauEstimated[,q] <- pmin(TauEstimated[,q], exp(maxcut))
      TauEstimated[,q] <- pmax(TauEstimated[,q], exp(mincut))
      TauPrim <- t(TauEstimated)
      if (SelfLoop == FALSE){
        div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated)-(TauPrim%*%TauEstimated);
        divMu1 <- TauPrim%*%X%*%TauEstimated
        divMu2 <- TauPrim%*%(1-X)%*%TauEstimated -(TauPrim%*%TauEstimated)
      } else {
        div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated);
        divMu1 <- TauPrim%*%X%*%TauEstimated;
        divMu2 <- TauPrim%*%(1-X)%*%TauEstimated;
      }
      ProdYX=Y*X;
      ProdY1X=Y*(1-X);
      Alpha <- CalculEst$AlphaEstimated;
      Variance <- CalculEst$VarianceEstimated;
      MuX0 <- CalculEst$MuEstimated2;
      MuX1 <- CalculEst$MuEstimated1;
      PI <- CalculEst$PIEstimated;
      H <- 1-PI;
      PI[,q] <- pmin(PI[,q], exp(maxcut));
      H[,q] <- pmin(H[,q], exp(maxcut));
      PI[,q] <- pmax(PI[,q], exp(mincut));
      H[,q] <- pmax(H[,q], exp(mincut));
      AlphaPrim <- t(Alpha);
      TauPrim <- t(TauEstimated);
      Bern1 <- TauPrim%*%X%*%TauEstimated ;
      Bern2 <-  TauPrim%*%(1-X)%*%TauEstimated ;
      Mutemp1 <- TauPrim%*%ProdYX%*%TauEstimated;
      Ysquare <- Y*Y
      Musquare1 <- MuX0 * MuX0
      Musquare2 <- MuX1 * MuX1
      alphaPart <- sum(TauEstimated%*%AlphaPrim);
      if (SelfLoop == FALSE){
        edgePart <- sum(Bern1*log(PI) + Bern2*log(H))-sum((TauPrim%*%TauEstimated)*log(H))- 1/2 *log(2*pi*Variance)* sum(div);
        covPart <-  sum((1/(2*Variance)) * (-(TauPrim %*% Ysquare %*% TauEstimated -2*MuX0*(TauPrim%*%Y%*%TauEstimated) + Musquare1*div) 
                    +(-2*MuX0*Mutemp1 + Musquare1*divMu1)-(-2*MuX1*Mutemp1 + Musquare2*divMu1)));           
        LogLikelihood <- alphaPart + edgePart + covPart;
     } else {
        edgePart <- sum(Bern1*log(PI) + Bern2*log(H))- 1/2 *log(2*pi*Variance)* sum(div);
        covPart  <- sum((1/(2*Variance)) * (-(TauPrim %*% Ysquare %*% TauEstimated -2*MuX0*(TauPrim%*%Y%*%TauEstimated) + Musquare1*div) 
                    +(-2*MuX0*Mutemp1 + Musquare1*divMu1)-(-2*MuX1*Mutemp1 + Musquare2*divMu1)));                     
        LogLikelihood <- alphaPart + edgePart + covPart;
     }

     TempICL <- LogLikelihood - (q-1)*log(nbNodes) -3*q*q*log(nbNodes*(nbNodes-1)/2) -log(nbNodes*(nbNodes-1)/2);
     if (is.nan(TempICL)==TRUE){
       TempICL<-(-exp(maxcut));
     }
     ICL <- c(ICL,TempICL);
     ratioAlpha <- c(ratioAlpha, alphaPart/LogLikelihood);
     ratioEdge  <- c(ratioEdge, edgePart/LogLikelihood);
     ratioCov   <- c(ratioCov, covPart/LogLikelihood);
     if (TempICL >= max(ICL[2:length(ICL)])){
       TempEst <- CalculEst;
     }
   }
   imaxICL <- which.max(ICL);
   print(imaxICL);
   maxICL <- ICL[imaxICL] # get Max
   result <- c(result,maxICL);
   resultEM[[q]] <- TempEst;
   maxRatioAlpha <- c(maxRatioAlpha,ratioAlpha[imaxICL]);
   maxRatioEdge <- c(maxRatioEdge,ratioEdge[imaxICL]);
   maxRatioCov <- c(maxRatioCov,ratioCov[imaxICL]);
  }
  if (Plot==TRUE){
  ##################### Representation #######################
    x11()
    result<-result[2:length(result)]
    maxi<-which(max(result)==result)
    plot(result,type="b",xlab="Number of classes",ylab="ICL value")
    points(maxi,result[maxi],col="red")
    title("ICL curve")
  ##############################################################
  }
  return(list(ICL=result[2:length(result)],
              EMestimation=resultEM,
              ratioAlpha=maxRatioAlpha,
              ratioEdge=maxRatioEdge,
              ratioCov=maxRatioCov));
}




## __________________________________________________________
##
## ICL
##
## INPUT :	Y		-> Symilarity matrix
##		X		-> Adjacency matrix
##		Qmax		-> Maximal number of classes
##		NbIteration	-> Number of iteration
##		Plot		-> Equal to TRUE if we want a 
##				   plot of ICL
##		SelfLoop	-> Equal to FALSE if the self loops 
##				   are not considered
##
## OUTPUT : 	Value of ICL criterion for each number of classes
##		
## => Plot of ICL curve
##
## Calculation of the ICL when value of Y does not 
## depends on the edges
## __________________________________________________________

ICLZ <-function(Y, X, Qmin, Qmax, NbIteration, loop = 10, Plot=TRUE, SelfLoop=FALSE){
   nbNodes <- dim(Y)[1];
   result<-0;
   resultEM<-vector(mode="list",length=((Qmax-Qmin)+1))
   TempEst<-0
   for (q in Qmin : Qmax){
     maxcut <- log(.Machine$double.xmax) - log(q);
     mincut <- log(.Machine$double.xmin);
     ICL <- 0;
     cat("Number of classes :",q,"\n");
     for (j in 1 : loop){
       cat("Number of loops :",j,"\n");
       Tau = t(rmultinom(size=1,rep(1/q,q),n=nbNodes));
       CalculEst = EMalgorithmZ(Tau,Y,X,NbIteration,Plot=FALSE,SelfLoop=SelfLoop);
       TauEstimated = CalculEst$TauEstimated;
       TauEstimated[,q] <- pmin(TauEstimated[,q], exp(maxcut));
       TauEstimated[,q] <- pmax(TauEstimated[,q], exp(mincut));

       TauPrim <- t(TauEstimated);
       if (SelfLoop == FALSE){
         div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated)-(TauPrim%*%TauEstimated)
       } else {
         div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated);
       }

       ProdYX<-Y*X;
       ProdY1X<-Y*(1-X);
       Alpha <- CalculEst$AlphaEstimated;
       Variance <- CalculEst$VarianceEstimated;
       MuX0 <- CalculEst$MuEstimated;
       PI <- CalculEst$PIEstimated;
       H <- 1-PI;
       PI[,q] <- pmin(PI[,q], exp(maxcut));
       H[,q] <- pmin(H[,q], exp(maxcut));
       PI[,q] <- pmax(PI[,q], exp(mincut));
       H[,q] <- pmax(H[,q], exp(mincut));

       AlphaPrim <- t(Alpha);
       TauPrim <- t(TauEstimated);
       Bern1 <- TauPrim%*%X%*%TauEstimated ;
       Bern2 <-  TauPrim%*%(1-X)%*%TauEstimated;
       Mutemp1 <- TauPrim%*%ProdYX%*%TauEstimated;
       Ysquare <- Y*Y;
       Musquare1 <- MuX0 * MuX0;
       if (SelfLoop == FALSE){
         LogLikelihood= sum(TauEstimated%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H))-sum((TauPrim%*%TauEstimated)*log(H)) 
          - 1/2 *log(2*pi*Variance)* sum(div) -sum((1/(2*Variance)) * (TauPrim %*% Ysquare %*% TauEstimated -2*MuX0*(TauPrim%*%Y%*%TauEstimated) + Musquare1*div));
       } else {
         LogLikelihood= sum(TauEstimated%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H)) 
          + 1/2 *log(2*pi*Variance)* sum(div) +sum((1/(2*Variance)) * (TauPrim %*% Ysquare %*% TauEstimated -2*MuX0*(TauPrim%*%Y%*%TauEstimated) + Musquare1*div));
       }
       TempICL <- LogLikelihood - (q-1)*log(nbNodes) -log(nbNodes*(nbNodes-1)/2)-2*q*q*log(nbNodes*(nbNodes-1)/2);
       if (is.nan(TempICL)==TRUE){
         TempICL<-(-exp(maxcut));
       }
       ICL  <- c(ICL,TempICL);
       if (TempICL>=max(ICL[2:length(ICL)])){
          TempEst<-CalculEst;
       }
     }
     maxICL<-max(ICL[2:length(ICL)])
     result <- c(result,maxICL)
     resultEM[[q]] <- TempEst;
    }

    if (Plot==TRUE){
    ##################### Representation #######################
      x11()
      result <- result[2:length(result)];
      maxi <- which(max(result)==result);
      plot(result,type="b",xlab="Number of classes",ylab="ICL value");
      points(maxi,result[maxi],col="red");
      title("ICL curve");
      ##############################################################
    }
  return(list(iclvalues=result[2:length(result)],EMestimation=resultEM))
}

## _________________________________________________________
##
## ICLNodes
##
## INPUT: Y		-> Symilarity matrix
##	  X		-> Adjacency matrix
##	  Qmax		-> Maximal number of classes
##	  NbIteration	-> Number of iteration
##	  Plot		-> Equal to TRUE if we want a 
##	 		   plot of ICL
##	  SelfLoop	-> Equal to FALSE if the self loops 
##			   are not considered
##
## OUTPUT : 	Value of ICL criterion for each number of classes
##		
## => Plot of ICL curve
##
## Calculation of the ICL when value of Y does not 
## depends on the edges
## __________________________________________________________

ICLNodes <-function(Y, X, Qmin, Qmax, NbIteration, loop = 10, Plot=TRUE, SelfLoop=FALSE){
   nbNodes <- dim(Y)[1];
   nbCov <- dim(Y)[2];
   result <- 0;
   TempEst <- 0;
   resultEM <- vector(mode="list",length=((Qmax-Qmin)+1));
   
   for (q in Qmin : Qmax) {
      maxcut <- log(.Machine$double.xmax) - log(q);
      mincut <- log(.Machine$double.xmin);
      cat("Number of classes :",q,"\n");
      ICL <- 0;
      for (j in 1 : loop){
        Norm <- 0;
        cat("Number of loops :",j,"\n");
        ## init with spectral clustering ##
        S <- buildSimilarity(Y);
        L1 <- LaplaceMatrix(S);
        L2 <- LaplaceMatrix(X);
        Lsum <- L1 + L2;
        sc <- spectralClustering(Lsum, q);
        Tau <- class.ind(sc$cluster); 
        #Tau <- t(rmultinom(size=1,rep(1/q,q),n=nbNodes));	
        
        CalculEst <- EMalgorithmNodes(Tau,Y,X,NbIteration,Plot=FALSE,SelfLoop=FALSE);
        TauEstimated <- CalculEst$TauEstimated;
        TauEstimated[,q] <- pmin(TauEstimated[,q], exp(maxcut));
        TauEstimated[,q] <- pmax(TauEstimated[,q], exp(mincut));
	
        TauPrim <- t(TauEstimated);
        div <- t(t(as.vector(colSums(TauEstimated))))%*%colSums(TauEstimated)-(TauPrim%*%TauEstimated);
        
        Alpha <- CalculEst$AlphaEstimated;
        VarCov <- CalculEst$VarCovEstimated;
        VarCovInv <- ginv(VarCov);
        detVar <- det(VarCov);
        PI <- CalculEst$PIEstimated;
        Mu <- CalculEst$Mu;
        H <- 1-PI;
        PI[,q] <- pmin(PI[,q], exp(maxcut));
        H[,q] <- pmin(H[,q], exp(maxcut));
        PI[,q] <- pmax(PI[,q], exp(mincut));
        H[,q] <- pmax(H[,q], exp(mincut));
	
        AlphaPrim <- t(Alpha);
        TauPrim <- t(TauEstimated);

        Bern1 <- TauPrim%*%X%*%TauEstimated;
        Bern2 <-  TauPrim%*%(1-X)%*%TauEstimated;
			
        for (i in 1:nbNodes){
          for (h in 1:q){
             tmpNorm <- TauEstimated[i,h]*(log(1/(2*pi^(nbNodes/2)*detVar^(1/2)))-1/2 *((Y[i,] - Mu[,h])%*%VarCovInv)%*%t(t((Y[i,] - Mu[,h]))));
             if (tmpNorm =='Inf'){
                tmpNorm <- mincut;
             }
             Norm <- c(Norm, tmpNorm);
             #print(Norm);
	   }
         }
			
         LogLikelihood <- sum(TauEstimated%*%AlphaPrim) + sum(Bern1*log(PI) + Bern2*log(H))-sum((TauPrim%*%TauEstimated)*log(H)) + sum(Norm);

         TempICL <- LogLikelihood - (q-1)*log(nbNodes)     #alpha
                        - nbCov*(nbCov-1)*log(nbNodes*(nbNodes-1)/2)  # zeta
                         -q*(q-1)*log(nbNodes*(nbNodes-1)/2)/2 #pi
                         - q*nbCov*log(nbNodes*(nbNodes-1)/2); #sigma
			
         if (is.nan(TempICL)==TRUE){TempICL<-(-exp(maxcut))}
           ICL  <- c(ICL,TempICL)
           if (TempICL>=max(ICL[2:length(ICL)])){
              TempEst <- CalculEst;
           }
          }
          maxICL <- max(ICL[2:length(ICL)]);
          result <- c(result,maxICL);
          resultEM[[q]] <- TempEst;
       }

     if (Plot==TRUE){
       ###################### Representation #######################
        x11();
        plot(result[2:length(result)]);
        title("ICL curve");
     }
    return(list(result[2:length(result)],EMestimation=resultEM));
}




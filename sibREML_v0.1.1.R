## Data format
## ID1 ID2 PHENO1 PHENO2 IBD_1 IBD2_2 ... IBD_k
## No missing data...
## No covariates other than intercept...
sibREML <- function(data,niter=50,tol=1e-6,verbose=TRUE,addIntercept=TRUE,niterEM=1,startValues=NULL){
  Np <- nrow(data)
  Y  <- as.matrix(data[,3:4]) # 3rd column
  if(addIntercept){
    GRM   <- cbind(C=1,data[,-(1:4)])
  }else{
    GRM   <- cbind(data[,-(1:4)])
  }
  nvars <- ncol(GRM)
  A     <- lapply(1:nvars, function(k){
    lapply(1:Np, function(i) matrix(c(1,GRM[i,k],GRM[i,k],1),2,2))
  })
  if(verbose){
    cat(paste0("Detected ",Np," pairs and ",nvars," variance components.\n"))
  }
  
  ## Generate summary statistics
  rp <- cor(Y[,1],Y[,2]); v1 <- var(Y[,1]); v2 <- var(Y[,2]); m1 <- mean(Y[,1]); m2 <- mean(Y[,2])
  cat(paste0("Pair-correlation ",round(rp,2)," | E[Y1] = ",round(m1,2), " - E[Y2] = ",round(m2,2),
             " | var[Y1] = ",round(v1,2), " - var[Y2] = ",round(v2,2),".\n"))
  
  ## Initialize
  if(is.null(startValues)){
    Vp <- 0.5 * (v1 + v2)
    VC <- rep(Vp/(nvars+1),nvars+1)
    Vg <- VC[1:nvars]
    Ve <- VC[1+nvars]
  }else{
    eVC <- nvars + 1
    nVC <- length(startValues)
    if(nVC != eVC){
      cat(paste("[Error] The number of input parameters ([`startValues=` option]) is ",
                nVC," and does not match the number of variance components expected (i.e.,",
                eVC,").\n"))
      cat("[Error] Returning NULL")
      return(NULL)
    }else{
      VC <- startValues
      Vg <- VC[1:nvars]
      Ve <- VC[1+nvars]
      Vp <- sum(Vg) + Ve
      pf <- paste0("Starting Values for variance components: Vg(1) = ",round(Vg[1],5))
      k  <- 2
      while(k<=nvars){
        pf <- c(pf,paste0(" - Vg(",k,") = ",round(Vg[k],5)))
        k  <- k + 1
      }
      pf <- c(pf,paste0(" - Ve = ",round(Ve,5)))
      if(verbose){
        cat(paste0(paste(pf,collapse = ""),".\n"))
      }
    }
  }

  Vm <- do.call("cbind",lapply(1:nvars,function(k) Vg[k] * GRM[,k]))
  J2 <- matrix(1,2,2)

  AI <- matrix(0,nrow=1+nvars,ncol=1+nvars)
  dL <- rep(0,1+nvars)
  
  eps <- 1
  
  ## History
  logLikHistory <- NULL
  ParamHistory  <- NULL
  ## Running AI-REML algorithm
  if(verbose){
    cat("Running AI-REML algorithm...\n")
  }
  for(it in 1:niter){
    SVm <- rowSums(Vm)
    V      <- lapply(1:Np,function(i){
      matrix(c(Vp,SVm[i],SVm[i],Vp),2,2) 
    })
    V_     <- lapply(V,function(v) matrix(c(v[2,2],-v[1,2],-v[2,1],v[1,1]),2,2)/(v[1,1]*v[2,2]-v[1,2]*v[2,1])) 
    SV_    <- sapply(V_, function(v_) sum(v_))
    xTV_x  <- sum(SV_)
    P      <- lapply(1:Np,function(i) V_[[i]] - 0.25 * SV_[i] * SV_[i] * J2  / xTV_x)
    #P      <- lapply(1:Np,function(i) V_[[i]] - (1/xTV_x) * crossprod(V_[[i]],crossprod(J2,V_[[i]])) )
    PA     <- lapply(1:nvars, function(k){
      lapply(1:Np, function(i) crossprod(P[[i]],A[[k]][[i]]) )
    })
    Py     <- lapply(1:Np, function(i){
      crossprod(P[[i]],Y[i,])
    })
    PAPy   <- lapply(1:nvars, function(k){
      lapply(1:Np, function(i){
       #PA[[k]][[i]]%*%Py[[i]] 
       crossprod(PA[[k]][[i]],Py[[i]])
    })})
    PPy    <- lapply(1:Np, function(i){
      crossprod(P[[i]],Py[[i]])
    })
    APy    <- lapply(1:nvars, function(k){
      lapply(1:Np, function(i){
        crossprod(A[[k]][[i]],Py[[i]])
      })
    })
    yPAPy  <- sapply(1:nvars, function(k){
      sum( sapply(1:Np, function(i){
        crossprod(APy[[k]][[i]],Py[[i]])
      }) )
    })
    yPPy   <- sum( sapply(1:Np, function(i) crossprod(Py[[i]]) ) )
    trPA   <- sapply(1:nvars, function(k){
      sum( sapply(1:Np, function(i) sum(diag(PA[[k]][[i]])) ) )
    })
    trP    <- sum( sapply(1:Np, function(i) sum(diag(P[[i]])) ) )
    
    ## Fixed effect: mu = (X'V-1X)-1 (X'V-1y)
    V_y     <- sapply(1:Np, function(i){
      sum( crossprod(V_[[i]],Y[i,]) )
    })
    xTV_y <- sum(V_y)
    mu    <- xTV_y / xTV_x
    
    ## Calculate log-likelihood
    logDet_V <- sum( sapply(1:Np, function(i) log(det(V[[i]]))  ) )
    yPy      <- sum( sapply(1:Np, function(i) crossprod(Y[i,],Py[[i]]) ) )
    logLik   <- -0.5 * (logDet_V + log(xTV_x) + yPy )
    
    ## Record history
    logLikHistory  <- c(logLikHistory,logLik)
    ParamHistory   <- rbind(ParamHistory,VC)
    if(it>1){
      eps <- abs( logLik - prev_logLik )
    }
    prev_logLik <- logLik
    pf <- paste0("[it=",it,"] - mu = ",round(mu,5))
    for(k in 1:nvars){
      pf <- c(pf,paste0(" - Vg(",k,") = ",round(Vg[k],5)))
    }
    pf <- c(pf,paste0(" - Ve = ",round(Ve,5),
                      " - logLik = ",round(logLik,5),
                      " - eps = ",format(eps,scientific=T,digits=5)))
    if(verbose){
      cat(paste0(paste(pf,collapse = ""),".\n"))
    }
    ## Calculate AI matrix
    for(k in 1:nvars){
      for(l in k:nvars){
        AI[k,l] <- AI[l,k] <- sum( sapply(1:Np,function(i) crossprod(APy[[k]][[i]],PAPy[[l]][[i]]) ) )
      }
      AI[k,nvars+1] <- AI[nvars+1,k] <- sum( sapply(1:Np,function(i) crossprod(APy[[k]][[i]],PPy[[i]]) ) )
    }
    AI[nvars+1,nvars+1] <- sum( sapply(1:Np,function(i) crossprod(Py[[i]],PPy[[i]]) ) )
    AI <- 0.5 * AI
    
    ## Calculate gradient
    for(k in 1:nvars){
      dL[k] <- -0.5*(trPA[k] - yPAPy[k])
    }
    dL[1+nvars] <- -0.5*(trP - yPPy)

    ## Update parameters
    if(it<=niterEM){ # EM iteration
      dx <- VC*VC*(2*dL)/Np
      VC <- VC + dx
      Vg <- VC[1:nvars]
      Ve <- VC[1+nvars]
    }else{
      dx <- solve(AI,dL)
      VC <- VC + dx
      Vg <- VC[1:nvars]
      Ve <- VC[1+nvars]
    }
    Vp <- sum(Vg) + Ve
    for(k in 1:nvars){
      Vm[,k] <- Vg[k] * GRM[,k]
    }
    
    if(eps<tol){
      break
    }
  }

  ## Prepare output  
  COV_VC <- solve(AI)
  SE_VC  <- sqrt(diag(COV_VC))
  Vp     <- sum(Vg) + Ve
  
  var_Vg          <- diag(COV_VC)[1:nvars]
  var_Vp          <- sum(COV_VC)
  var_Vg_total    <- sum(COV_VC[1:nvars,1:nvars])
  
  if(nvars>1){
    cov_Vp_Vg       <- rowSums(COV_VC[1:nvars,])
  }else{
    cov_Vp_Vg       <- sum(COV_VC[1:nvars,])
  }
  
  cov_Vp_Vg_total <- sum(COV_VC[1:nvars,])
  
  ## Heritability and SE for each component
  h2_g     <- Vg / Vp
  var_h2_g <- (h2_g^2)*( var_Vg/(Vg^2) - 2*cov_Vp_Vg/(Vg*Vp) + var_Vp/(Vp^2) )
  
  ## Heritability and SE for total genetic variance (including intercept)
  Vg_total <- sum(Vg)
  h2_total <- Vg_total / Vp
  var_h2_total <- (h2_total^2)*( var_Vg_total/(Vg_total^2) - 2*cov_Vp_Vg_total/(Vg_total*Vp) + var_Vp/(Vp^2) )
  
  ## If intercept included
  h2_noIntercept <- NULL
  if(addIntercept){
    Vg_noInt        <- sum(VC[2:nvars])
    var_Vg_noInt    <- sum(COV_VC[2:nvars,2:nvars])
    cov_Vp_Vg_noInt <- sum(COV_VC[2:nvars,])
    
    h2_noInt       <- Vg_noInt / Vp
    var_h2_noInt   <- (h2_noInt^2)*( var_Vg_noInt/(Vg_noInt^2) - 2*cov_Vp_Vg_noInt/(Vg_noInt*Vp) + var_Vp/(Vp^2) )
    se_h2_noInt    <- sqrt(var_h2_noInt)
    h2_noIntercept <- c(Estimate=h2_noInt,SE=se_h2_noInt)
  }
  colnames(ParamHistory) <- c(paste0("Varcomp",1:nvars),"Ve")
  ParamHistory <- cbind(Iter=1:nrow(ParamHistory),ParamHistory)
  
  logLikHistory <- cbind(Iter=1:length(logLikHistory),LogLik=logLikHistory)
  results  <- list(h2_g=rbind(Estimate=h2_g,SE=sqrt(var_h2_g)),
                   h2_total=c(Estimate=h2_total,SE=sqrt(var_h2_total)),
                   h2_noInt=h2_noIntercept,
                   FixedEffect=mu,
                   Vg=rbind(Estimate=Vg,SE=SE_VC[1:nvars]),
                   Ve=c(Estimate=Ve,SE=SE_VC[1+nvars]),
                   Vp=c(Estimate=Vp,SE=sqrt(var_Vp)),
                   VC=rbind(Estimate=VC,SE=SE_VC),
                   COV_VC=COV_VC,
                   logLik=logLik,
                   logLik_history=logLikHistory,
                   Param_history=ParamHistory)
  return(results)
}

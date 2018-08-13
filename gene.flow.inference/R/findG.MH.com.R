findG.MH.com <- function(H,G_adj,iter=10000,fixed_start=FALSE,g_init=1,q_init=1,
                      fixed_sd=FALSE,fixed_sd_start=FALSE,sdG_init=1/200,sdq_init=1/200,
                      sig2ep=(mean(H,na.rm=TRUE)/50)^2){
  # commute time version
  # uses commute time as if it were coalescence time
  # make sure matrix is square and find dimentions
  if(NROW(H) == NCOL(H)) n <- NROW(H)
  else stop("Matrix is not square. Please enter a symmetric matrix")
  
  # make sure G and H have the same dimensions
  if((NROW(H) != NROW(G_adj)) || (NCOL(H) != NCOL(G_adj))) stop("H and G_adj must have the same dimensions")
  
  #make sure matrix is symmetric up to computational error
  if(max(H-Matrix::t(H),na.rm=TRUE)/max(H,na.rm=TRUE) > 1e-10) stop("Matrix is not symmetric. Please enter a symmetric matrix")
  #make H exactly symmetric
  H <- (H+Matrix::t(H))/2
  #change H to dsCMatrix if it is not already
  #H <- as(H,"dsCMatrix") #don't need to do this
  #extract values into vector
  #h <- H@x
  h <- as.vector(H[upper.tri(H,diag=TRUE)])
  
  #return(h)
  #make sure sig2ep is the right length (either a scalar or same length as h)
  if((length(sig2ep) != 1) && (length(sig2ep) != length(h))){
    print(length(sig2ep))
    print(length(h))
    stop("sig2ep must be either a scalar or the same length as h")
  }
  #make sure G_adj is in sparse matrix format and add diagonal for structure matrix
  #G_adj <- as(G_adj,"dgCMatrix")
  #G_struct <- G_adj
  #diag(G_struct) <- 1
  #future: add check to make sure all non-zero entries are 1
  
  #find number of elements in g and q
  ng <- length(G_adj@x)
  nq <- n
  
  #set up priors for G and qma
  meanG <- 1;
  meanq <- 1;
  #convert to distribution parameters
  lamG <- 1/meanG
  lamq <- 1/meanq
  #decide variance for error on values of H
  #sig2ep <- (mean(H,na.rm=TRUE)/50)^2
  #set up parameters for the proposal distribution
  if(fixed_sd_start == FALSE){
    sdG <- rep(meanG/200,iter)
    sdq <- rep(meanq/200,iter)
  }
  else{
    sdG <- rep(sdG_init,iter)
    sdq <- rep(sdq_init,iter)
  }
  
  #initialize arrays for g and q
  g <- array(dim=c(iter,ng))
  q <- array(dim=c(iter,nq))
  
  #roll initial values
  g[1,] <- rexp(ng,lamG)
  #g[1,] <- g_init #use stating guess
  q[1,] <- rexp(nq,lamq)
  #q[1,] <- q_init
  
  #set values if starting point is fixed
  if(fixed_start == TRUE){
    g[1,] <- g_init
    q[1,] <- q_init
  }
  
  #initial value
  #find initial calculated H
  Hccalc <- HcCalc(g[1,],q[1,],G_adj,n)
  
  #this will break if H is not full
  #convert answer to dsCMatrix
  Hccalc <- (Hccalc+Matrix::t(Hccalc))/2
  #change H to dsCMatrix if it is not already
  #Hcalc <- as(Hcalc,"dsCMatrix") #unnecessary
  #extract values into vector
  #hcalc <- Hcalc@x
  hccalc <- as.vector(Hccalc[upper.tri(Hccalc,diag=TRUE)])
  
  #declare vector for log likelihoods
  lllh <- rep(NA,iter)
  lllhprimeg <- rep(NA,iter)
  lllhprimeq <- rep(NA,iter)
  lllh[1] <- loglikelihood_c(h,hccalc,sig2ep,g[1,],q[1,],lamG,lamq)
  #print(lllh)
  
  #or choose the starting value based on the best of some number of iid positions
  #if(fixed_start == FALSE){
  #  start <- mcStartSearch(h,G_adj,n,ng,nq,lamG,lamq,sig2ep,nguess=1000)
  #  g[1,] <- start$gmax
  #  q[1,] <- start$qmax
  #  lllh[1] <- start$lllh
  #}
  
  #indicator for number of rejections
  reject_g <- rep(0,iter)
  reject_q <- rep(0,iter)
  
  for(iteri in 2:iter){
    #propose a new set of parameters
    gprime <- abs(g[(iteri-1),] + rnorm(ng,0,sdG[iteri]))
    qprime <- abs(q[(iteri-1),] + rnorm(nq,0,sdq[iteri]))
    
    #above but allow to pass through negative
    #gprime <- g[(iteri-1),] + rnorm(ng,0,sdG[iteri])
    #qprime <- q[(iteri-1),] + rnorm(nq,0,sdq[iteri])
    
    #find calculated h from proposed move (only on g)
    Hccalcprimeg <- HcCalc(gprime,q[(iteri-1),],G_adj,n)
    #this will break if H is not full
    #convert answer to dsCMatrix
    Hccalcprimeg <- (Hccalcprimeg+Matrix::t(Hccalcprimeg))/2
    #change H to dsCMatrix if it is not already
    #Hcalcprimeg <- as(Hcalcprimeg,"dsCMatrix") unnecessary
    #extract values into vector
    #hcalcprimeg <- Hcalcprimeg@x
    hccalcprimeg <- as.vector(Hccalcprimeg[upper.tri(Hccalcprimeg,diag=TRUE)])
    
    #find likelihood of hcalcprime
    lllhprimeg[iteri] <- loglikelihood_c(h,hccalcprimeg,sig2ep,gprime,q[(iteri-1),],lamG,lamq)
    #print(paste0("lllhprime=",lllhprime))
    
    #determine next step according to MH algorithm
    if(lllhprimeg[iteri] > lllh[iteri-1]){
      g[iteri,] <- gprime
      q[iteri,] <- q[(iteri-1),]
      lllh[iteri] <- lllhprimeg[iteri]
    } else if(rbinom(1,1,exp(lllhprimeg[iteri]-lllh[iteri-1]))){ #accept at rate llhprime/llh when llhprime<llh
      g[iteri,] <- gprime
      q[iteri,] <- q[(iteri-1),]
      lllh[iteri] <- lllhprimeg[iteri]
    } else{                    #if not accepted, keep them the same for next iteration
      g[iteri,] <- g[(iteri-1),]
      q[iteri,] <- q[(iteri-1),]
      reject_g[iteri] <- 1 #set value to 1 if move was rejected
      lllh[iteri] <- lllh[iteri-1]
    }
    
    #check window every 100, then adjust proposal distribution variance as needed
    if(((iteri %% 100) == 0) && iteri != iter && fixed_sd == FALSE){
      if(sum(reject_g[iteri:(iteri-99)]) > 90){
        sdG[(iteri+1):iter] <- sdG[(iteri+1):iter]*.8
      } else if(sum(reject_g[iteri:(iteri-99)]) < 60){
        sdG[(iteri+1):iter] <- sdG[(iteri+1):iter]*1.2
      }
    }
    
    #find calculated h from proposed move (only on g)
    Hccalcprimeq <- HcCalc(g[iteri,],qprime,G_adj,n)
    #this will break if H is not full
    #convert answer to dsCMatrix
    Hccalcprimeq <- (Hccalcprimeq+Matrix::t(Hccalcprimeq))/2
    #change H to dsCMatrix if it is not already
    #Hcalcprimeq <- as(Hcalcprimeq,"dsCMatrix") #unecessary
    #extract values into vector
    #hcalcprimeq <- Hcalcprimeq@x
    hccalcprimeq <- as.vector(Hccalcprimeq[upper.tri(Hccalcprimeq,diag=TRUE)])
    
    #find likelihood of hcalcprime
    lllhprimeq[iteri] <- loglikelihood_c(h,hccalcprimeq,sig2ep,g[iteri,],qprime,lamG,lamq)
    #print(paste0("lllhprime=",lllhprime))
    
    #determine next step according to MH algorithm
    if(lllhprimeq[iteri] > lllh[iteri]){
      q[iteri,] <- qprime
      lllh[iteri] <- lllhprimeq[iteri]
    } else if(rbinom(1,1,exp(lllhprimeq[iteri]-lllh[iteri]))){ #accept at rate exp(lllhprime/lllh) when llhprime<llh
      q[iteri,] <- qprime
      lllh[iteri] <- lllhprimeq[iteri]
    } else{                    #if not accepted, keep them the same for next iteration
      q[iteri,] <- q[(iteri-1),]
      reject_q[iteri] <- 1 #set value to 1 if move was rejected
    }
    
    #check window every 100, then adjust proposal distribution variance as needed
    if(((iteri %% 100) == 0) && iteri != iter && fixed_sd == FALSE){
      if(sum(reject_q[iteri:(iteri-99)]) > 90){
        sdq[(iteri+1):iter] <- sdq[(iteri+1):iter]*.8
      } else if(sum(reject_q[iteri:(iteri-99)]) < 60){
        sdq[(iteri+1):iter] <- sdq[(iteri+1):iter]*1.2
      }
    }
    
  }
  
  #marginal means based on the last half of iterations
  g_marg <- colMeans(g[round(iter/2):iter,])
  q_marg <- colMeans(q[round(iter/2):iter,])
  
  #put marginal into matrix
  G <- G_adj
  G@x <- g_marg
  Matrix::diag(G) <- -Matrix::rowSums(G)
  
  #get rid of standard deviation vectors if they are constant
  if(fixed_sd == TRUE){
    sdG <- sdG[1]
    sdq <- sdq[1]
  }
  
  #find medians
  g_med <- matrixStats::colMedians(g)
  q_med <- matrixStats::colMedians(q)
  
  #find inter-quartile ranges
  g_iqr <- matrixStats::colIQRs(g)
  q_iqr <- matrixStats::colIQRs(q)
  
  return(list(G=G,g_marg=g_marg,q_marg=q_marg,g=g,q=q,g_med=g_med,q_med=q_med,lllh=lllh,lllhprimeg=lllhprimeg,lllhprimeq=lllhprimeq,
              reject_g=reject_g,reject_q=reject_q,sdG=sdG,sdq=sdq))
  
}

HcCalc <- function(g,q,G_adj,n){
  #declare generating matrix
  G <- G_adj
  
  #fill in transition rates
  G@x <- g
  #convert to non-sparse matrix
  G <- as.matrix(G) #this line appears to make code run 5x faster for 2x2 grid
  #fill diagonal
  Matrix::diag(G) <- -Matrix::rowSums(G)
  
  #declare hitting time matrix
  Hit <- matrix(0,nrow=n,ncol=n)
  
  #matrix of right sides
  b <- matrix(1,nrow=n,ncol=n)
  Matrix::diag(b) <- 0
  
  #inefficient
  for(i in 1:n){
    Hit[-i,i] <- solve(G[-i,-i],rep(-1,(n-1)))
  }
  
  #commute time
  C <- Hit + Matrix::t(Hit)
  
  #waiting times
  Q <- matrix(q,nrow=n,ncol=n) + matrix(q,nrow=n,ncol=n,byrow=TRUE)
  
  #coalescence approximation
  Hc <- C/4 + Q/2
  
  return(Hc)
}

loglikelihood_c <- function(h,hccalc,sig2ep,g,q,lamG,lamq){
  #hcalc is the calculated values of H from the current G and q
  #sig2ep is the variance on the error on the values of h
  
  #choose whether it bayesian or not (also change dU_an accordingly)
  #lllh <- -sum(lamG*g)-(1/(2*sig2ep))*sum((hccalc-h)^2,na.rm=TRUE) #bayesian
  lllh <- -sum(lamG*g)-sum(((hccalc-h)^2)/(2*sig2ep),na.rm=TRUE) #bayesian
  #lllh <- -(1/(2*sig2ep))*sum((hccalc-h)^2,na.rm=TRUE) #non-bayesian
  #lllh <- -sum(((hcalc-h)^2)/(2*sig2ep),na.rm=TRUE) #non-bayesian
  
  return(lllh)
}
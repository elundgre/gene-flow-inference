findG.HMC <- function(H,G_adj,const_coal=TRUE,iter=1000,fixed_start=FALSE,g_init=1,gam_init=1){
  #hamiltonian/hybrid monte carlo function
  #helped by https://theclevermachine.wordpress.com/2012/11/18/mcmc-hamiltonian-monte-carlo-a-k-a-hybrid-monte-carlo/
  #and http://www.cs.utoronto.ca/~radford/ham-mcmc-simple
  
  # make sure matrix is square and find dimentions
  if(NROW(H) == NCOL(H)) n <- NROW(H)
  else stop("Matrix is not square. Please enter a symmetric matrix")
  
  # make sure G and H have the same dimensions
  if((NROW(H) != NROW(G_adj)) || (NCOL(H) != NCOL(G))) stop("H and G_adj must have the same dimensions")
  
  #make sure matrix is symmetric up to computational error
  if(max(H-t(H),na.rm=TRUE)/max(H,na.rm=TRUE) > 1e-10) stop("Matrix is not symmetric. Please enter a symmetric matrix")
  #make H exactly symmetric
  H <- (H+t(H))/2
  #change H to dsCMatrix if it is not already
  #H <- as(H,"dsCMatrix") #don't need to do this
  #extract values into vector
  #h <- H@x
  h <- as.vector(H[upper.tri(H,diag=TRUE)])
  
  #make sure G_adj is in sparse matrix format and add diagonal for structure matrix
  G_adj <- as(G_adj,"dgCMatrix")
  #G_struct <- G_adj
  #diag(G_struct) <- 1
  #future: add check to make sure all non-zero entries are 1
  
  #find number of elements in g and gam
  ng <- length(G_adj@x)
  if(const_coal == TRUE) ngam <- 1 else ngam <- n
  
  #set up priors for G and gamma
  meanG <- 1;
  meangam <- 1;
  #convert to distribution parameters
  lamG <- 1/meanG
  lamgam <- 1/meangam
  #decide variance for error on values of H
  sig2ep <- (mean(H,na.rm=TRUE)/50)^2
  #set up parameters for the momentum
 # sdpG <- rep(meanG/200,iter)
 # sdpgam <- rep(meangam/200,iter)
  sdpG <- rep(1,iter)
  sdpgam <- rep(1,iter)
  
  #initialize arrays for g and gam
  g <- array(dim=c(iter,ng))
  gam <- array(dim=c(iter,ngam))
  
  #roll initial values
  g[1,] <- rexp(ng,lamG)
  #g[1,] <- g_init #use stating guess
  gam[1,] <- rexp(ngam,lamgam)
  #gam[1,] <- gam_init
  
  #set values if starting point is fixed
  if(fixed_start == TRUE){
    g[1,] <- g_init
    gam[1,] <- gam_init
  }
  
  #declare vector for log likelihoods
  lllh <- rep(NA,iter)
  lllhprime <- rep(NA,iter)
  #lllhprimeg <- rep(NA,iter)
  #lllhprimegam <- rep(NA,iter)
  #find first log likelihood
  Hl <- H_and_lllh(g[1,],gam[1,],G_adj,n,h,sig2ep,lamG,lamgam,const_coal)
  hcalc <- Hl$hcalc
  lllh[1] <- Hl$lllh
  #print(lllh)
  
  #or choose the starting value based on the best of some number of iid positions
  if(fixed_start == FALSE){
    start <- mcStartSearch(h,G_adj,n,ng,ngam,lamG,lamgam,sig2ep,const_coal,nguess=1000)
    g[1,] <- start$gmax
    gam[1,] <- start$gammax
    lllh[1] <- start$lllh
  }
  
  #indicator for number of rejections
  reject <- rep(0,iter)
  #reject_g <- rep(0,iter)
  #reject_gam <- rep(0,iter)
  
  #declare variables to store energies
  current_K <- rep(NA,iter)
  proposed_K <- rep(NA,iter)
  
  #step size for hamiltonian dynamics
  delta <- rep(1e-5,iter)
  #number of steps
  L <- 10
  
  #make array used in derivative of inverted matrix
  Deltas <- make_Deltas(G_adj,n,ng,ngam,const_coal)
  
  for(iteri in 2:iter){
    #put all parameters in a single vector
    gs <- c(g[(iteri-1),],gam[(iteri-1),])
    
    #sample random momentum
    p <- rep(NA,(ng+ngam))
    p[1:ng] <- rnorm(ng,0,sdpG[iteri])
    p[(ng+1):(ng+ngam)] <- rnorm(ngam,0,sdpgam[iteri])
    #save initial state
    current_p <- p
    
    #roll specific step size for this iteration
    deltatemp <- runif(1,delta[iteri]*.8,delta[iteri]*1.2)
    
    #for testing use
    if(FALSE){ #set to return after a specific step if you want
      #compare numerical and analytic methods after a certain step
      dU_n <- dU_num(lllh[iteri-1],gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,n,h,sig2ep,lamG,lamgam,delta,ng,ngam,const_coal)
      dU_a <- dU_an(lllh[iteri-1],gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,Deltas,n,h,hcalc,sig2ep,lamG,lamgam,ng,ngam,const_coal)
      
      d2U_a <- d2U_an(lllh[iteri-1],gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,Deltas,n,h,hcalc,sig2ep,lamG,lamgam,ng,ngam,const_coal)
      
      #test
      deltag <- 1e-5
      gs_new <- gs
      index1 <- 15
      gs_new[index1] <- gs_new[index1] + deltag
      Hl <- H_and_lllh(gs_new[1:ng],gs_new[(ng+1):(ng+ngam)],G_adj,n,h,sig2ep,lamG,lamgam,const_coal=TRUE)
      lllh_new <- Hl$lllh
      
      #print(-lllh_new + lllh[iteri-1])
      #print(-lllh_new + lllh[iteri-1] - deltag*dU_a[index1])
      #print(-lllh_new + lllh[iteri-1] - deltag*dU_a[index1] - (1/2)*(deltag^2)*d2U_a[index1,index1])
      
      print(eigen(d2U_a))
      
      return(list(dU_n=dU_n,dU_a=dU_a,d2U_an=d2U_an))
    }
    
    #use leap frog method to take L steps
    
    #first half step
    #p <- p - deltatemp/2*dU_num(lllh[iteri-1],gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,n,h,sig2ep,lamG,lamgam,delta,ng,ngam,const_coal)
    p <- p - deltatemp/2*dU_an(lllh[iteri-1],gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,
                                  Deltas,n,h,hcalc,sig2ep,lamG,lamgam,ng,ngam,const_coal)
    
    #full steps
    for(jL in 1:L){
      #full step for position
      gs <- gs + deltatemp*p
      #calculate new current log likelihood
      Hl <- H_and_lllh(gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,n,h,sig2ep,lamG,lamgam,const_coal)
      hcalc <- Hl$hcalc
      lllhstar <- Hl$lllh
      #full step for momentum
      #if(jL != L) p <- p - deltatemp*dU_num(lllhstar,gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,
      #                                  n,h,sig2ep,lamG,lamgam,delta,ng,ngam,const_coal)
      if(jL != L) p <- p - deltatemp*dU_an(lllhstar,gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,
                                              Deltas,n,h,hcalc,sig2ep,lamG,lamgam,ng,ngam,const_coal)
    }
    #last half step
    #p <- p - deltatemp/2*dU_num(lllhstar,gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,n,h,sig2ep,lamG,lamgam,delta,ng,ngam,const_coal)
    p <- p - deltatemp/2*dU_an(lllhstar,gs[1:ng],gs[(ng+1):(ng+ngam)],G_adj,
                                  Deltas,n,h,hcalc,sig2ep,lamG,lamgam,ng,ngam,const_coal)
    
    #record proposed log likelihood
    lllhprime[iteri] <- lllhstar
    
    #something about negating momemtum for symmetry even though it is discarded?
    
    #find rate at which it will be accepted
    current_U <- -lllh[(iteri-1)] #pretty sure there's a minus sign here?
    current_K[iteri] <- sum(current_p^2)/2
    proposed_U <- -lllhstar
    proposed_K[iteri] <- sum(p^2)/2
    
    #return(list(current_U,current_K,proposed_U,proposed_K))
    
    #choose whether to accept or reject
    if(runif(1) < exp(current_U-proposed_U+current_K[iteri]-proposed_K[iteri])){
      g[iteri,] <- gs[1:ng]
      gam[iteri,] <- gs[(ng+1):(ng+ngam)]
      lllh[iteri] <- lllhstar
    }
    else{
      g[iteri,] <- g[(iteri-1),]
      gam[iteri,] <- gam[(iteri-1),]
      lllh[iteri] <- lllh[(iteri-1)]
      reject[iteri] <- 1 #set value to 1 if move was rejected
    }
    
    #add momentum/step size change thing here?
    window <- 30 #window size
    if(((iteri %% window) == 0) && iteri != iter){
      if(sum(reject[iteri:(iteri-(window-1))]) >= round(.5*window)) delta[(iteri+1):iter] <- delta[iteri]*(2/3)
      else if(sum(reject[iteri:(iteri-9)]) <= round(.05*window)) delta[(iteri+1):iter] <- delta[iteri]*1.2
    }
    
  }
  
  #marginal ditributions based on the last half of iterations
  g_marg <- colMeans(g[round(iter/2):iter,])
  if(const_coal == TRUE) gam_marg <- mean(gam[round(iter/2):iter,])
  else gam_marg <- colMeans(gam[round(iter/2):iter,])
  
  #put marginal into matrix
  G <- G_adj
  G@x <- g_marg
  diag(G) <- -rowSums(G)
  
  return(list(G=G,g_marg=g_marg,gam_marg=gam_marg,g=g,gam=gam,lllh=lllh,lllhprime=lllhprime,reject=reject,
              current_K=current_K,proposed_K=proposed_K,gs=gs,delta=delta))
  
}

dU_num <- function(lllh_cur,g,gam,G_adj,n,h,sig2ep,lamG,lamgam,delta,ng,ngam,const_coal=TRUE){
  #step size to numerically find derivative
  eps <- 1e-5 #this should be small enough?
  
  dU <- rep(NA,(ng+ngam))
  
  for(i in 1:ng){
    #make new vector for shifted parameters
    gmod <- g
    gmod[i] <- g[i] + eps
    
    #calculate lllh from shifted
    Hl <- H_and_lllh(gmod,gam,G_adj,n,h,sig2ep,lamG,lamgam,const_coal)
    lllh_shift <- Hl$lllh

    #calculate numerical derivative
    dU[i] <- (-lllh_shift+lllh_cur)/eps
  }
  #now the same for gam
  for(i in 1:ngam){
    #make new vector for shifted parameters
    gammod <- gam
    gammod[i] <- gam[i] + eps
    
    #find new likelihood
    Hl <- H_and_lllh(g,gammod,G_adj,n,h,sig2ep,lamG,lamgam,const_coal)
    lllh_shift <- Hl$lllh
    
    #calculate numerical derivative
    dU[ng+i] <- (-lllh_shift+lllh_cur)/eps
  }
  
  return(dU)
}

dU_an <- function(lllh_cur,g,gam,G_adj,Deltas,n,h,hcalc,sig2ep,lamG,lamgam,ng,ngam,const_coal){
  
  dU <- rep(NA,(ng+ngam))
  
  #make and invert GG_c
  #declare generating matrix
  G <- G_adj
  
  #fill in transition rates
  G@x <- g
  #convert to non-sparse matrix
  G <- as.matrix(G) #this line appears to make code run 5x faster for 2x2 grid
  #fill diagonal
  diag(G) <- -rowSums(G)
  
  #now find the generating matrix for the product chain
  GG <- kronecker(G,diag(n)) + kronecker(diag(n),G)
  
  #make gamma into form where we can use it to make sub-markov chain
  if(const_coal==TRUE) Gam <- c(diag(rep(gam,n)))
  else Gam <- c(diag(gam))
  
  #turn generating matrix from product markov chain into one for sub-markov chain
  #to account for coalescence
  GG_c <- GG - diag(Gam)
  
  #invet GG_c
  GG_c_in <- solve(GG_c)
  
  #find dh for each g[i]
  
  #declare array for dh (shouldn't matter which index is which unless drop is true for subsetting?)
  dh <- array(dim=c((ng+ngam),length(h)))
  
  for(i in 1:(ng+ngam)){
    dH <- -GG_c_in %*% Deltas[[i]] %*% GG_c_in %*% rep(1,n^2)
    dH <- matrix(dH,n)
    
    dh[i,] <- dH[upper.tri(dH,diag=TRUE)]
    
    ##please don't let anyone see this
    #index <- 1
    #for(a in 1:n){
    #  for(b in 1:a){
    #    dh[i,index] <- dH[a,b]
    #    index <- index + 1
    #  }
    #}
  }
  
  #use lapply instead?
  #dh <- lapply()

  #find gradient of lllh w.r.t. hcalc
  del_lllh <- del_lllh(h,hcalc,sig2ep)
  
  for(i in 1:(ng+ngam)){
    #dU[i] <- del_lllh %*% dh[i,]
    dU[i] <- sum(del_lllh*dh[i,],na.rm=TRUE)
    #baysian part (do we really want this here?)
    if(i <= ng) dU[i] <- dU[i] +lamG
    else dU[i] <- dU[i] + lamgam
  }
  
  return(dU)
  
}

del_lllh <- function(h,hcalc,sig2ep){
  #gradient of lllh w.r.t. hcalc
  del_lllh <- -(1/sig2ep)*(hcalc-h)
  return(del_lllh)
}

make_Deltas <- function(G_adj,n,ng,ngam,const_coal){
  #this function determines how changing a parameter changes GG_c
  #slow, but only needs to run once, so it doesn't really matter
  
  #declare generating matrix
  G <- G_adj
  
  #don't need to fill in transition rates
  #convert to non-sparse matrix
  G <- as.matrix(G)
  #fill diagonal
  diag(G) <- -rowSums(G)
  
  #now find the generating matrix for the product chain
  GG <- kronecker(G,diag(n)) + kronecker(diag(n),G)
  
  #make gamma into form where we can use it to make sub-markov chain
  Gam <- c(diag(rep(1,n)))
  
  #turn generating matrix from product markov chain into one for sub-markov chain
  #to account for coalescence
  GG_c <- GG - diag(Gam)
  
  #Deltas <- array(dim=c(ng+ngam,n^2,n^2))
  Deltas <- as.list(rep(NA,(ng+ngam)))
  
  for(i in 1:ng){
    #say which value is changing
    gmod <- rep(1,ng)
    gmod[i] <- gmod[i] + 1
    #make matrix of which values change (super slow and stupid)
    Gmod <- G_adj
    #fill in matrix
    Gmod@x <- gmod
    #convert to non-sparse matrix
    Gmod <- as.matrix(Gmod)
    #fill diagonal
    diag(Gmod) <- -rowSums(Gmod)
    
    #make product
    GGmod <- kronecker(Gmod,diag(n)) + kronecker(diag(n),Gmod)
    
    #make gamma into form where we can use it to make sub-markov chain
    Gam <- c(diag(rep(1,n)))
    
    #turn generating matrix from product markov chain into one for sub-markov chain
    #to account for coalescence
    GG_cmod <- GGmod - diag(Gam)
    
    #find differences (values should be -1, 0, or 1)
    #Deltas[i,,] <- GG_cmod - GG_c
    Deltas[[i]] <- GG_cmod - GG_c
  }
  #now for the gammas
  for(i in 1:ngam){
    #initialize gamma
    gammod <- rep(1,ngam)
    gammod[i] <- gammod[i] + 1
    
    #can start after the kronecker product
    GGmod <- GG
    
    #make gamma into form where we can use it to make sub-markov chain
    if(const_coal==TRUE) Gam <- c(diag(rep(gammod,n)))
    else Gam <- c(diag(gammod))
    
    #turn generating matrix from product markov chain into one for sub-markov chain
    #to account for coalescence
    GG_cmod <- GGmod - diag(Gam)
    
    #find differences (values should be -1, 0, or 1)
    #Deltas[(ng+i),,] <- GG_cmod - GG_c
    Deltas[[(ng+i)]] <- GG_cmod - GG_c
  }
  
  return(Deltas)
}

H_and_lllh <- function(g,gam,G_adj,n,h,sig2ep,lamG,lamgam,const_coal=TRUE){
  #find Hcalc using HCalc
  Hcalc <- HCalc(g,gam,G_adj,n,const_coal)
  
  #this will break if H is not full
  #convert answer to dsCMatrix
  Hcalc <- (Hcalc+t(Hcalc))/2
  #change H to dsCMatrix if it is not already
  #Hcalc <- as(Hcalc,"dsCMatrix") #unnecessary
  #extract values into vector
  #hcalc <- Hcalc@x
  hcalc <- as.vector(Hcalc[upper.tri(Hcalc,diag=TRUE)])
  
  #find log likelihood
  lllh <- loglikelihood(h,hcalc,sig2ep,g,gam,lamG,lamgam)
  
  return(list(hcalc=hcalc,lllh=lllh))
}

HCalc <- function(g,gam,G_adj,n,const_coal=TRUE){
  #declare generating matrix
  G <- G_adj
  
  #fill in transition rates
  G@x <- g
  #convert to non-sparse matrix
  G <- as.matrix(G) #this line appears to make code run 5x faster for 2x2 grid
  #fill diagonal
  diag(G) <- -rowSums(G)
  
  #now find the generating matrix for the product chain
  GG <- kronecker(G,diag(n)) + kronecker(diag(n),G)
  
  #make gamma into form where we can use it to make sub-markov chain
  if(const_coal==TRUE) Gam <- c(diag(rep(gam,n)))
  else Gam <- c(diag(gam))
  
  #turn generating matrix from product markov chain into one for sub-markov chain
  #to account for coalescence
  GG_c <- GG - diag(Gam)
  
  #calculate coalescence times
  hcalc <- solve(GG_c,rep(-1,n^2))
  #write hcalc as matrix
  Hcalc <- matrix(hcalc,nrow=n)
  
  return(Hcalc)
}

loglikelihood <- function(h,hcalc,sig2ep,g,gam,lamG,lamgam){
  #hcalc is the calculated values of H from the current G and gam
  #sig2ep is the variance on the error on the values of h
  
  #disallow negative values
  if(!all(g>0,gam>0)) return(-Inf)
  
  #choose whether it bayesian or not (also change dU_an accordingly)
  lllh <- -sum(lamG*g)-sum(lamgam*gam)-(1/(2*sig2ep))*sum((hcalc-h)^2,na.rm=TRUE) #bayesian
  #lllh <- -(1/(2*sig2ep))*sum((hcalc-h)^2,na.rm=TRUE) #non-bayesian
  
  return(lllh)
}


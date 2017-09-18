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
  if(max(H-Matrix::t(H),na.rm=TRUE)/max(H,na.rm=TRUE) > 1e-10) stop("Matrix is not symmetric. Please enter a symmetric matrix")
  #make H exactly symmetric
  H <- (H+Matrix::t(H))/2
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
  diag(G) <- -Matrix::rowSums(G)
  
  return(list(G=G,g_marg=g_marg,gam_marg=gam_marg,g=g,gam=gam,lllh=lllh,lllhprime=lllhprime,reject=reject,
              current_K=current_K,proposed_K=proposed_K,gs=gs,delta=delta))
  
}


#mcmc.fns 
#file for auxiliary functions

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

mcStartSearch <- function(h,G_adj,n,ng,ngam,lamG,lamgam,sig2ep,const_coal,nguess=1000){
  #choose some number of random starting points and pick the best one
  #currently unsure about what to do about gamma
  
  #nguess <- 1000 #number of guesses to test
  
  #find mean of gamma
  #meangam <- mean(gam)
  
  #roll parameter arrays
  gs <- array(data=rexp(nguess*ng,rate=lamG),dim=c(nguess,ng))
  #roll gammas between .5 and 1.5 of current value
  #gams <- array(runif(nguess*ngam,meangam*.5,meangam*1.5),dim=c(nguess,ngam))
  #roll randomly
  gams <- array(rexp(nguess*ngam,1),dim=c(nguess,ngam))
  #declare likelihood array
  lllh <- rep(NA,nguess)
  
  for(i in 1:nguess){
    Hl <- H_and_lllh(gs[i,],gams[i,],G_adj,n,h,sig2ep,lamG,lamgam,const_coal)
    lllh[i] <- Hl$lllh
  }
  
  #find maximum
  index <- which.max(lllh)
  
  gmax <- gs[index,]
  gammax <- gams[index,]
  
  return(list(gmax=gmax,gammax=gammax,lllh=lllh[index]))
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
  
  #disallow negative values (only possible with HMC)
  if(!all(g>0,gam>0)) return(-Inf)
  
  #choose whether it bayesian or not (also change dU_an accordingly)
  lllh <- -sum(lamG*g)-sum(lamgam*gam)-(1/(2*sig2ep))*sum((hcalc-h)^2,na.rm=TRUE) #bayesian
  #lllh <- -(1/(2*sig2ep))*sum((hcalc-h)^2,na.rm=TRUE) #non-bayesian
  
  return(lllh)
}
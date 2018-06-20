#function for running mcmc stuff when you don't know the actual answer

run.mcmc <- function(width,height,g_known=TRUE,fixed_g=FALSE,g=1,coal_type=3,gam=1,const_coal=TRUE,noise=1/500,seed=sample(1000000000,1),
                     H=NA,h_se=NA,preburn_iter=1e5,burn_iter=1e5,iter=1e5,continue=FALSE,g_init,gam_init,sdG_init,sdgam_init,
                     noisy_H=FALSE,missing_H=FALSE,missing_loc,type="coal"){
  #coal_type: 1 is constant everywhere, 2 is iid, 3 is 1 everywhere, 4 is user specified
  #noise: st.dev for noise added to H relative to H. e.g. noise of .01 adds noise with stdev mean(H)*.01
  
  #require(Matrix)
  #require(igraph)
  
  #source(file.path(path,'mcmc.fns.R'))
  #source(file.path(path,'findG.HMC.R'))
  #source(file.path(path,'findG.MH.R'))
  #source(file.path(path,'findG.MH.com.R'))
  #source(file.path(path,'findG.nnls.R'))
  
  #save RNG state
  set.seed(seed)
  
  #number of states in the system
  n <- width*height
  
  #declare generating matrix
  G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height))) 
  G <- G_adj
  
  #find number of parameters
  ng <- length(G_adj@x)
  
  #find H based on g and gamma
  if(g_known == TRUE){
  
    #roll transition rates
    if(fixed_g == FALSE){
      g <- rexp(ng,1)
      #round up to nearest 10th
      g <- .1*ceiling(10*g)
    }
    
    #put rates in matrix
    G@x <- g
    
    #fill diagonal
    Matrix::diag(G) <- -Matrix::rowSums(G)
    
    #now find the generating matrix for the product chain
    GG <- Matrix::kronecker(G,Matrix::diag(n)) + Matrix::kronecker(Matrix::diag(n),G)
    
    #roll coalescence rates
    #1 is constant everywhere, 2 is iid, 3 is 1 everywhere, 4 is user specified
    if(coal_type == 1){
      gam <- rep(rexp(1,1),n)+.5
      gam <- .1*ceiling(10*gam)
      const_coal = TRUE
      ngam <- 1
    } else if(coal_type == 2){
      gam <- rexp(n,1)+.5
      gam <- .1*ceiling(10*gam)
      const_coal = FALSE
      ngam <- n
    } else if(coal_type == 3){
      gam <- rep(1,n)
      if(const_coal) ngam <- 1 else ngam <- n
    } else if(coal_type == 4){
      if(length(gam) != n) stop("invalid number of coalescence rates")
      if(const_coal) ngam <- 1 else ngam <- n
    }
    else{
      stop("invalid coalescence type")
    }
    
    #make gamma into form where we can use it to make sub-markov chain
    Gam <- c(Matrix::diag(gam))
    
    #turn generating matrix from product markov chain into one for sub-markov chain
    #to account for coalescence
    GG_c <- GG - Matrix::diag(Gam)
    
    #now solve for the coalescence times
    system.time(H <- Matrix::solve(GG_c,rep(-1,n^2)))
    
    #write H as matrix
    H <- matrix(H,nrow=n)
    H_infer <- H
    #print(H)
    
    #no scaling if g is known (may change later)
    scale <- 1
    
    #add noise to H
    sigep <- mean(H)*noise #stdev
    sig2ep <- sigep^2 #variance
    
    H_noise <- H + matrix(rnorm(n^2,mean=0,sd=sigep),n)
    #make it symmetric (replace lower triangle  with upper triangle)
    H_noise[lower.tri(H_noise)] <- t(H_noise)[lower.tri(H_noise)]
    if(noisy_H==TRUE) H_infer <- H_noise
  
  }
  else{
    H_infer <- H
    #use provided noise estimate (if provided)
    #if(h_se != NA){
    sigep <- h_se #stdev (standard error of mean in h)
    
    #scale input so coalescence parameter is of about order 1
    #coalescence should be approximately the number of states
    scale <- n/mean(H_infer)
    
    H_infer <- H_infer*scale
    sigep <- sigep*scale
    
    sig2ep <- sigep^2 #variance
    #}
    #else{ #use noise value if not
    #  sigep <- mean(H)*noise #stdev
    #  sig2ep <- sigep^2 #variance
    #}
  }
    
  #calculate G and gam from H (non-negative least squares solution)
  time <- system.time(ans2 <- findG.nnls(H,G_adj,const_coal))
  
  if(g_known == TRUE){
    print(paste0("Max difference between calculated and true G: ",max(abs(ans2$G-G))))
    print(paste0("Max difference between calculated and true gamma: ",max(abs(ans2$gam-gam))))
  
  
    #calculate G and gam from H_noise (non-negative least squares solution)
    system.time(ans3 <- findG.nnls(H_noise,G_adj,const_coal))
    
    print(max(abs(ans2$G-ans3$G)))
  }
  
  if(noisy_H==FALSE && missing_H==TRUE){ #make H2 with missing H values
    H2 <- H
    H2[missing_loc,] <- NA
    H2[,missing_loc] <- NA
    H_infer <- H2
  }
  
  if(noisy_H==TRUE && missing_H==TRUE){  #make H2_noise
    H2_noise <- H_noise
    H2_noise[missing_loc,] <- NA
    H2_noise[,missing_loc] <- NA
    H2_noise <- as.matrix(H2_noise)
    H_infer <- H2_noise
  }
  
  
  h <- as.vector(H[upper.tri(H,diag=TRUE)])
  if(g_known==TRUE){
    h_noise <- as.vector(H_noise[upper.tri(H_noise,diag=TRUE)])
  }
  
  if(continue==FALSE){
    if(type=="coal"){
      #pre-burn with not very strict likelihood
      time_pb <- system.time(pburn <- findG.MH(H_infer,G_adj,const_coal,iter=preburn_iter,sig2ep=sig2ep*10))
    }
    if(type=="com"){
      #pre-burn with not very strict likelihood
      time_pb <- system.time(pburn <- findG.MH.com(H_infer,G_adj,iter=preburn_iter,sig2ep=sig2ep*10))
    }
    if(type=="coal"){
      time_b <- system.time(burn <- findG.MH(H_infer,G_adj,const_coal,iter=burn_iter,sig2ep=sig2ep,
                                             fixed_start=TRUE,g_init=pburn$g[preburn_iter,],gam_init=pburn$gam[preburn_iter,]))
    }
    else if(type=="com"){
      time_b <- system.time(burn <- findG.MH.com(H_infer,G_adj,iter=burn_iter,sig2ep=sig2ep,
                                                 fixed_start=TRUE,g_init=pburn$g[preburn_iter,],q_init=pburn$q[preburn_iter,]))
    }
    maxval <- burn_iter
    minval <- round(burn_iter/2)
    g_init <- burn$g[maxval,]
    sdG_init <- mean(burn$sdG[minval:maxval])
    if(type=="coal"){
      gam_init <- burn$gam[maxval,]
      sdgam_init <- mean(burn$sdgam[minval:maxval])
    }
    else if(type=="com"){
      q_init <- burn$q[maxval,]
      sdq_init <- mean(burn$sdq[minval:maxval])
    }
  }
  
  if(type=="coal"){
    time <- system.time(ans <- findG.MH(H_infer,G_adj,const_coal,iter=iter,fixed_start=TRUE,g_init=g_init,gam_init=gam_init,
                                        fixed_sd=TRUE,fixed_sd_start=TRUE,sdG_init=sdG_init,sdgam_init=sdgam_init,sig2ep=sig2ep))
  }
  else if(type=="com"){
    time <- system.time(ans <- findG.MH.com(H_infer,G_adj,iter=iter,fixed_start=TRUE,g_init=g_init,q_init=q_init,
                                            fixed_sd=TRUE,fixed_sd_start=TRUE,sdG_init=sdG_init,sdq_init=sdq_init,sig2ep=sig2ep))
  }
  
  if(continue==FALSE) out <- list(ans=ans,pburn=pburn,burn=burn,seed=seed,sig2ep=sig2ep,
                                  time_pb=time_pb,time_b=time_b,time=time,iter=iter,H=H,h=h,H_infer=H_infer,scale=scale)
  else out <- list(ans=ans,seed=seed,sig2ep=sig2ep,time=time,iter=iter,H=H,h=h,H_infer,scale=scale)
  return(out)
}

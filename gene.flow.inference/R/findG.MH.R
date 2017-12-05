findG.MH <- function(H,G_adj,const_coal=TRUE,iter=10000,fixed_start=FALSE,g_init=1,gam_init=1,
                      fixed_sd=FALSE,fixed_sd_start=FALSE,sdG_init=1/200,sdgam_init=1/200,
                      sig2ep=(mean(H,na.rm=TRUE)/50)^2){
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
  #sig2ep <- (mean(H,na.rm=TRUE)/25)^2
  #set up parameters for the proposal distribution
  if(fixed_sd_start == FALSE){
    sdG <- rep(meanG/200,iter)
    sdgam <- rep(meangam/200,iter)
  }
  else{
    sdG <- rep(sdG_init,iter)
    sdgam <- rep(sdgam_init,iter)
  }
  
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
  
  #initial value
  #find initial calculated H
  Hcalc <- HCalc(g[1,],gam[1,],G_adj,n,const_coal)
  
  #this will break if H is not full
  #convert answer to dsCMatrix
  Hcalc <- (Hcalc+Matrix::t(Hcalc))/2
  #change H to dsCMatrix if it is not already
  #Hcalc <- as(Hcalc,"dsCMatrix") #unnecessary
  #extract values into vector
  #hcalc <- Hcalc@x
  hcalc <- as.vector(Hcalc[upper.tri(Hcalc,diag=TRUE)])
  
  #declare vector for log likelihoods
  lllh <- rep(NA,iter)
  lllhprimeg <- rep(NA,iter)
  lllhprimegam <- rep(NA,iter)
  lllh[1] <- loglikelihood(h,hcalc,sig2ep,g[1,],gam[1,],lamG,lamgam)
  #print(lllh)
  
  #or choose the starting value based on the best of some number of iid positions
  if(fixed_start == FALSE){
    start <- mcStartSearch(h,G_adj,n,ng,ngam,lamG,lamgam,sig2ep,const_coal,nguess=1000)
    g[1,] <- start$gmax
    gam[1,] <- start$gammax
    lllh[1] <- start$lllh
  }
  
  #indicator for number of rejections
  reject_g <- rep(0,iter)
  reject_gam <- rep(0,iter)
  
  for(iteri in 2:iter){
    #propose a new set of parameters
    gprime <- abs(g[(iteri-1),] + rnorm(ng,0,sdG[iteri]))
    gamprime <- abs(gam[(iteri-1),] + rnorm(ngam,0,sdgam[iteri]))
    
    #above but allow to pass through negative
    #gprime <- g[(iteri-1),] + rnorm(ng,0,sdG[iteri])
    #gamprime <- gam[(iteri-1),] + rnorm(ngam,0,sdgam[iteri])
    
    #find calculated h from proposed move (only on g)
    Hcalcprimeg <- HCalc(gprime,gam[(iteri-1),],G_adj,n,const_coal)
    #this will break if H is not full
    #convert answer to dsCMatrix
    Hcalcprimeg <- (Hcalcprimeg+Matrix::t(Hcalcprimeg))/2
    #change H to dsCMatrix if it is not already
    #Hcalcprimeg <- as(Hcalcprimeg,"dsCMatrix") unnecessary
    #extract values into vector
    #hcalcprimeg <- Hcalcprimeg@x
    hcalcprimeg <- as.vector(Hcalcprimeg[upper.tri(Hcalcprimeg,diag=TRUE)])
    
    #find likelihood of hcalcprime
    lllhprimeg[iteri] <- loglikelihood(h,hcalcprimeg,sig2ep,gprime,gam[(iteri-1),],lamG,lamgam)
    #print(paste0("lllhprime=",lllhprime))
    
    #determine next step according to MH algorithm
    if(lllhprimeg[iteri] > lllh[iteri-1]){
      g[iteri,] <- gprime
      gam[iteri,] <- gam[(iteri-1),]
      lllh[iteri] <- lllhprimeg[iteri]
    } else if(rbinom(1,1,exp(lllhprimeg[iteri]-lllh[iteri-1]))){ #accept at rate llhprime/llh when llhprime<llh
      g[iteri,] <- gprime
      gam[iteri,] <- gam[(iteri-1),]
      lllh[iteri] <- lllhprimeg[iteri]
    } else{                    #if not accepted, keep them the same for next iteration
      g[iteri,] <- g[(iteri-1),]
      gam[iteri,] <- gam[(iteri-1),]
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
    Hcalcprimegam <- HCalc(g[iteri,],gamprime,G_adj,n,const_coal)
    #this will break if H is not full
    #convert answer to dsCMatrix
    Hcalcprimegam <- (Hcalcprimegam+Matrix::t(Hcalcprimegam))/2
    #change H to dsCMatrix if it is not already
    #Hcalcprimegam <- as(Hcalcprimegam,"dsCMatrix") #unecessary
    #extract values into vector
    #hcalcprimegam <- Hcalcprimegam@x
    hcalcprimegam <- as.vector(Hcalcprimegam[upper.tri(Hcalcprimegam,diag=TRUE)])
    
    #find likelihood of hcalcprime
    lllhprimegam[iteri] <- loglikelihood(h,hcalcprimegam,sig2ep,g[iteri,],gamprime,lamG,lamgam)
    #print(paste0("lllhprime=",lllhprime))
    
    #determine next step according to MH algorithm
    if(lllhprimegam[iteri] > lllh[iteri]){
      gam[iteri,] <- gamprime
      lllh[iteri] <- lllhprimegam[iteri]
    } else if(rbinom(1,1,exp(lllhprimegam[iteri]-lllh[iteri]))){ #accept at rate exp(lllhprime/lllh) when llhprime<llh
      gam[iteri,] <- gamprime
      lllh[iteri] <- lllhprimegam[iteri]
    } else{                    #if not accepted, keep them the same for next iteration
      gam[iteri,] <- gam[(iteri-1),]
      reject_gam[iteri] <- 1 #set value to 1 if move was rejected
    }
    
    #check window every 100, then adjust proposal distribution variance as needed
    if(((iteri %% 100) == 0) && iteri != iter && fixed_sd == FALSE){
      if(sum(reject_gam[iteri:(iteri-99)]) > 90){
        sdgam[(iteri+1):iter] <- sdgam[(iteri+1):iter]*.8
      } else if(sum(reject_gam[iteri:(iteri-99)]) < 60){
        sdgam[(iteri+1):iter] <- sdgam[(iteri+1):iter]*1.2
      }
    }
    
  }
  
  #marginal means based on the last half of iterations
  g_marg <- colMeans(g[round(iter/2):iter,])
  if(const_coal == TRUE) gam_marg <- mean(gam[round(iter/2):iter,])
  else gam_marg <- colMeans(gam[round(iter/2):iter,])
  
  #put marginal into matrix
  G <- G_adj
  G@x <- g_marg
  Matrix::diag(G) <- -Matrix::rowSums(G)
  
  #get rid of standard deviation vectors if they are constant
  if(fixed_sd == TRUE){
    sdG <- sdG[1]
    sdgam <- sdgam[1]
  }
  
  #find medians
  g_med <- matrixStats::colMedians(g)
  gam_med <- matrixStats::colMedians(gam)
  
  return(list(G=G,g_marg=g_marg,gam_marg=gam_marg,g=g,gam=gam,g_med=g_med,gam_med=gam_med,lllh=lllh,lllhprimeg=lllhprimeg,
              lllhprimegam=lllhprimegam,reject_g=reject_g,reject_gam=reject_gam,sdG=sdG,sdgam=sdgam,H=H,h=h))
  
}
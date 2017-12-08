  #function for running many small graphs in parallel for commute vs coal comparison
  mult.small <- function(noise,gam,const_coal=TRUE,seed=sample(1000000000,1),symmetric,pb_iter_com=2e4,b_iter_com=8e4,iter_com=1e5,
                         pb_iter_coal=2e4,b_iter_coal=8e4,iter_coal=1e5){
  
  #require(Matrix)
  #require(igraph)
  #require(matrixStats)
    
  #source(paste0(path,'mcStartSearch.R'))
  #source(paste0(path,'findG.HMC.R'))
  #source(paste0(path,'findG.MH.R'))
  #source(paste0(path,'findG.MH.com.R'))
  #source(paste0(path,'findG.nnls.R'))
  #source(paste0(path,'run.mcmc.R'))
  
  #seed <- sample(1000000000,1) #now in top
  set.seed(seed)
  
  #if we are using a grid, specify the dimentions
  width <- 3
  height <- 2
  
  #number of states in the system
  n <- width*height
  
  G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height)))
  G <- G_adj
  #make into 
  if(symmetric==TRUE){
    G <- as(G,"symmetricMatrix")
  }
  
  #find number of parameters
  ng <- length(G_adj@x)
  
  #or make more random
  #roll transition rates
  if(symmetric==TRUE){
    g <- rexp(ng/2,1) #only roll half as many movement parameters if it is symmetric
  } else{
    g <- rexp(ng,1)
  }
  #round up to nearest 10th
  g <- .1*ceiling(10*g)
  
  #put values in matrix
  G@x <- g
  G <- as(G,"dgCMatrix")
  g <- G@x   #get all values out
  
  #fill diagonal
  Matrix::diag(G) <- -Matrix::rowSums(G)
  
  gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
  igraph::E(gr)$curved <- TRUE
  plot(gr,layout=as.matrix(expand.grid(1:width,1:height)),edge.width=g,edge.label=paste0("g",1:ng,"=",round(g*1000)/1000),
       main="3x2 Symmetric Graph")
  
  #set.seed(sample(1000000,1))
  
  #expand gamma for use in function
  if(length(gam)==1){
    gam <- rep(gam,n)
  }
  if(length(gam) != n){
    stop("Supplied gamma must be a scalar or a vector of length width x heght")
  }
  
  a3x2_com <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=gam,const_coal=const_coal,noise=noise,seed=seed,
                        preburn_iter=pb_iter_com,burn_iter=b_iter_com,iter=iter_com,noisy_H=TRUE,type="com")
  
  #boxplot(a3x2_com$ans$g[10*(1:(iter_com/10)),order(g)],outline=FALSE,main="3x2 Symmetric Graph (Commute Time Inference)",
  #        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  #points(g[order(g)],pch=19,col=3)
  
  
  a3x2_coal <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=gam,const_coal=const_coal,noise=noise,seed=seed,
                         preburn_iter=pb_iter_coal,burn_iter=b_iter_coal,iter=iter_coal,noisy_H=TRUE,type="coal")
  
  #boxplot(a3x2_coal$ans$g[10*(1:(iter_coal/10)),order(g)],outline=FALSE,main="3x2 Symmetric Graph (Coalescence Time Inference)",
  #        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  #points(g[order(g)],pch=19,col=3)

  return(list(g=g,gam=gam,G=G,seed=seed,noise=noise,a3x2_com=a3x2_com,a3x2_coal=a3x2_coal,
              width=width,height=height,n=n,G_adj=G_adj,pb_iter_com=pb_iter_com,b_iter_com=b_iter_com,iter_com=iter_com,
              pb_iter_coal=pb_iter_coal,b_iter_coal=b_iter_coal,iter_coal=iter_coal))
  }
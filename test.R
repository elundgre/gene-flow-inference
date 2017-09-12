#testing
#enter file path (put / at end)
path <- '~/Documents/gene-flow-inference/'
source(paste0(path,'run.mcmc.R'))

require(Matrix)
require(igraph)
require(matrixStats)

seed <- sample(1000000000,1)
set.seed(seed)

#if we are using a grid, specify the dimentions
width <- 3
height <- 2

#number of states in the system
n <- width*height

G_adj <- get.adjacency(make_lattice(c(width,height)))
G <- G_adj
#make into 
G <- as(G,"symmetricMatrix")

#find number of parameters
ng <- length(G_adj@x)

#roll transition rates
#g <- sample(c(rep(.5,4),rep(2,3)))

#or make more random
#roll transition rates
g <- rexp(ng/2,1)
#round up to nearest 10th
g <- .1*ceiling(10*g)

#put values in matrix
G@x <- g
G <- as(G,"dgCMatrix")
g <- G@x   #get all values out

gr <- graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
E(gr)$curved <- TRUE
plot(gr,layout=as.matrix(expand.grid(1:width,1:height)),edge.width=g,edge.label=paste0("g",1:ng,"=",round(g*1000)/1000),
     main="3x2 Symmetric Graph")

preburn_iter <- 2e4
burn_iter <- 4e4
iter <- 1e5

a3x2s_com <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=rep(1,n),const_coal=TRUE,noise=1/500,seed=seed,
                      preburn_iter=preburn_iter,burn_iter=burn_iter,iter=iter,noisy_H=TRUE,type="com",path=path)

boxplot(a3x2s_com$ans$g[10*(1:(iter/10)),order(g)],outline=FALSE,main="3x2 Symmetric Graph (Commute Time Inference)",
        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
points(g[order(g)],pch=19,col=3)

preburn_iter_coal <- 2e4
burn_iter_coal <- 4e4
iter_coal <- 1e5

a3x2s_coal <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=rep(1,n),const_coal=TRUE,noise=1/500,seed=seed,
                       preburn_iter=preburn_iter_coal,burn_iter=burn_iter_coal,iter=iter_coal,noisy_H=TRUE,type="coal",path=path)

boxplot(a3x2s_coal$ans$g[10*(1:(iter_coal/10)),order(g)],outline=FALSE,main="3x2 Symmetric Graph (Coalescence Time Inference)",
        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
points(g[order(g)],pch=19,col=3)

a3x2s_coal_hmc <- findG.HMC(a3x2s_coal$H,G_adj,const_coal=TRUE,iter=10000,fixed_start=FALSE)

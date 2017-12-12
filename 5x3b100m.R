#script for 5x3 with barriers, noise equal to 1/100 mean of H, with missing values and unknown coalescence

require(parallel)
require(gene.flow.inference)

#file path for output
#out_path <- "~/Dropbox/code/gene-flow-inference/plots/"
#or use current path
out_path <- getwd()

seed <- 129503301 #chosen from random.org uniform 1 to 1000000000
set.seed(seed)

#if we are using a grid, specify the dimentions
width <- 5
height <- 3

#number of states in the system
n <- width*height

#find adjacency matrix
G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height)))
G <- G_adj

#find number of parameters
ng <- length(G_adj@x)

#roll transition rates
g <- rexp(ng,1)

#add barriers (for 5x3 grid)
g[c(17,5,21,8,37,24,40,28)] <- 0

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
     main="5x3 Graph with Barriers")

#specify coalescence type
const_coal <- FALSE
#specify coalescence rates
gam <- rep(1,n)

pb_iter_coal <- 1e6
b_iter_coal <- 4e6
iter_coal <- 5e6

#specify amount of noise
noise <- 1/100

#specify that there are missing locations
missing_H <- TRUE
#specify locations
missing_loc <- c(2,5,7,15)

a_coal <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=gam,const_coal=const_coal,noise=noise,seed=seed,
                   preburn_iter=pb_iter_coal,burn_iter=b_iter_coal,iter=iter_coal,noisy_H=TRUE,type="coal",
                   missing_H = missing_H,missing_loc=missing_loc)

save.image("5x3bm.Rdata")

boxplot(a_coal$ans$g[1000*(1:iter_coal/1000),order(g)],outline=FALSE,
        main="5x3 Graph with Barriers, Noise 1/100, 4 Missing Values",
        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
points(g[order(g)],pch=19,col=3)

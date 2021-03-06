#4x4 asymmetric biased
require(gene.flow.inference)

seed <- sample(1000000000,1)
set.seed(seed)

#if we are using a grid, specify the dimentions
width <- 4
height <- 4

#number of states in the system
n <- width*height

G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height)))
G <- G_adj

#find number of parameters
ng <- length(G_adj@x)

#specify values of g
g <- rep(.5,ng)
g[c(1,4,7,12,16,20,26,30,34,40,43,46,2,5,8,10,13,17,21,24,27,31,35,38)] <- 2

#put rates in matrix
G@x <- g

#fill diagonal
diag(G) <- -Matrix::rowSums(G)

#now find the generating matrix for the product chain
GG <- Matrix::kronecker(G,diag(n)) + Matrix::kronecker(diag(n),G)

#specify gam
gam <- rep(1,n)

#make gamma into form where we can use it to make sub-markov chain
Gam <- c(diag(gam))

#turn generating matrix from product markov chain into one for sub-markov chain
#to account for coalescence
GG_c <- GG - diag(Gam)

#now solve for the coalescence times
H <- Matrix::solve(GG_c,rep(-1,n^2))
#write H as matrix
H <- matrix(H,nrow=n)

#solve for commute times
Hc <- HcCalc(g,Matrix::diag(H),G_adj,n)

plot(H,Hc,main="Commute Time + Diversity vs Coalescence Time (Biased Asymmetric Graph)",
     xlab="Coalescence Time",ylab="Commute Time + Diversity")
abline(0,1) #add y=x line

gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
plot(gr,layout=as.matrix(expand.grid(1:width,1:height)),edge.width=g,edge.label=paste0("g",1:ng,"=",round(g*1000)/1000),
     main="3x2 Symmetric Graph")

preburn_iter <- 1e5
burn_iter <- 1e6
iter <- 4e6

#set.seed(sample(1000000,1))

a4x4asb_com <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=gam,const_coal=TRUE,noise=1/500,seed=seed,
                       preburn_iter=preburn_iter,burn_iter=burn_iter,iter=iter,noisy_H=TRUE,type="com")

boxplot(a4x4asb_com$ans$g[10*(1:(iter/10)),order(g)],outline=FALSE,main="4x4 Biased Asymmetric Graph (Commute Time Inference)",
        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
points(g[order(g)],pch=19,col=3)

preburn_iter_coal <- 1e5
burn_iter_coal <- 1e6
iter_coal <- 4e6

a4x4asb_coal <- run.mcmc(width,height,fixed_g=TRUE,g=g,coal_type=4,gam=gam,const_coal=TRUE,noise=1/500,seed=seed,
                        preburn_iter=preburn_iter_coal,burn_iter=burn_iter_coal,iter=iter_coal,noisy_H=TRUE,type="coal")

boxplot(a4x4asb_coal$ans$g[10*(1:(iter_coal/10)),order(g)],outline=FALSE,main="4x4 Biased Asymmetric (Coalescence Time Inference)",
        names=paste0("g",order(g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
points(g[order(g)],pch=19,col=3)

#put everything in list
a4x4asb <- list(coal=a4x4asb_coal,com=a4x4asb_com,H=H,Hc=Hc,G_adj=G_adj,g=g,gam=gam,burn_iter=burn_iter,burn_iter_coal=burn_iter_coal,
                iter=iter,iter_coal=iter_coal,preburn_iter_coal=preburn_iter_coal,preburn_iter=preburn_iter,seed=seed,width=width,height=height)

save(a4x4asb,file="a4x4asb.RData")

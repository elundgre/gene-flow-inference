#script for inference on tortoise data
mylib <- "/home/cmb-00/plr/elundgre/R_packages/"
.libPaths(c(.libPaths(),mylib))

require(gene.flow.inference)
require(Matrix)

#select which watershed resolution to use (even number 2 through 8)
#WBD_res <- 8

#import adjacency matrix
#G_adj <- read.csv(paste0("WBD",WBD_res,"_adjacency.csv"))

#use handpicked?
G_adj <- read.csv("handpicked_WBD8_adjacency.csv")

#remove first column
G_adj <- G_adj[,-1]
#make into sparse matrix
G_adj <- as(as.matrix(G_adj),"dgCMatrix")

#find number of movement parameters
ng <- length(G_adj@x)

#read in locations of tortoises
tort_locs <- read.csv("watershed_assignments.csv")
#remove first tortoise
tort_locs <- tort_locs[-1,]
#subset out the ones we care about (columns 2-5 of tort_locs contains locations at resolutions 2,4,6,8)
#tort_loc <- tort_locs[,c(1,WBD_res/2+1)] 
#use handpicked?
tort_loc <- tort_locs[,c(1,7)] 

#read in pairwise genetic distances
pairwise_pis <- read.csv("all_angsd_snps.pwp.csv")
#remove first tortoise
pairwise_pis[which(pairwise_pis$etort1=="etort-1"),] <- NA

#remove self comparisons
pairwise_pis$pi[which(pairwise_pis$etort1==pairwise_pis$etort2)] <- NA

#find names of individuals
tort_names <- unique(tort_loc[,1])

#find number of individuals
n_ind <- length(tort_names)

#find number of locations (columns 2-5 of tort_locs contains locations at resolutions 2,4,6,8)
n <- length(unique(tort_loc[,2]))

#find number of torts in each location
n_ind_loc <- vector(length=n)
for(i in 1:n){
  n_ind_loc[i] <- length(which(tort_loc[,2]==i))
}

#create full (individual to individual) genetic distance matrix from pairwise_pis
dist <- matrix(nrow=n_ind,ncol=n_ind)
for(i in 1:n_ind){
  #fill in below diagonal
  dist[i:n_ind,i] <- pairwise_pis$pi[which(pairwise_pis$etort1==tort_names[i])]
}
#fill in upper triangle
dist[upper.tri(dist)] <- t(dist)[upper.tri(dist)]

#initialize matrices and helping variables
Hs <- matrix(nrow=n,ncol=n) #mean genetic distance matrix (location to location)
Hs_se <- matrix(nrow=n,ncol=n) #matrix for standard error of the mean for values of Hs
i_inds <- list() #list for indices of individuals from 1st location
j_inds <- list() #list for indices of individuals from 2nd location
#inds <- list() #list for sorting i and j indices
dist_ij <- list() #list for subsets of full (genome to genome) genetic distance
k <- 0

#fill in lower triangle of Hs and Hs_se
for(i in 1:n){
  for(j in 1:i){ 
    k <- k+1
    #find indices (and distances) of both genomes of the individuals in the ith and jth locations
    i_inds[[k]] <- which(tort_loc[,2]==i)
    j_inds[[k]] <- which(tort_loc[,2]==j)
    #concatenate columns and sort each row independently
    #inds[[k]] <- t(apply(cbind(i_inds[[i]],j_inds[[j]]),1,sort))
    #subset of pairwise pis for individuals in i vs j
    #dist_ij[[k]] <- pairwise_pis$pi[intersect(which(pairwise_pis$etort1==tort_loc[inds[[k]][,1],1]),
    #                                          which(pairwise_pis$etort2==tort_loc[inds[[k]][,2],1]))]
    #intersect(which(pairwise_pis$etort1 %in% tort_loc[i_inds[[1]],1]),
    #          which(pairwise_pis$etort2 %in% tort_loc[j_inds[[1]],1]))
    #subset matrix
    dist_ij[[k]] <- dist[i_inds[[k]],j_inds[[k]]]
    #remove data on and above diag if this is a self comparison
    if(i==j) dist_ij[[k]][upper.tri(dist_ij[[k]],diag=TRUE)] <- NA
    #means for each pair
    Hs[i,j] <- mean(dist_ij[[k]],na.rm=TRUE)
    #standard error of the mean (#divide by sqrt of min # of ind in either location)
    Hs_se[i,j] <- sd(dist_ij[[k]],na.rm=TRUE)/sqrt(min(n_ind_loc[i],n_ind_loc[j]))
    #print(length(dist_ij[[k]]))
    #print(length(locsi),locsj)
  }
}

#make it symmetric (replace upper triangle  with lower triangle)
Hs[upper.tri(Hs)] <- t(Hs)[upper.tri(Hs)]
Hs_se[upper.tri(Hs_se)] <- t(Hs_se)[upper.tri(Hs_se)]

#make H into vector
hs <- as.vector(Hs[upper.tri(Hs,diag=TRUE)])
#make standard error into a vector
hs_se <- as.vector(Hs_se[upper.tri(Hs_se,diag=TRUE)])

#number of iterations for mcmc
preburn_iter=1e2
burn_iter=3e2
iter=8e4

seed <- sample(1000000,1)

#save workspace
save(list=ls(),file=paste0("tort_mcmc_hand.RData"))

a <- run.mcmc(G_adj_known=TRUE,G_adj=G_adj,g_known=FALSE,const_coal=FALSE,H=Hs,h_se=hs_se,seed=seed,
              preburn_iter=preburn_iter,burn_iter=burn_iter,iter=iter,noisy_H=FALSE,type="coal")

#save workspace again
save(list=ls(),file=paste0("tort_mcmc_hand.RData"))

#boxplot(a$ans$g[1000*(1:(iter/1000)),],outline=FALSE,main=paste0("Posterior Distributions: "),
#        names=paste0("g",1:ng),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)

#boxplot(a$ans$gam[1000*(1:(iter/1000)),],outline=FALSE,main=paste0("Posterior Distributions: "),
#        names=paste0("gam",1:n),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)

#plot(a$ans$lllh[100*(1:(iter/100))],type='l')

#for(i in 1:ng){
#  plot(a$ans$g[1000*(1:(iter/1000)),i],type='l')
#}

#for(i in 1:n){
#  plot(a$ans$gam[1000*(1:(iter/1000)),i],type='l')
#}

#h_fit <- H_and_lllh(a$ans$g_med,a$ans$gam_med,G_adj,n,hs,hs_se,1,1,const_coal = FALSE)
#plot(hs,h_fit$hcalc/a$scale)
#abline(0,1)

#h_fit <- H_and_lllh(a$ans$g[which.max(a$ans$lllh),],a$ans$gam[which.max(a$ans$lllh),],G_adj,n,hs,hs_se,1,1,const_coal = FALSE)
#plot(hs,h_fit$hcalc/a$scale)
#abline(0,1)

#gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
#igraph::E(gr)$curved <- TRUE
#igraph::V(gr)$color <- rainbow(13)
#plot(gr,layout=gCentroid(sub_wbd8, byid=TRUE)@coords,edge.width=a$ans$g_med,
#     edge.label=paste0("g",1:ng,"=",round(a$ans$g_med*1000)/1000),main=paste0("Graph Structure"))

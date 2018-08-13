#do inference on a 5x3 with barriers grid using data from slim

require(gene.flow.inference)
filepath <- getwd()

#read individual positions into table
positions <- read.table("positions3.txt",header=TRUE)
#positions$x <- 5-positions$x #accidentally made is backwards on early run

#read all genomes
genomes <- t(as.matrix(do.call(cbind,strsplit(scan("ms3.txt",skip=3,what='numeric'),""))))

#specify dimensions of grid
width <- 5
height <- 3

#number of states in the system
n <- width*height

#find adjacency matrix
G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height)))
G <- G_adj

#find number of elements in G
ng <- length(G_adj@x)

#number of individuals per location
ni <- 50
#initialize vector for ids of individuals
indiv <- vector(length=ni*n)

#specify area bounds
xmin <- 0
xmax <- 5
ymin <- 0
ymax <- 3
#find size of each location square
loc_w <- (xmax-xmin)/width #width of each location square
loc_h <- (ymax-ymin)/height #height of each location square

#initialize matrix for centers of locations
centers <- matrix(nrow=n,ncol=2)

#initialize list for indices of individuals near center of each location
indiv <- list()
#initialize vector for indices of sampled individuals (from all locations)
indiv_s <- rep(NA,ni*n)

#loop through each location left to right, bottom to top
#sample ni individuals near the center of each location
for(j in 1:height){
  for(i in 1:width){
    #index for locations
    k <- (j-1)*width+i
    #find center of each location
    centers[k,] <- c(i*loc_w-loc_w/2,j*loc_h-loc_h/2)
    #find individuals in the middle quarter (middle half in each dimension) #is there a better way to do this?
    indiv[[k]] <- intersect(intersect(which(positions$x<centers[k,1]+loc_w/4),which(positions$x>centers[k,1]-loc_w/4)), 
                            intersect(which(positions$y<centers[k,2]+loc_h/4),which(positions$y>centers[k,2]-loc_h/4)))
    if(length(indiv[[k]]) > ni){ #if there are at least as many individuals as we want from each location
      indiv_s[((k-1)*ni+1):(k*ni)] <- indiv[[k]][1:ni] #same ni of the total individuals (order is already randomized)
    }
    else{
      indiv_s[((k-1)*ni+1):((k-1)*ni+1+length(indiv[[k]]))] <- indiv[[k]] #if there are fewer than ni, just take all of them
      print(paste(c("Only",length[[k]],"individuals in square",k)))
    }
    #indiv[(j*(i-1)+i):(j*(i-1)+i+(ni-1))] <- order(sqrt((positions$x-(i-.5))^2+(positions$y-(j-.5))^2))[1:ni]
  }
}

#find indices of genomes of sampled individuals (2 genomes per individuals since it is diploid)
gen_ind <- as.vector(rbind((2*indiv_s-1),2*indiv_s)) #ex: if indiv_s <- c(2,3,5), this returns c(3,4,5,6,9,10)
#subset out genomes from chosen individuals
genomes_s <- genomes[gen_ind,]

#may want to delete unselected genomes to free memory
rm(genomes)
gc()

#calculate distance matrix
#find genome-genome distance
gen_dist <- dist(genomes_s,method="manhattan",diag=FALSE) #may take a while
#convert to matrix
dist <- as.matrix(gen_dist)

#initialize matrices and helping variables
Hs <- matrix(nrow=n,ncol=n) #mean genetic distance matrix (location to location)
Hs_se <- matrix(nrow=n,ncol=n) #matrix for standard error of the mean for values of Hs
dist_ij <- list() #list for subsets of full (genome to genome) genetic distance
k <- 0

#fill in lower triangle of Hs and Hs_se
for(i in 1:n){ #i, j, and k are different here as compared to above
  for(j in 1:i){ 
    k <- k+1
    #find indices (and distances) of both genomes of the individuals in the ith and jth locations
    #subset of matrix for each pair of locations
    dist_ij[[k]] <- dist[((i-1)*2*ni+1):(i*2*ni),((j-1)*2*ni+1):(j*2*ni)] #ex: if ni=50,i=1,j=3, then subsets dist[[k]]=dist[1:100,201:300]
    #remove data on and above diag if this is a self comparison
    if(i==j) dist_ij[[k]][upper.tri(dist_ij[[k]],diag=TRUE)] <- NA
    #means for each pair
    Hs[i,j] <- mean(dist_ij[[k]],na.rm=TRUE)
    #standard error of the mean (#divide by sqrt of min # of ind in either location)
    Hs_se[i,j] <- sd(dist_ij[[k]],na.rm=TRUE)/sqrt(min(length(indiv[[i]]),length(indiv[[j]]),ni))
    #print(length(dist_ij[[k]]))
    #print(length(locsi),locsj)
  }
}
#make it symmetric (replace upper triangle  with lower triangle)
Hs[upper.tri(Hs)] <- t(Hs)[upper.tri(Hs)]
Hs_se[upper.tri(Hs_se)] <- t(Hs_se)[upper.tri(Hs_se)]
Hs
cairo_pdf(filename=paste0(filepath,"/H_matrix",".pdf"),width=5,height=5)
  image(Hs,xaxt='n',yaxt='n')
dev.off()
#make H into vector
hs <- as.vector(Hs[upper.tri(Hs,diag=TRUE)])
#make standard error into a vector
hs_se <- as.vector(Hs_se[upper.tri(Hs_se,diag=TRUE)])

#do some plotting

#view histograms
hist(positions$x)
hist(positions$y)

#2d histogram (full)
#k <- MASS::kde2d(positions$x,positions$y,n=c(50,30))
#cairo_pdf(filename=paste0(filepath,"/pop_heatmap_full.pdf"),width=5,height=3)
#  image(k)
#dev.off()
#2d histogram (subset)
#ks <- MASS::kde2d(positions[start+1:ind,2],positions[start+1:ind,3],n=c(50,30))
#image(ks)
#png(filename=paste0(filepath,"/ind_locs_sub",start,".png"),width=1200,height=800)
#  plot(positions[start+1:ind,2],positions[start+1:ind,3])
#dev.off()

#make landscape (current way)
mapValues = matrix(rep(1.0,1500),nrow=30,ncol=50);
mapValues[20+1,10:29+1] = 0;
mapValues[9+1,20:39+1] = 0;
mapValues = t(mapValues[,(ncol(mapValues)-1):0]);


cairo_pdf(filename=paste0(filepath,"/ind_locs_5x3b.pdf"),width=7,height=5)
  image(0:49/(49/5),0:29/(29/3),mapValues,xlim=c(0,5),ylim=c(0,3),
        xlab="Position (x)",ylab="Position (y)") #plot landscape
  points(positions[,2:3],pch=46) #all individuals (small points)
  for(k in 1:n){
    points(positions[indiv[[k]][1:ni],2:3],pch=20) #larger points
  }
  #plot boxes
  abline(h=1); abline(h=2); abline(v=1); abline(v=2); abline(v=3); abline(v=4)
dev.off()


#find distances between all individuals (repeat each position twice because they are diploid)
pos_dist <- dist(positions[rep(indiv_s,each=2),2:3],method="euclidian",diag=FALSE)

png(filename=paste0(filepath,"/ibd_5x3b",".png"),width=1200,height=800)
  plot(as.vector(pos_dist),as.vector(gen_dist),pch=16,col=adjustcolor('black', 0.02),cex.lab=1.2,cex.main=2,
       main="Genetic Distance vs Geographic Distance",xlab="Geographic Distance",ylab="Genetic Distance")
  abline(lsfit(as.vector(pos_dist),as.vector(gen_dist)),col="blue")
dev.off()
  
#use gene flow inference method
preburn_iter_coal <- 1e6
burn_iter_coal <- 3e6
iter_coal <- 4e6

#set seed
seed <- sample(1000000000,1)
set.seed(seed)

save.image(paste0(filepath,"/workspace.RData"))

system.time(a5x3b <- run.mcmc(width,height,g_known=FALSE,const_coal=TRUE,H=Hs,h_se=hs_se,seed=seed,
                                  preburn_iter=preburn_iter_coal,burn_iter=burn_iter_coal,iter=iter_coal,noisy_H=FALSE,type="coal"))

save.image(paste0(filepath,"/workspace.RData"))

boxplot(a5x3b$ans$g[10*(1:(iter_coal/10)),order(a5x3b$ans$g_med)],outline=FALSE,main="5x3 Graph with Barrier (Coalescence Time Inference)",
        names=paste0("g",order(a5x3b$ans$g_med)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)

cairo_pdf(filename=paste0(filepath,"/5x3grid",".pdf"),width=8,height=7)
  par(mar=(c(0,0,2,0))) #margins all zero except top
  gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
  igraph::E(gr)$curved <- TRUE
  plot(gr,layout=as.matrix(expand.grid(1:width,1:height)),edge.width=2*a5x3b$ans$g_med,
       edge.label=paste0("g",1:ng,"=",round(a5x3b$ans$g_med*1000)/1000),main="5x3 Graph with Barriers")
dev.off()
  
#plot trace
index <- 5
plot(c(a5x3b$pburn$g[,index],a5x3b$burn$g[,index],a5x3b$ans$g[,index]),type='l')


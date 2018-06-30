#function for doing inference on a slim output (or any ms output with continuout space position data)

slim.inference <- function(name="",fname="",xlim=c(0,1),ylim=c(0,1),width=4,height=4,ni=50,sample_area=.5,
                           seed=sample(1000000000,1),preburn_iter=1e6,burn_iter=3e6,iter=4e6,inpath=getwd(),outpath=getwd(),
                           pos_name="positions3.txt",gen_name="ms3.txt",landscape_matrix=NULL){
  #name should be readable for plots
  #fname should not have any spaces
  #ni is the number of individuals sampled per location
  #sample_area is the fraction of the square in each dimension from which you sample individuals
  #   #ex. if sample_area is .5 and the sqaure limits are from 0 to 0.5 in x and 0.5 to 1.0 in y,
  #   #then sample individuals from .125 to .375 in x and .625 and .875 in y
  
  #set rng seed
  set.seed(seed)
  
  #read individual positions into table
  positions <- read.table(paste(inpath,pos_name,sep="/"),header=TRUE)
  
  #find total number of individuals
  t_ind <- length(positions[,1])
  
  #find number of elements in each genome
  gen_length <- length(do.call(cbind,strsplit(scan(paste(inpath,gen_name,sep="/"),skip=3,what='numeric',nlines=1),"")))
  
  #declare list for different sections of the file
  genomes <- list()
  
  #find cutoffs
  cutoffs <- round((0:5)*2*t_ind/5) #2 genomes per individual
  
  #read in values
  for(i in 1:5){
    genomes[[i]] <- 
      Matrix::t(Matrix::Matrix(
        as.numeric(do.call(cbind,strsplit(
          scan(paste(inpath,gen_name,sep="/"),skip=3+cutoffs[i],what='numeric',nlines=(cutoffs[i+1]-cutoffs[i])),""))),nrow=gen_length))
    gc() #clean up
  }
  #combine matrices
  genomes <- do.call(Matrix::rBind,genomes)
  
  #number of states in the system
  n <- width*height
  
  #find adjacency matrix
  G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height)))
  G <- G_adj
  
  #find number of elements in G
  ng <- length(G_adj@x)
  
  #initialize vector for ids of individuals
  indiv <- vector(length=ni*n)
  
  #specify area bounds
  xmin <- xlim[1]
  xmax <- xlim[2]
  ymin <- ylim[1]
  ymax <- ylim[2]
  #find size of each location square
  loc_w <- (xmax-xmin)/width #width of each location square
  loc_h <- (ymax-ymin)/height #height of each location square
  
  #initialize matrix for centers of locations
  centers <- matrix(nrow=n,ncol=2)
  
  #initialize list for indices of individuals near center of each location
  indiv <- list()
  #initialize list for indices of sampled individuals (from all locations)
  indiv_s_l <- list()
  
  #loop through each location left to right, bottom to top
  #sample ni individuals near the center of each location
  for(j in 1:height){
    for(i in 1:width){
      #index for locations
      k <- (j-1)*width+i
      #find center of each location
      centers[k,] <- c(i*loc_w-loc_w/2,j*loc_h-loc_h/2)
      #find individuals in the middle quarter (middle half in each dimension) #is there a better way to do this?
      indiv[[k]] <- intersect(intersect(which(positions$x<centers[k,1]+loc_w/(2/sample_area)),
                                        which(positions$x>centers[k,1]-loc_w/(2/sample_area))), 
                              intersect(which(positions$y<centers[k,2]+loc_h/(2/sample_area)),
                                        which(positions$y>centers[k,2]-loc_h/(2/sample_area))))
      if(length(indiv[[k]]) > ni){ #if there are at least as many individuals as we want from each location
        #indiv_s[((k-1)*ni+1):(k*ni)] <- indiv[[k]][1:ni] #same ni of the total individuals (order is already randomized)
        indiv_s_l[[k]] <- indiv[[k]][1:ni] #same ni of the total individuals (order is already randomized)
      }
      else{
        #indiv_s[((k-1)*ni+1):((k-1)*ni+1+length(indiv[[k]]))] <- indiv[[k]] #if there are fewer than ni, just take all of them
        indiv_s_l[[k]] <- indiv[[k]]
        print(paste(c(length(indiv[[k]]),"individuals in location",k)))
      }
    }
  }
  
  #unlist indices of sampled individual
  indiv_s <- unlist(indiv_s_l)
  
  #find starting indices for which genomes correspond to each location
  start_ind <- rep(0,(n+1)) #declare vector
  start_ind[1] <- 1 #first location starts at beginning
  for(k in 2:(n+1)){ #also include "start" point for after the end of the last location for later convenience
    start_ind[k] <- start_ind[k-1]+2*length(indiv_s_l[[k-1]]) #multiply by 2 since it is diploid
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
      #dist_ij[[k]] <- dist[((i-1)*2*ni+1):(i*2*ni),((j-1)*2*ni+1):(j*2*ni)] #ex: if ni=50,i=1,j=3, then subsets dist[[k]]=dist[1:100,201:300]
      dist_ij[[k]] <- dist[start_ind[i]:(start_ind[i+1]-1),start_ind[j]:(start_ind[j+1]-1)]
      #remove data on and above diag if this is a self comparison
      if(i==j) dist_ij[[k]][upper.tri(dist_ij[[k]],diag=TRUE)] <- NA
      #means for each pair
      Hs[i,j] <- mean(dist_ij[[k]],na.rm=TRUE)
      #standard error of the mean (#divide by sqrt of min # of ind in either location)
      #Hs_se[i,j] <- sd(dist_ij[[k]],na.rm=TRUE)/sqrt(min(length(indiv[[i]]),length(indiv[[j]]),ni))
      Hs_se[i,j] <- sd(dist_ij[[k]],na.rm=TRUE)/sqrt(min(length(indiv_s_l[[i]]),length(indiv_s_l[[j]])))
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
  
  #do some plotting
  
  #plot H matrix
  cairo_pdf(filename=paste0(outpath,"/H_matrix_",fname,".pdf"),width=6,height=6)
    image(Hs,xaxt='n',yaxt='n',main=paste0("Genetic Distance Matrix (H): ",name))
  dev.off()
  
  #plot individuals on landscape
  cairo_pdf(filename=paste0(outpath,"/ind_locs_",fname,".pdf"),width=7,height=(5*(height/width)+2))
  #make plot title
  main <- paste0("Individual Locations: ",name)
    if(is.matrix(landscape_matrix) == TRUE){
      image(seq(xmin,xmax,length.out=nrow(landscape_matrix)),seq(ymin,ymax,length.out=ncol(landscape_matrix)),
            mapValues,xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xlab="Position (x)",ylab="Position (y)",main=main) #plot landscape
      points(positions[,2:3],pch=46) #all individuals (small points)
    } else{
      plot(positions[,2:3],pch=46,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xaxs='i',yaxs='i', #make plot end exactly at limits
           xlab="Position (x)",ylab="Position (y)",main=main) #all individuals (small points)
    }
    for(k in 1:n){
      points(positions[indiv[[k]][1:ni],2:3],pch=20) #larger points for sampled individuals
    }
    #plot grid lines
    for(i in 0:(width)) abline(v=(xmin+i*loc_w))
    for(i in 0:(height)) abline(h=(ymin+i*loc_h))
  dev.off()
  
  #find distances between all individuals (repeat each position twice because they are diploid)
  pos_dist <- dist(positions[rep(indiv_s,each=2),2:3],method="euclidian",diag=FALSE)
  
  png(filename=paste0(outpath,"/ibd_",fname,".png"),width=1200,height=800)
    plot(as.vector(pos_dist),as.vector(gen_dist),pch=16,col=adjustcolor('black', 0.02),cex.lab=1.2,cex.main=2,
         main=paste0("Genetic Distance vs Geographic Distance: ",name),xlab="Geographic Distance",ylab="Genetic Distance")
    abline(lsfit(as.vector(pos_dist),as.vector(gen_dist)),col="blue")
  dev.off()
  
  #save workspace
  #save(ls(),list=ls(),file=paste0(outpath,"/",fname,".RData"))
  save(list=ls(),file=paste0(outpath,"/",fname,".RData"))
  
  #run mcmc
  system.time(a <- run.mcmc(width,height,g_known=FALSE,const_coal=TRUE,H=Hs,h_se=hs_se,seed=seed,
                               preburn_iter=preburn_iter,burn_iter=burn_iter,iter=iter,noisy_H=FALSE,type="coal"))
  
  save(list=ls(),file=paste0(outpath,"/",fname,".RData"))
  
  cairo_pdf(filename=paste0(outpath,"/posterior_dists_",fname,".pdf"),width=8,height=7)
    boxplot(a$ans$g[10*(1:(iter/10)),order(a$ans$g_med)],outline=FALSE,main=paste0("Posterior Distributions: ",fname),
            names=paste0("g",order(a$ans$g_med)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  dev.off()
    
  cairo_pdf(filename=paste0(outpath,"/grid_",fname,".pdf"),width=8,height=7)
    par(mar=(c(0,0,2,0))) #margins all zero except top
    gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
    igraph::E(gr)$curved <- TRUE
    plot(gr,layout=as.matrix(expand.grid(1:width,1:height)),edge.width=2*a$ans$g_med,
         edge.label=paste0("g",1:ng,"=",round(a$ans$g_med*1000)/1000),main=paste0("Graph Structure: ",name))
  dev.off()
  
}
#functions for plotting data and inference that has already been completed

plot.ind.locs <- function(name,positions,indiv_s,width,heigth,loc_w,loc_h,landscape_matrix=NA,xlim,ylim){
  #make plot title
  main <- paste0("Individual Locations: ",name)
  if(is.matrix(landscape_matrix) == TRUE){
    image(seq(xlim[1],xlim[2],length.out=nrow(landscape_matrix)),seq(ylim[1],ylim[2],length.out=ncol(landscape_matrix)),
          landscape_matrix,xlim=xlim,ylim=ylim,
          xlab="Position (x)",ylab="Position (y)",main=main) #plot landscape
    points(positions[,2:3],pch=46) #all individuals (small points)
  } else{
    plot(positions[,2:3],pch=46,xlim=xlim,ylim=ylim,xaxs='i',yaxs='i', #make plot end exactly at limits
         xlab="Position (x)",ylab="Position (y)",main=main) #all individuals (small points)
  }
  points(positions[indiv_s,2:3],pch=20) #larger points for sampled individuals
  #plot grid lines
  for(i in 0:(width)) abline(v=(xmin+i*loc_w))
  for(i in 0:(height)) abline(h=(ymin+i*loc_h))
}

plot.H.matrix <- function(Hs,name){
  image(Hs,xaxt='n',yaxt='n',main=paste0("Genetic Distance Matrix (H): ",name))
}

plot.ibd <- function(pos_dist,gen_dist,name){
  plot(as.vector(pos_dist),as.vector(gen_dist),pch=16,col=adjustcolor('black', 0.02),cex.lab=1.2,cex.main=2,
       main=paste0("Genetic Distance vs Geographic Distance: ",name),xlab="Geographic Distance",ylab="Genetic Distance")
  #abline(lsfit(as.vector(pos_dist),as.vector(gen_dist)),col="blue")
}

plot.posteriors <- function(g_ts,g_med,name,ord=1:length(g_med),colors=1,subdiv=100){
  boxplot(g_ts[subdiv*(1:(length(g_ts[,1])/subdiv)),ord],outline=FALSE,boxcol=colors,medcol=colors,whiskcol=colors,staplecol=colors,
          main=paste0("Posterior Distributions: ",name),
          names=paste0("g",ord),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
}

plot.grid <- function(g_med,width,height,G_adj_known=FALSE,G_adj=NA,ng=length(G_adj@x),name,mar=c(0,0,2,0)){
  par(mar=mar) #margins all zero except top default
  if(G_adj_known==FALSE && width%%1==0 && height%%1==0){ #make sure width and height are integers
    #declare generating matrix
    G_adj <- igraph::get.adjacency(igraph::make_lattice(c(width,height))) 
    ng <- length(G_adj@x)
  }
  gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
  igraph::E(gr)$curved <- TRUE
  plot(gr,layout=as.matrix(expand.grid(1:width,1:height)),edge.width=2*g_med,
       edge.label=paste0("g",1:ng,"=",round(g_med*1000)/1000),main=paste0("Graph Structure: ",name))
}
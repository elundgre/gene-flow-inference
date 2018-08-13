#4x4 coalescence vs commute plots

a <- list(a4x4u,a4x4s,a4x4as,a4x4asb)

layout(matrix(1:4,4))
par(mar=c(2.5, 2.7, 1, 0.5), mgp=c(1.4, 0.5,0))

plot(a[[1]]$H,a[[1]]$Hc,main="Commute Time + Diversity vs Coalescence Time (Uniform Graph)",
     xlab="Coalescence Time",ylab="Commute Time + Diversity")
abline(0,1) #add y=x line

plot(a[[2]]$H,a[[2]]$Hc,main="Commute Time + Diversity vs Coalescence Time (Symmetric Graph)",
     xlab="Coalescence Time",ylab="Commute Time + Diversity")
abline(0,1) #add y=x line

plot(a[[3]]$H,a[[3]]$Hc,main="Commute Time + Diversity vs Coalescence Time (Asymmetric Graph)",
     xlab="Coalescence Time",ylab="Commute Time + Diversity")
abline(0,1) #add y=x line

plot(a[[4]]$H,a[[4]]$Hc,main="Commute Time + Diversity vs Coalescence Time (Biased Asymmetric Graph)",
     xlab="Coalescence Time",ylab="Commute Time + Diversity")
abline(0,1) #add y=x line

layout(matrix(1:8,4,byrow=TRUE))
par(mar=c(3.2, 3, 2, 1), mgp=c(2, 0.6,0))

#uniform graph
boxplot(a[[1]]$com$ans$g[100*(1:(a[[1]]$iter/100)),order(a[[1]]$g)],outline=FALSE,main="4x4 Uniform Graph (Commute Time Inference)",
        names=paste0("g",order(a[[1]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,6),cex.axis=0.6)
points(a[[1]]$g[order(a[[1]]$g)],pch=19,col=3)

boxplot(a[[1]]$coal$ans$g[100*(1:(a[[1]]$iter/100)),order(a[[1]]$g)],outline=FALSE,main="4x4 Uniform Graph (Coalescence Time Inference)",
        names=paste0("g",order(a[[1]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,6),cex.axis=0.6)
points(a[[1]]$g[order(a[[1]]$g)],pch=19,col=3)

#symmetric graph
boxplot(a[[2]]$com$ans$g[100*(1:(a[[2]]$iter/100)),order(a[[2]]$g)],outline=FALSE,main="4x4 Symmetric Graph (Commute Time Inference)",
        names=paste0("g",order(a[[2]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,9),cex.axis=0.6)
points(a[[2]]$g[order(a[[2]]$g)],pch=19,col=3)

boxplot(a[[2]]$coal$ans$g[100*(1:(a[[2]]$iter/100)),order(a[[2]]$g)],outline=FALSE,main="4x4 Symmetric Graph (Coalescence Time Inference)",
        names=paste0("g",order(a[[2]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,9),cex.axis=0.6)
points(a[[2]]$g[order(a[[2]]$g)],pch=19,col=3)

#asymmetric graph
boxplot(a[[3]]$com$ans$g[100*(1:(a[[3]]$iter/100)),order(a[[3]]$g)],outline=FALSE,main="4x4 Asymmetric Graph (Commute Time Inference)",
        names=paste0("g",order(a[[3]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,6),cex.axis=0.6)
points(a[[3]]$g[order(a[[3]]$g)],pch=19,col=3)

boxplot(a[[3]]$coal$ans$g[100*(1:(a[[3]]$iter/100)),order(a[[3]]$g)],outline=FALSE,main="4x4 Asymmetric Graph (Coalescence Time Inference)",
        names=paste0("g",order(a[[3]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,6),cex.axis=0.6)
points(a[[3]]$g[order(a[[3]]$g)],pch=19,col=3)

#biased asymmetric graph
boxplot(a[[4]]$com$ans$g[100*(1:(a[[4]]$iter/100)),order(a[[4]]$g)],outline=FALSE,main="4x4 Biased Asymmetric Graph (Commute Time Inference)",
        names=paste0("g",order(a[[4]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,16),cex.axis=0.6)
points(a[[4]]$g[order(a[[4]]$g)],pch=19,col=3)

boxplot(a[[4]]$coal$ans$g[100*(1:(a[[4]]$iter/100)),order(a[[4]]$g)],outline=FALSE,main="4x4 Biased Asymmetric Graph (Coalescence Time Inference)",
        names=paste0("g",order(a[[4]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2,ylim=c(0,16),cex.axis=0.6)
points(a[[4]]$g[order(a[[4]]$g)],pch=19,col=3)


#plot grid structures

layout(matrix(1:4,2))
par(mar=c(2,3,2,3))

gr <- igraph::graph_from_adjacency_matrix(a[[1]]$G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
plot(gr,layout=as.matrix(expand.grid(1:a[[1]]$width,1:a[[1]]$height)),edge.width=3*a[[1]]$g,edge.label=paste0("g",1:length(a[[1]]$g),"=",round(a[[1]]$g*1000)/1000),
     main="4x4 Uniform Graph",asp=0)

gr <- igraph::graph_from_adjacency_matrix(a[[2]]$G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
plot(gr,layout=as.matrix(expand.grid(1:a[[2]]$width,1:a[[2]]$height)),edge.width=3*a[[2]]$g,edge.label=paste0("g",1:length(a[[2]]$g),"=",round(a[[2]]$g*1000)/1000),
     main="4x4 Symmetric Graph",asp=0)

gr <- igraph::graph_from_adjacency_matrix(a[[3]]$G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
plot(gr,layout=as.matrix(expand.grid(1:a[[3]]$width,1:a[[3]]$height)),edge.width=3*a[[3]]$g,edge.label=paste0("g",1:length(a[[3]]$g),"=",round(a[[3]]$g*1000)/1000),
     main="4x4 Asymmetric Graph",asp=0)

gr <- igraph::graph_from_adjacency_matrix(a[[4]]$G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
plot(gr,layout=as.matrix(expand.grid(1:a[[4]]$width,1:a[[4]]$height)),edge.width=3*a[[4]]$g,edge.label=paste0("g",1:length(a[[4]]$g),"=",round(a[[4]]$g*1000)/1000),
     main="4x4 Biased Asymmetric Graph",asp=0)


#find errors
require(matrixStats)

e4x4u_com <- mean(abs(a[[1]]$g-colMedians(a[[1]]$com$ans$g)))
print(e4x4u_com)

e4x4u_coal <- mean(abs(a[[1]]$g-colMedians(a[[1]]$coal$ans$g)))
print(e4x4u_coal)

e4x4s_com <- mean(abs(a[[2]]$g-colMedians(a[[2]]$com$ans$g)))
print(e4x4s_com)

e4x4s_coal <- mean(abs(a[[2]]$g-colMedians(a[[2]]$coal$ans$g)))
print(e4x4s_coal)

e4x4as_com <- mean(abs(a[[3]]$g-colMedians(a[[3]]$com$ans$g)))
print(e4x4as_com)

e4x4as_coal <- mean(abs(a[[3]]$g-colMedians(a[[3]]$coal$ans$g)))
print(e4x4as_coal)

e4x4asb_com <- mean(abs(a[[4]]$g-colMedians(a[[4]]$com$ans$g)))
print(e4x4asb_com)

e4x4asb_coal <- mean(abs(a[[4]]$g-colMedians(a[[4]]$coal$ans$g)))
print(e4x4asb_coal)

#many runs of small coal vs com

require(parallel)
require(gene.flow.inference)

#file path for output
out_path <- "~/Dropbox/code/gene-flow-inference/plots/"
#or use current path
#outpath <- getwd()

#number of replicates
reps <- 1

values <- c(rep(1/1000,reps),rep(1/200,reps),rep(1/100,reps),rep(1/50,reps))
seeds <- sample(1000000000,reps)

time_m <- system.time(mult <- mcmapply(mult.small,noise=values,const_coal=TRUE,seed=rep(seeds,4),mc.cores=4,SIMPLIFY=FALSE))
save.image(paste0(out_path,"3x2_mult6.RData"))

#png(file=paste0(out_path,"3x2_plots/mult6box%02d.png"),width=8*144, height=4*144, res=144, pointsize=10)
for(i in 1:(4*reps)){
  boxplot(mult[[i]]$a3x2s_com$ans$g[10*(1:(mult[[i]]$a3x2s_com$iter/10)),order(mult[[i]]$g)],outline=FALSE,
          main=paste("3x2 Symmetric Graph (Commute Time Inference)",i),
          names=paste0("g",order(mult[[i]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  points(mult[[i]]$g[order(mult[[i]]$g)],pch=19,col=3)

  boxplot(mult[[i]]$a3x2s_coal$ans$g[10*(1:(mult[[i]]$a3x2s_coal$iter/10)),order(mult[[i]]$g)],outline=FALSE,
          main=paste("3x2 Symmetric Graph (Coalescence Time Inference)",i),
          names=paste0("g",order(mult[[i]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  points(mult[[i]]$g[order(mult[[i]]$g)],pch=19,col=3)
  
  boxplot(mult[[i]]$a3x2s_coal$ans$gam[10*(1:(mult[[i]]$a3x2s_coal$iter/10)),1],outline=FALSE,
          main=paste("3x2 Symmetric Graph (Coalescencec Time Inference) gamma",i),
          names=paste0("gam",1),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  points(rep(1,n),pch=19,col=3)
  
  boxplot(mult[[i]]$a3x2s_com$ans$q[10*(1:(mult[[i]]$a3x2s_com$iter/10)),order(diag(mult[[i]]$a3x2s_com$H))],outline=FALSE,
          main=paste("3x2 Symmetric Graph (Commute Time Inference) q",i),
          names=paste0("q",order(diag(mult[[i]]$a3x2s_com$H))),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  points(diag(mult[[i]]$a3x2s_com$H)[order(diag(mult[[i]]$a3x2s_com$H))],pch=19,col=3)
}
#dev.off()

#plot log likelihoods
png(file=paste0(out_path,"3x2_plots/mult5lllh%02d.png"),width=8*144, height=4*144, res=144, pointsize=10)
for(i in 1:100){
  plot(mult[[i]]$a3x2s_com$pburn$lllh,type='l',main=paste("com preburn lllh",i))
  plot(mult[[i]]$a3x2s_com$burn$lllh,type='l',main=paste("com burn lllh",i))
  plot(mult[[i]]$a3x2s_com$ans$lllh,type='l',main=paste("com lllh",i))
  plot(mult[[i]]$a3x2s_coal$pburn$lllh,type='l',main=paste("coal preburn lllh",i))
  plot(mult[[i]]$a3x2s_coal$burn$lllh,type='l',main=paste("coal burn lllh",i))
  plot(mult[[i]]$a3x2s_coal$ans$lllh,type='l',main=paste("coal lllh",i))
}
dev.off()

coal_diff <- rep(0,100)
com_diff <- rep(0,100)
hc_diff <- rep(0,100)

for(i in 1:100){
  coal_diff[i] <- mean(abs(mult[[i]]$g-colMedians(mult[[i]]$a3x2s_coal$ans$g)))
  com_diff[i] <- mean(abs(mult[[i]]$g-colMedians(mult[[i]]$a3x2s_com$ans$g)))
}

hist(coal_diff[1:25],n=25)
hist(com_diff)

boxplot(cbind(coal_diff[1:25],coal_diff[26:50],coal_diff[51:75],coal_diff[76:100],
              com_diff[1:25],com_diff[26:50],com_diff[51:75],com_diff[76:100]))

for(i in 1:100){
  hc_diff[i] <- mean(abs(mult[[i]]$a3x2s_com$H_infer-HcCalc(colMedians(mult[[i]]$a3x2s_com$ans$g),
                                                         colMedians(mult[[i]]$a3x2s_com$ans$q),G_adj,n)))
}

plot(hc_diff,com_diff)

#find average interquartile range width
coal_iqrs <- rep(0,100)
com_iqrs <- rep(0,100)
for(i in 1:100){
  coal_iqrs[i] <- mean(colIQRs(mult[[i]]$a3x2s_coal$ans$g))
  com_iqrs[i] <- mean(colIQRs(mult[[i]]$a3x2s_com$ans$g))
}

boxplot(cbind(coal_iqrs[1:25],coal_iqrs[26:50],coal_iqrs[51:75],coal_iqrs[76:100],
              com_iqrs[1:25],com_iqrs[26:50],com_iqrs[51:75],com_iqrs[76:100]))

plot(coal_diff[1:25],com_diff[1:25])
plot(coal_diff[26:50],com_diff[26:50])
plot(coal_diff[51:75],com_diff[51:75])
plot(coal_diff[76:100],com_diff[76:100])


#compare error to proposal distribution step size
sdG_all_coal <- rep(0,100)
for(i in 1:100){
  sdG_all_coal[i] <- mult[[i]]$a3x2s_coal$ans$sdG
}
plot(sdG_all_coal,coal_diff)

#try burn-in with larger sigma
time4 <- system.time(ans4 <- findG_MH4(mult[[15]]$a3x2s_coal$H_infer,G_adj,const_coal=TRUE,iter=20000,
                                       sig2ep=(1/50)*mean(mult[[15]]$a3x2s_coal$H_infer)))
time4 <- system.time(ans4b <- findG_MH4(mult[[15]]$a3x2s_coal$H_infer,G_adj,const_coal=TRUE,iter=80000,fixed_start=TRUE,
                                       g_init=ans4$g[20000,],gam_init=ans4$gam[20000,],sig2ep=(1/1000)*mean(mult[[15]]$a3x2s_coal$H_infer)))

plot(colMeans(abs(mult[[15]]$g-t(ans4$g))),type='l')
plot(ans4$lllh,type='l')
plot(colMeans(abs(mult[[15]]$g-t(ans4b$g))),type='l')
plot(ans4b$lllh,type='l')

#mcmapply version
seeds <- sample(1000000000,25)
mult6 <- mcmapply(mult_small,noise=values,const_coal=TRUE,seed=rep(seeds,4),mc.cores=6,SIMPLIFY=FALSE)
save.image("~/Documents/3x2_mult6.RData")

for(i in 1:2){
  boxplot(test[[i]]$a3x2s_com$ans$g[10*(1:(iter/10)),order(test[[i]]$g)],outline=FALSE,main=paste("3x2 Symmetric Graph (Commute Time Inference)",i),
          names=paste0("g",order(mult[[i]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  points(test[[i]]$g[order(test[[i]]$g)],pch=19,col=3)
  
  boxplot(test[[i]]$a3x2s_coal$ans$g[10*(1:(iter_coal/10)),order(test[[i]]$g)],outline=FALSE,
          main=paste("3x2 Symmetric Graph (Coalescence Time Inference)",i),
          names=paste0("g",order(test[[i]]$g)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
  points(test[[i]]$g[order(test[[i]]$g)],pch=19,col=3)
}

for(i in 2){
  plot(test[[i]]$a3x2s_com$pburn$lllh,type='l',main=paste("com preburn lllh",i))
  plot(test[[i]]$a3x2s_com$burn$lllh,type='l',main=paste("com burn lllh",i))
  plot(test[[i]]$a3x2s_com$ans$lllh,type='l',main=paste("com lllh",i))
  plot(test[[i]]$a3x2s_coal$pburn$lllh,type='l',main=paste("coal preburn lllh",i))
  plot(test[[i]]$a3x2s_coal$burn$lllh,type='l',main=paste("coal burn lllh",i))
  plot(test[[i]]$a3x2s_coal$ans$lllh,type='l',main=paste("coal lllh",i))
}

coal_diff_test <- rep(0,6)
com_diff_test <- rep(0,6)

for(i in 1:6){
  coal_diff_test[i] <- mean(abs(test[[i]]$g-colMedians(test[[i]]$a3x2s_coal$ans$g)))
  com_diff_test[i] <- mean(abs(test[[i]]$g-colMedians(test[[i]]$a3x2s_com$ans$g)))
}

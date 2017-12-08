#script for running many small sims with varying noise on symmetric grids

require(parallel)
require(gene.flow.inference)

#file path for output
#out_path <- "~/Dropbox/code/gene-flow-inference/plots/"
#or use current path
out_path <- getwd()

#number of replicates
reps <- 25

noise_values <- c(rep(1/1000,reps),rep(1/200,reps),rep(1/100,reps),rep(1/50,reps))
seeds <- sample(1000000000,reps)

time_m <- system.time(mult <- mcmapply(mult.small,noise=noise_values,gam=1,const_coal=TRUE,symmetric=TRUE,
                                       seed=rep(seeds,4),mc.cores=6,SIMPLIFY=FALSE))
save.image(paste0(out_path,"mult_noise_sym.RData"))

coal_diff <- rep(0,100)
com_diff <- rep(0,100)
hc_diff <- rep(0,100)

for(i in 1:100){
  coal_diff[i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
}

hist(coal_diff[1:25],n=25)
hist(com_diff)

boxplot(cbind(coal_diff[1:25],coal_diff[26:50],coal_diff[51:75],coal_diff[76:100],
              com_diff[1:25],com_diff[26:50],com_diff[51:75],com_diff[76:100]))

#find mean absolute value differences between H and Hc
for(i in 1:100){
  hc_diff[i] <- mean(abs(mult[[i]]$a3x2_com$H_infer-HcCalc(mult[[i]]$a3x2_com$ans$g_med,mult[[i]]$a3x2_com$ans$q_med,G_adj,mult[[i]]$n)))
}

plot(hc_diff,com_diff)

#find average interquartile range width
coal_iqrs <- rep(0,100)
com_iqrs <- rep(0,100)
for(i in 1:100){
  coal_iqrs[i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqrs[i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_iqrs[1:25],coal_iqrs[26:50],coal_iqrs[51:75],coal_iqrs[76:100],
              com_iqrs[1:25],com_iqrs[26:50],com_iqrs[51:75],com_iqrs[76:100]))

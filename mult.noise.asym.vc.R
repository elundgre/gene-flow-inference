#different noise levels, asymmetric, unknown if coalescence rates are same everywhere

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

time_m <- system.time(mult <- mcmapply(mult.small,noise=noise_values,gam=1,const_coal=FALSE,symmetric=FALSE,
                                       seed=rep(seeds,4),mc.cores=6,SIMPLIFY=FALSE))
save.image(paste0(out_path,"mult_noise_asym_vc.RData"))

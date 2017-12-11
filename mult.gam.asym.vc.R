#different coalescence rates, asymmetric, unknown if coalescence rates are same everywhere
require(parallel)
require(gene.flow.inference)

#file path for output
#out_path <- "~/Dropbox/code/gene-flow-inference/plots/"
#or use current path
out_path <- getwd()

#number of replicates
reps <- 25

gam_values <- c(rep(10,reps),rep(1,reps),rep(0.1,reps))
seeds <- sample(1000000000,reps)

time_m <- system.time(mult <- mcmapply(mult.small,noise=1/200,gam=gam_values,const_coal=FALSE,symmetric=FALSE,
                                       seed=rep(seeds,4),mc.cores=6,SIMPLIFY=FALSE))
save.image("mult_gam_asym_vc.RData")

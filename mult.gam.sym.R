#different coalescence rates, symmetric

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

time_m <- system.time(mult <- mcmapply(mult.small,noise=1/200,gam=gam_values,const_coal=TRUE,symmetric=TRUE,
                                       seed=rep(seeds,4),mc.cores=6,SIMPLIFY=FALSE))
save.image("mult_gam_sym.RData")

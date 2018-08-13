#script for running slim.inference
mylib <- "/home/cmb-00/plr/elundgre/R_packages/"
.libPaths(c(.libPaths(),mylib))

require(gene.flow.inference)

name <- "Biased Migration"
fname <- "bias_2x2"

slim.inference(name=name,fname=fname,xlim=c(0,1),ylim=c(0,1),width=2,height=2,ni=200,sample_area = .75,
               preburn_iter=1e5,burn_iter = 1e6,iter = 1e6,pos_name="positions3.txt",gen_name="ms3.txt")

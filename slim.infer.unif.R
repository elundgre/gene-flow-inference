#script for running slim.inference
mylib <- "/home/cmb-00/plr/elundgre/R_packages/"
.libPaths(c(.libPaths(),mylib))

require(gene.flow.inference)

name <- "Uniform Migration"
fname <- "unif_4x4_1_2"

slim.inference(name=name,fname=fname,xlim=c(0,1),ylim=c(0,1),width=4,height=4,ni=50,sample_area = .75,
               preburn_iter=1e6,burn_iter = 3e6,iter = 4e6,pos_name="positions3.txt",gen_name="ms3.txt",ind_shift=50)

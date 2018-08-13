#script for running slim.inference
mylib <- "/home/cmb-00/plr/elundgre/R_packages/"
.libPaths(c(.libPaths(),mylib))

require(gene.flow.inference)

name <- "Landscape with Barriers (Commute Inference)"
fname <- "5x3b_com"

mapValues = matrix(rep(1.0,300*500),nrow=300,ncol=500);
mapValues[95:105,201:400] = 0;
mapValues[196:206,101:300] = 0;
mapValues = t(mapValues[,(ncol(mapValues)):1]);

slim.inference(name=name,fname=fname,xlim=c(0,5),ylim=c(0,3),width=5,height=3,ni=50,sample_area = .75,
               preburn_iter=1e6,burn_iter = 3e6,iter = 4e6,pos_name="positions3.txt",gen_name="ms3.txt",
               landscape_matrix = mapValues,type="com")

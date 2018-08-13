#script for plotting slim stuff
#move to the directory with output from the slim.infer.unif.R (uniform square case)
require(gene.flow.inference)

#filepaths
inpath <- getwd()
outpath <- getwd()
#outpath <- "~/Documents/data/slim/new/plots/"

#load workspace

#specify name and file name
name <- "Uniform Migration 1"
fname <- "unif_1"

cairo_pdf(filename=paste0(outpath,"/ind_locs_",fname,".pdf"),width=7,height=(5*(height/width)+2))
  plot.ind.locs(name=name,positions=positions,indiv_s=indiv_s,width=width,heigth=height,
                loc_w=loc_w,loc_h=loc_h,xlim=xlim,ylim=ylim)
dev.off()

cairo_pdf(filename=paste0(outpath,"/H_matrix_",fname,".pdf"),width=6,height=6)
  plot.H.matrix(Hs,name)
dev.off()

png(filename=paste0(outpath,"/ibd_",fname,".png"),width=1200,height=800)
  plot.ibd(pos_dist=pos_dist,gen_dist=gen_dist,name=name)
dev.off()

#specify order
g_ur <- c(3,6,9,11,14,18,22,15,19,23,25,28,32,36,29,33,37,39,41,44,47,42,45,48)
g_dl <- c(1,4,7,2,5,8,10,12,16,20,13,17,21,24,26,30,34,27,31,35,38,40,43,46)
ord <- c(3,1,6,4,9,7,11,2,14,5,18,8,22,10,15,12,19,16,23,20,25,13,28,17,32,21,36,24,
         29,26,33,30,37,34,39,27,41,31,44,35,47,38,42,40,45,43,48,46)
#colors <- rep(c(rgb(0,0,0),rgb(0,1,0)),length(ord)/2)

cairo_pdf(filename=paste0(outpath,"/post_dists_",fname,".pdf"),width=8,height=7)
  plot.posteriors(g_ts=a$ans$g,g_med=a$ans$g_med,name=name)
dev.off()

cairo_pdf(filename=paste0(outpath,"/grid_",fname,".pdf"),width=8,height=7)
  plot.grid(g_med=a$ans$g_med,width=width,height=height,G_adj=G_adj,ng=ng,name=name)
dev.off()


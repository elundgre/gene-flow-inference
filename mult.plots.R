#script for making plots of all the 3x2 multi-runs

#declare lists for the mean absolute value errors of the posterior medians
coal_diff <- list()
com_diff <- list()
coal_iqr <- list()
com_iqr <- list()

#index variable for which sim we are doing
sim <- 1

#load noisy symmetric
load("mult_noise_sym.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[1]][1:25],coal_diff[[1]][26:50],coal_diff[[1]][51:75],coal_diff[[1]][76:100],
              com_diff[[1]][1:25],com_diff[[1]][26:50],com_diff[[1]][51:75],com_diff[[1]][76:100]),
        main="Symmetric Graph, Varying Noise",ylab="Mean absolute error",
        names=c(paste0("coal,",c(1/1000,1/200,1/100,1/50)),paste0("com,",c(1/1000,1/200,1/100,1/50))),las=2,ylim=c(0,3),
        cex.main=.8,cex.lab=.8,cex.axis=.8)

#clear old run values
rm(mult)

#next
sim <- 2
load("mult_noise_asym.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[2]][1:25],coal_diff[[2]][26:50],coal_diff[[2]][51:75],coal_diff[[2]][76:100],
              com_diff[[2]][1:25],com_diff[[2]][26:50],com_diff[[2]][51:75],com_diff[[2]][76:100]),
        main="Asymmetric Graph, Varying Noise",ylab="Mean absolute error",
        names=c(paste0("coal,",c(1/1000,1/200,1/100,1/50)),paste0("com,",c(1/1000,1/200,1/100,1/50))),las=2,ylim=c(0,3),
        cex.main=.8,cex.lab=.8,cex.axis=.8)
#clear old run values
rm(mult)

#next
sim <- 3
load("mult_noise_sym_vc.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[3]][1:25],coal_diff[[3]][26:50],coal_diff[[3]][51:75],coal_diff[[3]][76:100],
              com_diff[[3]][1:25],com_diff[[3]][26:50],com_diff[[3]][51:75],com_diff[[3]][76:100]),
        main="Symmetric Graph, Varying Noise, Multiple Coalescence Parameters",ylab="Mean absolute error",
        names=c(paste0("coal,",c(1/1000,1/200,1/100,1/50)),paste0("com,",c(1/1000,1/200,1/100,1/50))),las=2,ylim=c(0,3),
        cex.main=.8,cex.lab=.8,cex.axis=.8)
#clear old run values
rm(mult)

#next
sim <- 4
load("mult_noise_asym_vc.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[4]][1:25],coal_diff[[4]][26:50],coal_diff[[4]][51:75],coal_diff[[4]][76:100],
              com_diff[[4]][1:25],com_diff[[4]][26:50],com_diff[[4]][51:75],com_diff[[4]][76:100]),
        main="Asymmetric Graph, Varying Noise, Multiple Coalescence Parameters",ylab="Mean absolute error",
        names=c(paste0("coal,",c(1/1000,1/200,1/100,1/50)),paste0("com,",c(1/1000,1/200,1/100,1/50))),las=2,ylim=c(0,3),
        cex.main=.8,cex.lab=.8,cex.axis=.8)
#clear old run values
rm(mult)

#next
sim <- 5
load("mult_gam_sym.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[5]][1:25],coal_diff[[5]][26:50],coal_diff[[5]][51:75],
              com_diff[[5]][1:25],com_diff[[5]][26:50],com_diff[[5]][51:75]),
        main="Symmetric Graph, Varying Coalescence Rate",ylab="Mean absolute error",
        names=c(paste0("coal,",c(10,1,0.1)),paste0("com,",c(10,1,0.1))),las=2,ylim=c(0,3),cex.main=.8,cex.lab=.8,cex.axis=.8)

#clear old run values
rm(mult)

#next
sim <- 6
load("mult_gam_asym.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[6]][1:25],coal_diff[[6]][26:50],coal_diff[[6]][51:75],
              com_diff[[6]][1:25],com_diff[[6]][26:50],com_diff[[6]][51:75]),
        main="Asymmetric Graph, Varying Coalescence Rate",ylab="Mean absolute error",
        names=c(paste0("coal,",c(10,1,0.1)),paste0("com,",c(10,1,0.1))),las=2,ylim=c(0,3),cex.main=.8,cex.lab=.8,cex.axis=.8)

#clear old run values
rm(mult)

#next
sim <- 7
load("mult_gam_sym_vc.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[7]][1:25],coal_diff[[7]][26:50],coal_diff[[7]][51:75],
              com_diff[[7]][1:25],com_diff[[7]][26:50],com_diff[[7]][51:75]),
        main="Symmetric Graph, Varying Coalescence Rate, Multiple Coalescence Parameters",ylab="Mean absolute error",
        names=c(paste0("coal,",c(10,1,0.1)),paste0("com,",c(10,1,0.1))),las=2,ylim=c(0,3),cex.main=.8,cex.lab=.8,cex.axis=.8)

#clear old run values
rm(mult)

#next
sim <- 8
load("mult_gam_asym_vc.RData")
#initialize
coal_diff[[sim]] <- rep(0,100)
com_diff[[sim]] <- rep(0,100)
coal_iqr[[sim]] <- rep(0,100)
com_iqr[[sim]] <- rep(0,100)

for(i in 1:length(coal_diff[[sim]])){
  coal_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_coal$ans$g_med))
  com_diff[[sim]][i] <- mean(abs(mult[[i]]$g-mult[[i]]$a3x2_com$ans$g_med))
  coal_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_coal$ans$g))
  com_iqr[[sim]][i] <- mean(matrixStats::colIQRs(mult[[i]]$a3x2_com$ans$g))
}

boxplot(cbind(coal_diff[[8]][1:25],coal_diff[[8]][26:50],coal_diff[[8]][51:75],
              com_diff[[8]][1:25],com_diff[[8]][26:50],com_diff[[8]][51:75]),
        main="Asymmetric Graph, Varying Coalescence Rate, Multiple Coalescence Parameters",ylab="Mean absolute error",
        names=c(paste0("coal,",c(10,1,0.1)),paste0("com,",c(10,1,0.1))),las=2,ylim=c(0,3),cex.main=.8,cex.lab=.8,cex.axis=.8)

#clear old run values
rm(mult)

#make plots
layout(matrix(1:4,2))
par(mar=c(3.5, 2.7, 1, 0.5), mgp=c(1.5, 0.6,0))
par(cex.main=.9)

#mean absolute errors, noise plots

boxplot(cbind(coal_diff[[1]][1:25],coal_diff[[1]][26:50],coal_diff[[1]][51:75],coal_diff[[1]][76:100],
              com_diff[[1]][1:25],com_diff[[1]][26:50],com_diff[[1]][51:75],com_diff[[1]][76:100]),
        main="Symmetric Graph, Varying Noise",ylab="Mean absolute error",ylim=c(0,3),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_diff[[2]][1:25],coal_diff[[2]][26:50],coal_diff[[2]][51:75],coal_diff[[2]][76:100],
              com_diff[[2]][1:25],com_diff[[2]][26:50],com_diff[[2]][51:75],com_diff[[2]][76:100]),
        main="Asymmetric Graph, Varying Noise",ylab="Mean absolute error",las=2,ylim=c(0,3),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_diff[[3]][1:25],coal_diff[[3]][26:50],coal_diff[[3]][51:75],coal_diff[[3]][76:100],
              com_diff[[3]][1:25],com_diff[[3]][26:50],com_diff[[3]][51:75],com_diff[[3]][76:100]),
        main="Symmetric Graph, Varying Noise, Multiple Coalescence Parameters",ylab="Mean absolute error",
        ylim=c(0,3),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_diff[[4]][1:25],coal_diff[[4]][26:50],coal_diff[[4]][51:75],coal_diff[[4]][76:100],
              com_diff[[4]][1:25],com_diff[[4]][26:50],com_diff[[4]][51:75],com_diff[[4]][76:100]),
        main="Asymmetric Graph, Varying Noise, Multiple Coalescence Parameters",ylab="Mean absolute error",
        ylim=c(0,3),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

#mean absolute errors, gamma plots

boxplot(cbind(coal_diff[[5]][1:25],coal_diff[[5]][26:50],coal_diff[[5]][51:75],
              com_diff[[5]][1:25],com_diff[[5]][26:50],com_diff[[5]][51:75]),
        main="Symmetric Graph, Varying Coalescence Rate",ylab="Mean absolute error",
        ylim=c(0,3),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_diff[[6]][1:25],coal_diff[[6]][26:50],coal_diff[[6]][51:75],
              com_diff[[6]][1:25],com_diff[[6]][26:50],com_diff[[6]][51:75]),
        main="Asymmetric Graph, Varying Coalescence Rate",ylab="Mean absolute error",ylim=c(0,3),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_diff[[7]][1:25],coal_diff[[7]][26:50],coal_diff[[7]][51:75],
              com_diff[[7]][1:25],com_diff[[7]][26:50],com_diff[[7]][51:75]),
        main="Symmetric Graph, Varying Coalescence Rate, Multiple Coalescence Parameters",ylab="Mean absolute error",
        ylim=c(0,3),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_diff[[8]][1:25],coal_diff[[8]][26:50],coal_diff[[8]][51:75],
              com_diff[[8]][1:25],com_diff[[8]][26:50],com_diff[[8]][51:75]),
        main="Asymmetric Graph, Varying Coalescence Rate, Multiple Coalescence Parameters",ylab="Mean absolute error",
        ylim=c(0,3),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

#inter-quartile ranges, noise plots

boxplot(cbind(coal_iqr[[1]][1:25],coal_iqr[[1]][26:50],coal_iqr[[1]][51:75],coal_iqr[[1]][76:100],
              com_iqr[[1]][1:25],com_iqr[[1]][26:50],com_iqr[[1]][51:75],com_iqr[[1]][76:100]),
        main="Symmetric Graph, Varying Noise",ylab="Mean Interquartile Range",ylim=c(0,1.7),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_iqr[[2]][1:25],coal_iqr[[2]][26:50],coal_iqr[[2]][51:75],coal_iqr[[2]][76:100],
              com_iqr[[2]][1:25],com_iqr[[2]][26:50],com_iqr[[2]][51:75],com_iqr[[2]][76:100]),
        main="Asymmetric Graph, Varying Noise",ylab="Mean Interquartile Range",ylim=c(0,1.7),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_iqr[[3]][1:25],coal_iqr[[3]][26:50],coal_iqr[[3]][51:75],coal_iqr[[3]][76:100],
              com_iqr[[3]][1:25],com_iqr[[3]][26:50],com_iqr[[3]][51:75],com_iqr[[3]][76:100]),
        main="Symmetric Graph, Varying Noise, Multiple Coalescence Parameters",ylab="Mean Interquartile Range",
        ylim=c(0,1.7),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_iqr[[4]][1:25],coal_iqr[[4]][26:50],coal_iqr[[4]][51:75],coal_iqr[[4]][76:100],
              com_iqr[[4]][1:25],com_iqr[[4]][26:50],com_iqr[[4]][51:75],com_iqr[[4]][76:100]),
        main="Asymmetric Graph, Varying Noise, Multiple Coalescence Parameters",ylab="Mean Interquartile Range",
        ylim=c(0,1.7),xaxt='n')
abline(v=4.5, lty=3, col='red')
axis(1, at=1:8, rep(c(1/1000,1/200,1/100,1/50), 2), las=2)
axis(1, at=c(2.5, 6.5), line=1.9, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

#inter-quartile ranges, coalescence rate plots

boxplot(cbind(coal_iqr[[5]][1:25],coal_iqr[[5]][26:50],coal_iqr[[5]][51:75],
              com_iqr[[5]][1:25],com_iqr[[5]][26:50],com_iqr[[5]][51:75]),
        main="Symmetric Graph, Varying Coalescence Rate",ylab="Mean Interquartile Range",
        ylim=c(0,1.7),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_iqr[[6]][1:25],coal_iqr[[6]][26:50],coal_iqr[[6]][51:75],
              com_iqr[[6]][1:25],com_iqr[[6]][26:50],com_iqr[[6]][51:75]),
        main="Asymmetric Graph, Varying Coalescence Rate",ylab="Mean Interquartile Range",
        ylim=c(0,1.7),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_iqr[[7]][1:25],coal_iqr[[7]][26:50],coal_iqr[[7]][51:75],
              com_iqr[[7]][1:25],com_iqr[[7]][26:50],com_iqr[[7]][51:75]),
        main="Symmetric Graph, Varying Coalescence Rate, Multiple Coalescence Parameters",ylab="Mean Interquartile Range",
        ylim=c(0,1.7),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

boxplot(cbind(coal_iqr[[8]][1:25],coal_iqr[[8]][26:50],coal_iqr[[8]][51:75],
              com_iqr[[8]][1:25],com_iqr[[8]][26:50],com_iqr[[8]][51:75]),
        main="Asymmetric Graph, Varying Coalescence Rate, Multiple Coalescence Parameters",ylab="Mean Interquartile Range",
        ylim=c(0,1.7),xaxt='n')
abline(v=3.5, lty=3, col='red')
axis(1, at=1:6, rep(c(10,1,0.1), 2), las=2)
axis(1, at=c(2, 5), line=1.2, labels=c("coalescence", "commute"), tick=FALSE, lwd=0)

#means of means
for(i in 1:4){
  print(paste(i,"coal_diff",mean(coal_diff[[i]][1:25])))
  print(paste(i,"coal_iqr",mean(coal_iqr[[i]][1:25])))
}
mean(coal_diff[[1]][1:25])


load("~/Documents/data/slim/new/2d/sim_2d.txt_010279/Unif.4x4.RData")
hs_1 <- hs
load("~/Documents/data/slim/new/2d/sim_2d.txt_010279/unif_4x4_1_2.RData")
hs_2 <- hs
load("~/Documents/data/slim/new/2d/sim_2d.txt_011671/Unif.4x4.RData")
hs_3 <- hs
load("~/Documents/data/slim/new/2d/sim_2d.txt_015640/Unif.4x4.RData")
hs_4 <- hs


plot(hs_1,hs_2,pch=1,col=1)
points(hs_1,hs_3,pch=2,col=2)
points(hs_1,hs_4,pch=3,col=4)

matplot(hs_1,cbind(hs_2,hs_3,hs_4),pch=c(1,2,3),col=c(1,2,4),xlab="H Values, Sim 1, Set 1",ylab="Alternative H Values",
        main="Comparison of Different Genetic Distance Matrices")
legend(500,1300,legend=c("Sim 1, Set 2","Sim 2","Sim 3"),pch=c(1,2,3),col=c(1,2,4))
abline(0,1)

mean(abs(hs_1-hs_2)/hs_1)
mean(abs(hs_1-hs_3)/hs_1)
mean(abs(hs_1-hs_4)/hs_1)

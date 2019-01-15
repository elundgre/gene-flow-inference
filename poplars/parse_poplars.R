library(rgdal)
library(sp)
library(rgeos)
library(spdep)
library(colorspace)
library(maps)

# Download data from https://datadryad.org/resource/doi:10.5061/dryad.7s848
#   and convert linebreaks

adjacency <- matrix(0, nrow=9, ncol=9)
adjlist <- list( c(1,6), c(1,7), c(6,7), c(2,5), c(7,2), c(2,4), c(4,8), c(3,4), c(3,8), c(8,9), c(9,3), c(2,8), c(2,3), c(2,6), c(3,5), c(6,5) )
for (adj in adjlist) {
    adjacency[adj[1], adj[2]] <- adjacency[adj[2], adj[1]] <- 1
}


# Divide into groups
pdf('poplar_groups.pdf', width=12, height=12)
    pops <- read.csv("Populus_metadata.csv", header=TRUE, skip=1, stringsAsFactors=FALSE)
    pops$Latitude <- as.numeric(pops$Latitude)
    pops$Longitude <- as.numeric(pops$Longitude)
    pops$Elevation <- as.numeric(pops$Elevation)
    pops <- subset(pops, !is.na(pops$Latitude))
    pops$Species <- factor(pops$Species)

    pop_coords <- SpatialPoints(as.matrix(pops[,c("Longitude", "Latitude")]),
                                proj4string=CRS("+proj=longlat"))

    north_separator <- list(a = 147.6,
                            b = 0.7076)
    south_separator <- list(a = 143.5,
                            b = 0.7076)

    sep_fun <- function (lat, lon, ab) {
        # is above line or not?
        return(lat - ab$b * lon > ab$a)
    }

    pops$groups <- NA
    #
    this_group <- "East balsamifera"
    pops$groups[pops$Species == "Populus balsamifera"
               & pops$Longitude > (-120)] <- this_group
    #
    this_group <- "Alaska balsamifera"
    pops$groups[pops$Species == "Populus balsamifera"
               & pops$Longitude < (-141)] <- this_group
    #
    this_group <- "Northern balsamifera"
    pops$groups[pops$Species == "Populus balsamifera"
               & pops$Longitude > (-141)
               & pops$Longitude < (-120)] <- this_group
    #
    this_group <- "US trichocarpa"
    pops$groups[pops$Species == "Populus trichocarpa"
               & pops$Latitude < 47.5] <- this_group
    #
    this_group <- "South CA trichocarpa"
    pops$groups[pops$Species == "Populus trichocarpa"
               & (pops$Latitude > 47.5
                  & pops$Latitude < 52.2)
               | (pops$Longitude < -125.7
                  & !sep_fun(pops$Latitude, pops$Longitude, south_separator))] <- this_group
    #
    this_group <- "Inland trichocarpa"
    pops$groups[pops$Species == "Populus trichocarpa"
               & pops$Longitude > -125.7
               & pops$Longitude < -120
               & pops$Latitude > 52.2] <- this_group
    #
    this_group <- "Central CA trichocarpa"
    pops$groups[pops$Species == "Populus trichocarpa"
               & pops$Longitude < -125.7
               & sep_fun(pops$Latitude, pops$Longitude, south_separator)
               & ! sep_fun(pops$Latitude, pops$Longitude, north_separator)
               & pops$Latitude > 52.2] <- this_group
    #
    this_group <- "Northern CA trichocarpa"
    pops$groups[pops$Species == "Populus trichocarpa"
               & pops$Longitude < -125.7
               & sep_fun(pops$Latitude, pops$Longitude, north_separator)] <- this_group
    #
    this_group <- "Northeast balsamifera"
    pops$groups[pops$Species == "Populus balsamifera"
               & pops$Latitude > 60
               & pops$Longitude > (-120)] <- this_group
    # REMOVE the furthest-east group
    this_group <- NA
    pops$groups[pops$Longitude > -100] <- this_group

    pops$groups <- factor(pops$groups)

    centroids <- data.frame(Longitude = tapply(pops$Longitude, pops$groups, mean),
                            Latitude = tapply(pops$Latitude, pops$groups, mean),
                            group = levels(pops$groups))


    cols <- rainbow(10)
    plot(pop_coords, pch=21, cex=2, 
         col=adjustcolor(c("blue", "red")[as.numeric(pops$Species)], 0.75),
         bg=cols[as.numeric(pops$groups)],
    lwd=3)
    text(pop_coords@coords, labels=as.numeric(pops$groups))
    points(centroids[,1:2], pch=21, col=adjustcolor(cols[1:nlevels(pops$groups)], 0.25), cex=5)
    text(centroids[,1:2], labels=1:nlevels(pops$groups), col=cols[1:nlevels(pops$groups)], cex=5)
    map(add=TRUE)

    segments(x0=centroids[sapply(adjlist, '[', 1), 1],
             x1=centroids[sapply(adjlist, '[', 2), 1],
             y0=centroids[sapply(adjlist, '[', 1), 2],
             y1=centroids[sapply(adjlist, '[', 2), 2],
             lty=3, col='blue')

    abline(v = -120)
    abline(h = 60)
    abline(v = -141)
    abline(h = 47.5)
    abline(a = north_separator$a,
           b = north_separator$b)
    abline(v=-125.7)
    abline(h=52.2)
    abline(a = south_separator$a,
           b = south_separator$b)
    if (construct_data) {
        title(main="conStruct drainages")
    } else {
        title(main="full dataset")
    }
dev.off()

########
# compute pairwise distances
########

con <- pipe("tail -n +10 FOR_DRYAD_TandB_fwd_june2013_434.txt | sed -e 's/|[^\t]*//g'", open="r")
genotypes <- read.delim(con, sep="\t", stringsAsFactors=FALSE)[,-1]
close(con)

samp.indices <- match(colnames(genotypes)[-1], gsub("-", ".", pops$Accession))
stopifnot(all(!is.na(samp.indices)))
pops <- pops[samp.indices,]

# remove snps with >= 5% missing data, removing 2314 and retaining 30756
genotypes[genotypes == ""] <- NA
genotypes <- genotypes[rowMeans(is.na(genotypes)) < 0.05,]
# Look at missing data by individual
pmiss <- colMeans(is.na(genotypes)[,-1])
plot(sort(pmiss), col=as.numeric(pops$Species)[order(pmiss)])
with(pops, plot(Longitude, Latitude, cex=100*pmiss, col=Species))

# remove furthest east indiv and indiv with 15% missing
remove_these <- (pops$Longitude > -100) | (pmiss > 0.1)
pops <- pops[!remove_these,]
genotypes <- genotypes[,c(TRUE, !remove_these)]

geno1 <- do.call(cbind, lapply(genotypes[-1], substr, 1, 1))
geno2 <- do.call(cbind, lapply(genotypes[-1], substr, 2, 2))

D <- matrix(NA, nrow=ncol(genotypes)-1, ncol=ncol(genotypes)-1)
for (i in 1:ncol(D)) {
    D[i, i] <- mean(geno1[,i] != geno2[,i], na.rm=TRUE)
    for (j in (i:ncol(D))[-1]) {
        D[i, j] <- D[j, i] <- (mean(geno1[,i] != geno1[,j], na.rm=TRUE) 
                               + mean(geno1[,i] != geno2[,j], na.rm=TRUE)
                               + mean(geno2[,i] != geno1[,j], na.rm=TRUE)
                               + mean(geno2[,i] != geno2[,j], na.rm=TRUE))/4
    }
}
rownames(D) <- colnames(D) <- pops$Accession


# Look at IBD
dists <- sqrt(outer(pops$Latitude, pops$Latitude, "-")^2 
              + outer(pops$Longitude, pops$Longitude, "-")^2)
ut <- upper.tri(dists, diag=TRUE) 
cols <- D
cols[] <- NA
cols[(pops$Species[row(dists)] == "Populus trichocarpa") & (pops$Species[col(dists)] == "Populus trichocarpa")] <- adjustcolor('blue', 0.5)
cols[(pops$Species[row(dists)] == "Populus trichocarpa") & (pops$Species[col(dists)] != "Populus trichocarpa")] <- adjustcolor('red', 0.5)
cols[(pops$Species[row(dists)] != "Populus trichocarpa") & (pops$Species[col(dists)] == "Populus trichocarpa")] <- adjustcolor('red', 0.5)
cols[(pops$Species[row(dists)] != "Populus trichocarpa") & (pops$Species[col(dists)] != "Populus trichocarpa")] <- adjustcolor('green', 0.5)

png(file="populus_ibd.png", width=8*288, height=6*288, res=288)
plot(dists[ut], D[ut], col=cols[ut], pch=20, cex=0.5)
legend("topleft", pch=1, col=c("blue", "red", "green"), legend=c("tri", "tri-bal", "bal"))
dev.off()

#######
# output
#######

stopifnot(all(pops$Accession == gsub(".", "-", colnames(genotypes)[-1], fixed=TRUE)))

write.table(pops[,c("Accession", "Species", "Latitude", "Longitude", "groups")], file="populus_info.tsv", row.names=FALSE)
write.table(D, file="populus_genetic_distance.tsv")

library(rgdal)
library(sp)
library(rgeos)
library(spdep)
library(colorspace)
library(maps)


adjacency <- matrix(0, nrow=9, ncol=9)
adjlist <- list( c(1,6), c(1,7), c(6,7), c(2,5), c(7,2), c(2,4), c(4,8), c(3,4), c(3,8), c(8,9), c(9,3), c(2,8), c(2,3), c(2,6), c(3,5), c(6,5) )
for (adj in adjlist) {
    adjacency[adj[1], adj[2]] <- adjacency[adj[2], adj[1]] <- 1
}


# Divide into groups
pdf('poplar_groups.pdf', width=12, height=12)
poplist <- lapply(c('construct'=TRUE, 'all'=FALSE), function (construct_data) {
    if (construct_data) {
        load("poplar.data.Robj")
        pops <- data.frame(poplar.data$coords)
        colnames(pops) <- c("Longitude", "Latitude")
        pops$Species <- poplar.data$sp.ID
        pops$sample_size <- c(10, 40, 10, 34, 17, 50, 14, 12, 4, 32, 10, 9, 17, 29, 8, 26, 7, 26, 13, 8, 10, 23, 4, 5, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    } else {
        pops <- read.csv("Populus_metadata.csv", header=TRUE, skip=1, stringsAsFactors=FALSE)
        pops$Latitude <- as.numeric(pops$Latitude)
        pops$Longitude <- as.numeric(pops$Longitude)
        pops$Elevation <- as.numeric(pops$Elevation)
        pops <- subset(pops, !is.na(pops$Latitude))
        pops$Species <- factor(pops$Species)
    }

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
    return(list(pops=pops, centroids=centroids))
})
dev.off()


# compute pairwise distances
groups <- poplist$construct$pops$groups
load("poplar.data.Robj")
n <- tapply(poplar.data$sample.sizes, groups, sum)
P <- matrix(NA, nrow=nlevels(groups), ncol=ncol(poplar.data$freqs))
species <- ifelse(grepl("trichocarpa", levels(groups)), "trichocarpa", "balsamifera")
for (k in 1:nlevels(groups)) {
    ut <- (!is.na(groups)) & (groups == levels(groups)[k])
    P[k,] <- colSums(poplar.data$freqs[ut,,drop=FALSE] * poplar.data$sample.sizes[ut]) / n[k]
}
stopifnot(all(!is.na(P)))
D <- tcrossprod(P, 1-P) + tcrossprod(1-P, P)
diag(D) <- 2 * rowSums(P * (1-P)) * n / (n - 1)
D <- D / ncol(P)

centroids <- poplist$all$centroids
dists <- sqrt(outer(centroids[,1], centroids[,1], "-")^2 
              + outer(centroids[,1], centroids[,1], "-")^2)

# not so much IBD?
ut <- upper.tri(dists, diag=TRUE) # & (species[row(dists)] == "trichocarpa") & (species[col(dists)] == "trichocarpa")
plot(dists[ut], D[ut])

if (FALSE) {
    # check isolation by distance in the unaggregated data
    ut <- (poplar.data$sp.ID == "Populus trichocarpa")
    P <- poplar.data$freqs[ut,]
    n <- poplar.data$sample.sizes[ut]
    D <- tcrossprod(P, 1-P) + tcrossprod(1-P, P)
    diag(D) <- 2 * rowSums(P * (1-P)) * n / (n - 1)
    D <- D / ncol(P)
    centroids <- poplar.data$coords[ut,]
    dists <- sqrt(outer(centroids[,1], centroids[,1], "-")^2 
                  + outer(centroids[,1], centroids[,1], "-")^2)
    plot(dists[upper.tri(dists)],
         D[upper.tri(dists)])
}

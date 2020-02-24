# # -------------------------------------------------------------------------------
# FIGURES FOR PAPER
# EX 1 - ELECTION2005 DATASET
# EX 2 - DIAMONDS DATASET
# EX 3 - AIRQUALITY DATASET
# EX 4 - LESMIS DATASET
# EX 5 - CLASSICS FROM GUTENBERG
# -------------------------------------------------------------------------------

# CHANGE FOLDER TO DIRECTORY OF SCRIPTS AND DATA IN THE LINE BELOW
folder <- paste( getwd(), "/Supp_Mat/", sep="")
source(paste(folder, "/Functions_For_Paper.R", sep=""))

# IF DOBIN IS NOT INSTALLED, PLEASE RUN THE NEXT LINE
# devtools::install_github("sevvandi/dobin")

library("dobin")
library("OutliersO3")
library("graphics")
library("ggplot2")
library("mbgraphic")
library("gridExtra")
library("SOMbrero")
library("igraph")

# -------------------------------------------------------------------------------
# EX 1 - ELECTION2005 DATASET
# -------------------------------------------------------------------------------
data <- Election2005[, c(6, 10, 17, 28)]
names(data) <- c("Area", "Population_density", "Birthrate", "Car_ownership")

O3y <- O3prep(data, method=c("HDo", "PCS", "BAC", "adjOut", "DDC", "MCD"))
O3y1 <- O3plotM(O3y)
O3y1$gO3

# We use norm=3 because many methods in O3prep apart of 'HDo' use normalization methods similar to Median-IQR. 

out <- dobin(data, frac=0.9, norm=3) 
kk <- min(ceiling(dim(data)[1]/10),25)
knn_dist <- FNN::knn.dist(out$coords[, 1:2], k = kk)
knn_dist <- knn_dist[ ,kk]
ord <- order(knn_dist, decreasing=TRUE)
ord[1:8]
labls <- paste("X", ord[1:8], sep="")
out$vec

labs <- rep("norm", length(ord))
labs[ord[1:8]] <- "O3out"
df <- as.data.frame(out$coords[, 1:2])
colnames(df) <- c("DC1", "DC2")
df2 <- df[ord[1:8], ]
ggplot(df, aes(x=DC1,y=DC2)) + geom_point(aes(shape=labs, color=labs), size=2 ) + geom_text(data=df2, aes(DC1, DC2, label = labls), nudge_x = 0.5) + theme_bw()


# -------------------------------------------------------------------------------
# EX 2 - DIAMONDS DATASET
# -------------------------------------------------------------------------------
data(diamonds, package="ggplot2")
data <- diamonds[1:5000, c(1, 5, 6, 8:10)]
pPa <- O3prep(data, k1=5, method=c("HDo", "PCS", "adjOut"), tolHDo = 0.001, tolPCS=0.001, toladj=0.001, boxplotLimits=10)
pPx <- O3plotM(pPa)
pPx$gO3x + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))


out <- dobin(data, frac=0.9, norm=3) 
kk <- min(ceiling(dim(data)[1]/10),25)
knn_dist <- FNN::knn.dist(out$coords[, 1:3], k = kk)
knn_dist <- knn_dist[ ,kk]
ord <- order(knn_dist, decreasing=TRUE)
ord[1:4]
labls <- paste("X", ord[1:4], sep="")
out$vec

labs <- rep("norm", length(ord))
labs[ord[1:4]] <- "O3out"
df <- as.data.frame(out$coords[, 1:2])
colnames(df) <- c("DC1", "DC2")
df2 <- df[ord[1:4], ]
ggplot(df, aes(x=DC1,y=DC2)) + geom_point(aes(shape=labs, color=labs), size=2 ) + geom_text(data=df2, aes(DC1, DC2, label = labls), nudge_x = 0.7, nudge_y = -0.07) + theme_bw()



# -------------------------------------------------------------------------------
# EX 3 - AIRQUALITY DATASET
# -------------------------------------------------------------------------------
data(airquality)
colnames(airquality)
airquality2 <- airquality

# replace NA with the col mean
for(i in 1:4){
  dd <- airquality[, i]
  inds <- which(is.na(dd))
  airquality2[inds, i] <- mean(dd, na.rm=TRUE)
}

O3y <- O3prep(airquality2[, 1:4], method=c("HDo", "PCS", "BAC", "adjOut", "DDC", "MCD"), tolHDo=0.00005, tolBAC=10^(-20))
O3y1 <- O3plotM(O3y)
O3y1$gO3x  + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
par(mfrow=c(1,1))

out <- dobin(airquality2[, 1:4], frac=0.9, norm=3) 
kk <- min(ceiling(dim(data)[1]/10),25)
knn_dist <- FNN::knn.dist(out$coords[, 1:2], k = kk)
knn_dist <- knn_dist[ ,kk]
ord <- order(knn_dist, decreasing=TRUE)
ord[1]
labls <- paste("X", ord[1], sep="" )
out$vec
colnames(airquality2[, 1:4])

tt <- 1
labs <- rep("norm", length(ord))
labs[ord[1:tt]] <- "O3out"
df <- as.data.frame(out$coords[, 1:2])
colnames(df) <- c("DC1", "DC2")
df2 <- df[ord[1:tt], ]
ggplot(df, aes(x=DC1,y=DC2)) + geom_point(aes(shape=labs, color=labs), size=2 ) + geom_text(data=df2, aes(DC1, DC2, label = labls), nudge_x = 0.2) + theme_bw()


# -------------------------------------------------------------------------------
# EX 4 - LESMIS DATASET
# -------------------------------------------------------------------------------
data(lesmis)
igraph.options(plot.layout=layout.circle, vertex.size=10)
plot(lesmis,vertex.label.color="black", vertex.label.dist=1, vertex.edge.dist=20, vertex.size=4, vertex.label.cex=0.7)

# Data preparation
centrality <- alpha_centrality(lesmis, alpha=0.9)
# https://igraph.org/r/doc/alpha_centrality.html

trans <- transitivity(lesmis, type="local", vids=V(lesmis))
trans[which(is.nan(trans))] <- 0
# https://igraph.org/r/doc/transitivity.html

close <- closeness(lesmis)
between <- betweenness(lesmis)
deg <- degree(lesmis)
neigh <- knn(lesmis)
neigh <- neigh$knn
pgrk <- page_rank(lesmis)
pgrk <- pgrk$vector

dat <- cbind.data.frame(centrality,trans, close, between, deg, neigh, pgrk)
dat <- as.matrix(dat)
O3s <- O3prep(dat, method=c("HDo",  "BAC", "MCD"))
O3s1 <- O3plotM(O3s, caseNames=V(lesmis)$label )
O3s1$gO3x + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))

out <- dobin(dat, frac=0.9, norm=3) 
kk <- min(ceiling(dim(dat)[1]/10),25)
kdist <- FNN::knn.dist(out$coords[ ,1:3], k=kk)
kdist25 <- kdist[ ,kk]
ord <- order(kdist25, decreasing=TRUE)
pts <- ord[1:5]
V(lesmis)$label[pts]
out$vec[,1]
colnames(dat)

tt <- 5
labs <- rep("norm", length(ord))
labs[ord[1:tt]] <- "O3out"
df <- as.data.frame(out$coords[, 1:2])
colnames(df) <- c("DC1", "DC2")
df2 <- df[ord[1:tt], ]
ggplot(df, aes(x=DC1,y=DC2)) + geom_point(aes(shape=labs, color=labs), size=2 ) + geom_text(data=df2, aes(DC1, DC2, label = V(lesmis)$label[ord[1:tt]]), nudge_x = -2, nudge_y=0.3) + theme_bw()


# -------------------------------------------------------------------------------
# EX 5 - CLASSICS FROM GUTENBERG
# -------------------------------------------------------------------------------
# folder <- paste( getwd(), "/data/", sep="")
load(file= paste(folder, "book_vectors.rda", sep=""))
book_titles <- books[ ,1]
books <- books[ ,-1] # removing title

# Use PCA as a change of coordinates
pca <- prcomp(books, scale=TRUE, center=TRUE) 
# Use all PCA components
X <- pca$x
nbooks <- dim(pca$x)[2]

colnames(X) <- c(paste("V", 1:nbooks, sep=""))
O3y <- O3prep(X, k1=20, method=c("HDo"), tolHDo = 0.01)
O3y1 <- O3plotT(O3y, caseNames = strtrim(book_titles,20))
O3y1$gO3 + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))


kk <- min(ceiling(dim(books)[1]/10),25)
out <- dobin(books, frac=0.95, norm=1)  
knndist <- FNN::knn.dist(out$coords[,1:3], k=kk)
knndist <- knndist[,kk]
ord <- order(knndist, decreasing=TRUE)
book_titles[ord[1]]


tt <- 1
labs <- rep("norm", length(ord))
labs[ord[1:tt]] <- "O3out"
df <- as.data.frame(out$coords[, 1:2])
colnames(df) <- c("DC1", "DC2")
df2 <- df[ord[1:tt], ]
ggplot(df, aes(x=DC1,y=DC2)) + geom_point(aes(shape=labs, color=labs), size=2 ) + geom_label(data=df2, aes(DC1, DC2, label = book_titles[ord[1:tt]]), nudge_y = -3, nudge_x=-3) +coord_equal() + theme_bw()



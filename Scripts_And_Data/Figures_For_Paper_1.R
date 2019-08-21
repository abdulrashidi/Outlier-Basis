# -------------------------------------------------------------------------------
# FIGURES FOR PAPER
# FIG 1 - X AND Y SPACES
# EXP 1 - 2 NORMAL DISTRIBUTIONS WITH ONE MOVING OUT
# EXP 2 - 3 NORMAL DISTRIBUTIONS - BIMODAL -  WITH ONE MOVING INTO THE TROUGH
# EXP 3 - UNIFORM DISTRIBUTION IN 20 DIMENSIONS - OUTLIER AT 0.9 IN i DIMENSIONS
# -------------------------------------------------------------------------------

# CHANGE FOLDER TO DIRECTORY OF SCRIPTS AND DATA IN THE LINE BELOW
folder <- paste( getwd(), "/paper/", sep="")
source(paste(folder, "Functions_For_Paper.R", sep=""))

# To install dobin
devtools::install_github("sevvandi/dobin")
library("dobin")
library("ggplot2")
library("IsolationForest")
library("DMwR")
library("ggplot2")
library("pROC")
library("tidyr")

# -------------------------------------------------------------------------------
# FIG 1 - X AND Y SPACES
# -------------------------------------------------------------------------------
set.seed(1)
x1 <- rnorm(100)
x2 <- rnorm(100)
X <- cbind(x1, x2)
X <- rbind(X, c(2,5))
X <- as.data.frame(X)
labs <- c(rep("black",100), "red")
ggplot(X, aes(x1,x2)) + geom_point(color =labs) + theme_bw()
out <- dobin(X)
Y <- out$Y
Ypairs <- out$Ypairs
labs <- rep("black", dim(Ypairs)[1])
labs[which(Ypairs %in% 101)] <- "red"
ggplot(as.data.frame(Y), aes(x1,x2)) + geom_point(color=labs) + theme_bw() + xlab("y1") + ylab("y2")

# -------------------------------------------------------------------------------
# EXP 1 - 2 NORMAL DISTRIBUTIONS WITH ONE MOVING OUT
# -------------------------------------------------------------------------------

values <- rep(0, 10)
pp <- 10

roc_u1 <- roc_sc <- roc_x1 <- matrix(0, nrow=pp, ncol=10)
knn_x1 <- knn_u1 <- knn_sc <- matrix(0, nrow=pp, ncol=10)
pca_knn <-  pca_lof <- matrix(0, nrow=pp, ncol=10)
ifor_u1 <- ifor_sc <- ifor_pca <-  matrix(0, nrow=pp, ncol=10)


kkk <- min(max( ceiling(405/50),2), 25)
set.seed(123)
for(kk in 1:pp){
  x2 <- rnorm(405)
  x3 <- rnorm(405)
  x4 <- rnorm(405)
  x5 <- rnorm(405)
  x6 <- rnorm(405)
  x1_1 <- rnorm(400)
  
  for(i in 1:10){
    mu2 <- (i-1)*0.5
    x1_2 <- rnorm(5, mean=mu2, sd=0.2)
    x1 <- c(x1_1, x1_2)
    X <- cbind(x1,x2,x3,x4,x5,x6)
    labs <- c(rep(0,400), rep(1,5))
    
    # NEW SCALING - START
    # TAKE 5 - knncomp1
    out <- dobin(X, frac=0.95, norm=1, vis=FALSE)
    XSc <- out$coords[ ,1:3]
    
    # PCA
    pca <- prcomp(X, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:3]
    

    # Test Min-Max Scaling - LOF
    Xu <- apply(X, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    roc_u1[kk,i] <- roc_obj1$auc
    
    # Test New Scaling - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    roc_sc[kk,i] <- roc_obj2$auc
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    

    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    knn_u1[kk,i] <- roc_obj3$auc
    
    # Test New Scaling - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    knn_sc[kk,i] <- roc_obj4$auc
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu)
    #evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_u1[kk, i] <- roc_obj$auc
    
    
    # Test New Scaling - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc)
    #evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_sc[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_pca[kk, i] <- roc_obj$auc
    
  }
}

lab <- c("mu", "AUC")
gg_title <- "Area under ROC - LOF "
mu <- seq(0, 4.5, by =0.5)
df_lof <- draw_three_curves(roc_u1, roc_sc, pca_lof,  mu, lab, gg_title, 3)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN"
df_knn <- draw_three_curves(knn_u1, knn_sc, pca_knn,  mu, lab, gg_title, 3)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest"
df_ifor <- draw_three_curves(ifor_u1, ifor_sc, ifor_pca,  mu, lab, gg_title, 3)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"


gg_title <- "Area under ROC for LOF, KNN and iForest"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab(expression(mu)) + ylab(lab[2]) + facet_grid(.~Method)   + theme_bw()



# -------------------------------------------------------------------------------
# EXP 2 - 3 NORMAL DISTRIBUTIONS - BIMODAL -  WITH ONE MOVING INTO THE TROUGH
# -------------------------------------------------------------------------------
values <- rep(0, 10)
pp <- 10

roc_u1 <- roc_sc <- roc_x1 <- matrix(0, nrow=pp, ncol=10)
knn_x1 <- knn_u1 <- knn_sc <- matrix(0, nrow=pp, ncol=10)
pca_knn <-  pca_lof <- matrix(0, nrow=pp, ncol=10)
ifor_u1 <- ifor_sc <- ifor_pca <-  matrix(0, nrow=pp, ncol=10)

kkk <- min(max( ceiling(805/50),2), 25)
set.seed(123)
for(kk in 1:pp){
  x2 <- rnorm(405)
  x3 <- rnorm(405)
  x4 <- rnorm(405)
  x5 <- rnorm(405)
  x6 <- rnorm(405)
  x1_1 <- rnorm(mean = 5, 400)
  
  for(i in 1:10){
    mu2 <-  5 - (i-1)*0.5
    x1_2 <- rnorm(5, mean=mu2, sd=0.2)
    x1 <- c(x1_1, x1_2)
    X1 <- cbind(x1,x2,x3,x4,x5,x6)
    X2 <- cbind(-1*x1_1,x2[1:400],x3[1:400],x4[1:400],x5[1:400],x6[1:400])
    X <- rbind(X1, X2)
    labs <- c(rep(0,400), rep(1,5), rep(0,400))
    

    # TAKE 4 - knncomp
    out <- dobin(X, frac=0.95, norm=1, vis=FALSE)
    XSc <- out$coords[ ,1:3]
    
    # PCA
    pca <- prcomp(X, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:3]
    

    # Test Min-Max Scaling - LOF
    Xu <- apply(X, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    roc_u1[kk,i] <- roc_obj1$auc
    
    # Test New Scaling - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    roc_sc[kk,i] <- roc_obj2$auc
    
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    knn_u1[kk,i] <- roc_obj3$auc
    
    # Test New Scaling - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    knn_sc[kk,i] <- roc_obj4$auc
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu)
    #evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_u1[kk, i] <- roc_obj$auc
    
    
    # Test New Scaling - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc)
    #evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_sc[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_pca[kk, i] <- roc_obj$auc
    
  }
}

lab <- c("mu", "AUC")
gg_title <- "Area under ROC - LOF "
mu <- seq(0, 4.5, by =0.5)
df_lof <- draw_three_curves(roc_u1, roc_sc, pca_lof,  mu, lab, gg_title, 3)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN"
df_knn <- draw_three_curves(knn_u1, knn_sc, pca_knn,  mu, lab, gg_title, 3)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest"
df_ifor <- draw_three_curves(ifor_u1, ifor_sc, ifor_pca,  mu, lab, gg_title, 3)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"


gg_title <- "Area under ROC for LOF, KNN and iForest"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab(expression(mu)) + ylab(lab[2]) + facet_grid(.~Method)   + theme_bw()


# -------------------------------------------------------------------------------
# EXP 3 - UNIFORM DISTRIBUTION IN 20 DIMENSIONS - OUTLIER AT 0.9 IN i DIMENSIONS
# -------------------------------------------------------------------------------

dd <- 20
pp <- 10
roc_u1 <- roc_sc <-  matrix(0, nrow=pp, ncol=dd)
knn_u1 <- knn_sc <- matrix(0, nrow=pp, ncol=dd)
pca_knn <-  pca_lof <- matrix(0, nrow=pp, ncol=dd)
ifor_u1 <- ifor_sc <- ifor_pca <-  matrix(0, nrow=pp, ncol=dd)

kkk <- min(max( ceiling(500/50),2), 25)
set.seed(1)

for(kk in 1:pp){
  x <- runif(500)
  labs <- c(rep(0, 499), 1)
  
  for(i in 1:dd){
    x <- cbind(x, runif(500))
  }
  
  for(i in 1:dd){
    x[500, 1:i] <- rep(0.9, i)
    Xu <- apply(x, 2, unitize_1)


    # TAKE 4 - knncomp
    out <- dobin(x, frac=0.95, norm=1, vis=FALSE)
    XSc <- out$coords[ ,1:10]
    
    # PCA
    pca <- prcomp(x, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:10]
    
    
    # Test Min-Max Scaling - LOF
    Xu <- apply(x, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    roc_u1[kk,i] <- roc_obj1$auc
    
    # Test New Scaling - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    roc_sc[kk,i] <- roc_obj2$auc
    
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    knn_u1[kk,i] <- roc_obj3$auc
    
    # Test New Scaling - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    knn_sc[kk,i] <- roc_obj4$auc
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu, ntree=20)
    #evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_u1[kk, i] <- roc_obj$auc
    
    
    # Test New Scaling - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc, ntree=20, colSamp = TRUE, nColSamp = 5, colWeight = 10:1  )
    #evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_sc[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x, ntree=20, colSamp = TRUE, nColSamp = 5, colWeight = 10:1  )
    #evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    ifor_pca[kk, i] <- roc_obj$auc
    
  }
}

lab <- c("Outlier Dimension", "AUC")
gg_title <- "Area under ROC - LOF "
mu <- seq(1, 20, by =1)
df_lof <- draw_three_curves(roc_u1, roc_sc, pca_lof,  mu, lab, gg_title, 10)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN"
df_knn <- draw_three_curves(knn_u1, knn_sc, pca_knn,  mu, lab, gg_title, 10)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest"
df_ifor <- draw_three_curves(ifor_u1, ifor_sc, ifor_pca,  mu, lab, gg_title, 10)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"


gg_title <- "Area under ROC for LOF, KNN and iForest"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab("Outlier Dimension") + ylab(lab[2]) + facet_grid(.~Method)   + theme_bw()
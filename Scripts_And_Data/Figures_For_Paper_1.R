# -------------------------------------------------------------------------------
# FIGURES FOR PAPER
# FIG 1 - X AND Y SPACES
# EXP 1 - 2 NORMAL DISTRIBUTIONS WITH ONE MOVING OUT
# EXP 2 - 3 NORMAL DISTRIBUTIONS - BIMODAL -  WITH ONE MOVING INTO THE TROUGH
# EXP 3 - UNIFORM DISTRIBUTION IN 20 DIMENSIONS - OUTLIER AT 0.9 IN i DIMENSIONS
# -------------------------------------------------------------------------------
# CHANGE FOLDER TO DIRECTORY OF SCRIPTS AND DATA IN THE LINE BELOW
folder <- paste(getwd(), "/Supp_Mat/", sep="")
source(paste(folder, "Functions_For_Paper.R", sep=""))

# IF DOBIN IS NOT INSTALLED, PLEASE RUN THE NEXT LINE
# devtools::install_github("sevvandi/dobin")

library("dobin")
library("ggplot2")
library("IsolationForest")
library("DMwR")
library("pROC")
library("tidyr")
library("ICS")

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
out <- dobin(X, norm=3)
Y <- out$Y
Ypairs <- out$Ypairs
labs <- rep("black", dim(Ypairs)[1])
labs[which(Ypairs %in% 101)] <- "red"
ggplot(as.data.frame(Y), aes(x1,x2)) + geom_point(color=labs)  + theme_bw() + xlab("y1") + ylab("y2")



# -------------------------------------------------------------------------------
# EXP 1 - 2 NORMAL DISTRIBUTIONS WITH ONE MOVING OUT
# -------------------------------------------------------------------------------

values <- rep(0, 10)
pp <- 10

u1_lof <- u1_knn <- u1_ifor <- matrix(0, nrow=pp, ncol=10)
pca_lof <- pca_knn <- pca_ifor <- matrix(0, nrow=pp, ncol=10)
pca2_lof <- pca2_knn <- pca2_ifor <- matrix(0, nrow=pp, ncol=10)
cov4_lof <- cov4_knn <- cov4_ifor <- matrix(0, nrow=pp, ncol=10)
dob_lof <- dob_knn <- dob_ifor <- matrix(0, nrow=pp, ncol=10)

  
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
    # ---------------------------------------------
    # COORDINATES 
    # DOBIN 
    out <- dobin(X, frac=0.95, norm=1)
    XSc <- out$coords[ ,1:3]
    
    # PCA 
    pca <- prcomp(X, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:3]
    pca_x2 <- pca$x[ ,4:6]

    # ICS package - cov4
    icscov4 <- cov4(X)
    cov4_x <- solve(icscov4)%*%t(X)
    cov4_x <- t(cov4_x[1:3,])
    
    # ---------------------------------------------
    
    # ---------------------------------------------
    # LOF
    
    # Test Min-Max Scaling - LOF
    Xu <- apply(X, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    u1_lof[kk,i] <- roc_obj1$auc
    
    
    # Test DOBIN - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    dob_lof[kk,i] <- roc_obj2$auc
    

    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    
    # Test PCA 2 - LOF
    lof3 <- lofactor(pca_x2, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca2_lof[kk,i] <- roc_obj$auc
    

    # Test COV4 - LOF
    lof4 <- lofactor(cov4_x, k=kkk)
    roc_obj <- roc(labs, lof4, direction="<")
    cov4_lof[kk,i] <- roc_obj$auc
    
    # ---------------------------------------------

    # ---------------------------------------------
    # KNN
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    u1_knn[kk,i] <- roc_obj3$auc
    
    
    # Test DOBIN - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    dob_knn[kk,i] <- roc_obj4$auc


    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    
    # Test PCA 2 - KNN
    knndist3 <- FNN::knn.dist(pca_x2, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca2_knn[kk,i] <- roc_obj$auc
    

    # Test COV4 - KNN
    knndist4 <- FNN::knn.dist(cov4_x, k=kkk)
    knndist4 <- knndist4[ ,kkk]
    roc_obj <- roc(labs, knndist4, direction="<")
    cov4_knn[kk,i] <- roc_obj$auc
    
    # ---------------------------------------------
    
    # ---------------------------------------------
    # isolation forest
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu)
    #evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    u1_ifor[kk, i] <- roc_obj$auc
    
    
    # Test DOBIN - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc)
    #evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    dob_ifor[kk, i] <- roc_obj$auc
    
 
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    pca_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca 2 - isolationForest
    pca_x2 <- as.data.frame(pca_x2)
    tr <-IsolationTrees(pca_x2)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x2,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    pca2_ifor[kk, i] <- roc_obj$auc
    

    # Test COV4 - isolationForest
    cov4_x <- as.data.frame(cov4_x)
    tr <-IsolationTrees(cov4_x)
    #evaluate anomaly score
    as <-AnomalyScore(cov4_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    cov4_ifor[kk, i] <- roc_obj$auc

    # ---------------------------------------------
  }
}

lab <- c("mu", "AUC")
gg_title <- "Area under ROC - LOF "
mu <- seq(1, 10, by =1)
df_lof <- draw_five_curves(u1_lof, dob_lof, pca_lof,  pca2_lof, cov4_lof,  mu, lab, gg_title, 3)

df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN"
df_knn <- draw_five_curves(u1_knn, dob_knn,  pca_knn,  pca2_knn, cov4_knn,  mu, lab, gg_title, 3)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest"
df_ifor <- draw_five_curves(u1_ifor,  dob_ifor, pca_ifor,  pca2_ifor, cov4_ifor, mu, lab, gg_title, 3)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"

gg_title <- "Experiment 1"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5)  + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab("Iteration") + ylab(lab[2]) + facet_grid(.~Method)  + scale_x_discrete(limits=1:10) +  ylim(0.25, 1)  + theme_bw()



# -------------------------------------------------------------------------------
# EXP 2 - 3 NORMAL DISTRIBUTIONS - BIMODAL -  WITH ONE MOVING INTO THE TROUGH
# -------------------------------------------------------------------------------
values <- rep(0, 10)
pp <- 10


u1_lof <- u1_knn <- u1_ifor <- matrix(0, nrow=pp, ncol=10)
pca_lof <- pca_knn <- pca_ifor <- matrix(0, nrow=pp, ncol=10)
pca2_lof <- pca2_knn <- pca2_ifor <- matrix(0, nrow=pp, ncol=10)
cov4_lof <- cov4_knn <- cov4_ifor <- matrix(0, nrow=pp, ncol=10)
dob_lof <- dob_knn <- dob_ifor <- matrix(0, nrow=pp, ncol=10)

kkk <- min(max( ceiling(805/50),2), 25)
set.seed(123)
for(kk in 1:pp){

  X <- data.frame(
    x1 = c(rnorm(405,mean=5), rnorm(400, mean=-5)),
    x2 = rnorm(805),
    x3 = rnorm(805),
    x4 = rnorm(805),
    x5 = rnorm(805),
    x6 = rnorm(805)
  )
  labs <- c(rep(0,400), rep(1,5), rep(0,400))
  
  for(i in 1:10){
    mu2 <-  5 - (i-1)*0.5
    x1_2 <- rnorm(5, mean=mu2, sd=0.2)
    X[401:405, 1] <- x1_2

    #----------------------------------------------------
    # COORDINATES 
    
    # DOBIN
    out <- dobin(X, frac=0.95, norm=1)
    XSc <- out$coords[ ,1:3]
    
    # PCA
    pca <- prcomp(X, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:3]
    pca_x2 <- pca$x[ ,4:6]
    

    # ICS package - cov4
    icscov4 <- cov4(X)
    cov4_x <- solve(icscov4)%*%t(X)
    cov4_x <- t(cov4_x[1:3,])
    #----------------------------------------------------
    
    #----------------------------------------------------
    # LOF
    
    # Test Min-Max Scaling - LOF
    Xu <- apply(X, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    u1_lof[kk,i] <- roc_obj1$auc
    
    # Test DOBIN - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    dob_lof[kk,i] <- roc_obj2$auc
    
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    
    # Test PCA2 - LOF
    lof3 <- lofactor(pca_x2, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca2_lof[kk,i] <- roc_obj$auc

        
    # Test COV4 - LOF
    lof4 <- lofactor(cov4_x, k=kkk)
    roc_obj <- roc(labs, lof4, direction="<")
    cov4_lof[kk,i] <- roc_obj$auc
    #----------------------------------------------------
    
    
    #----------------------------------------------------
    # KNN
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    u1_knn[kk,i] <- roc_obj3$auc
    
    
    # Test DOBIN - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    dob_knn[kk,i] <- roc_obj4$auc
    
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc

    
    # Test PCA2 - KNN
    knndist3 <- FNN::knn.dist(pca_x2, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca2_knn[kk,i] <- roc_obj$auc
    
        
    # Test COV4 - KNN
    knndist4 <- FNN::knn.dist(cov4_x, k=kkk)
    knndist4 <- knndist4[ ,kkk]
    roc_obj <- roc(labs, knndist4, direction="<")
    cov4_knn[kk,i] <- roc_obj$auc
    #----------------------------------------------------
    
    #----------------------------------------------------
    # iForest
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu)
    # evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    u1_ifor[kk, i] <- roc_obj$auc
    
    
    # Test DOBIN - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc)
    #evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    dob_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    pca_ifor[kk, i] <- roc_obj$auc

    
    # Test pca2 - isolationForest
    pca_x2 <- as.data.frame(pca_x2)
    tr <-IsolationTrees(pca_x2)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x2,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    pca2_ifor[kk, i] <- roc_obj$auc
        

    # Test COV4 - isolationForest
    cov4_x <- as.data.frame(cov4_x)
    tr <-IsolationTrees(cov4_x)
    #evaluate anomaly score
    as <-AnomalyScore(cov4_x,tr)
    # show anomaly score
    roc_obj <- roc(labs, as$outF, direction="<")
    cov4_ifor[kk, i] <- roc_obj$auc
    
  }
}

lab <- c("mu", "AUC")
gg_title <- "Area under ROC - LOF - EXP2 "
mu <- seq(1, 10, by =1)
df_lof <- draw_five_curves(u1_lof, dob_lof, pca_lof,  pca2_lof, cov4_lof,  mu, lab, gg_title, 3)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN - EXP2"
df_knn <- draw_five_curves(u1_knn, dob_knn,  pca_knn,  pca2_knn, cov4_knn, mu, lab, gg_title, 3)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest - EXP2"
df_ifor <- draw_five_curves(u1_ifor,  dob_ifor, pca_ifor,  pca2_ifor, cov4_ifor, mu, lab, gg_title, 3)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"

gg_title <- "Experiment 2"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab("Iteration") + ylab(lab[2]) + facet_grid(.~Method) + scale_x_discrete(limits=1:10) + ylim(0.25, 1) +   theme_bw()


# -------------------------------------------------------------------------------
# EXP 3 - UNIFORM DISTRIBUTION IN 20 DIMENSIONS - OUTLIER AT 0.9 IN i DIMENSIONS
# -------------------------------------------------------------------------------

dd <- 20
pp <- 10

u1_lof <- u1_knn <- u1_ifor <- matrix(0, nrow=pp, ncol=20)
pca_lof <- pca_knn <- pca_ifor <- matrix(0, nrow=pp, ncol=20)
pca2_lof <- pca2_knn <- pca2_ifor <- matrix(0, nrow=pp, ncol=20)
cov4_lof <- cov4_knn <- cov4_ifor <- matrix(0, nrow=pp, ncol=20)
dob_lof <- dob_knn <- dob_ifor <- matrix(0, nrow=pp, ncol=20)


kkk <- min(max( ceiling(500/50),2), 25)
set.seed(1)
sdids <- sample( 1:9999999, pp)

for(kk in 1:pp){
  set.seed(sdids[kk])
  x <- runif(500)
  for(j in 1:dd){
    x <- cbind(x, runif(500))
  }
  labs <- c(rep(0, 499), 1)
  
  for(i in 1:dd){   
    x[500, 1:i] <- rep(0.9, i)
    Xu <- apply(x, 2, unitize_1)
    
    # -------------------------------------------------------
    # COORDINATES
    
    # DOBIN
    out <- dobin(x, frac=0.95, norm=1)
    XSc <- out$coords[ ,1:10]
    
    # PCA
    pca <- prcomp(x, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:10]
    pca_x2 <- pca$x[ ,11:20]
    

    # ICS package - cov4
    icscov4 <- cov4(x)
    cov4_x <- solve(icscov4)%*%t(x)
    cov4_x <- t(cov4_x[1:10,])
    # -------------------------------------------------------

    # -------------------------------------------------------
    # LOF
    
    # Test Min-Max Scaling - LOF
    Xu <- apply(x, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    u1_lof[kk,i] <- roc_obj1$auc
    
    # Test DOBIN - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    dob_lof[kk,i] <- roc_obj2$auc
    
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    
    # Test PCA2 - LOF
    lof3 <- lofactor(pca_x2, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca2_lof[kk,i] <- roc_obj$auc
    

    # Test COV4 - LOF
    lof4 <- lofactor(cov4_x, k=kkk)
    roc_obj <- roc(labs, lof4, direction="<")
    cov4_lof[kk,i] <- roc_obj$auc
    # -------------------------------------------------------
    
    # -------------------------------------------------------
    # KNN
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    u1_knn[kk,i] <- roc_obj3$auc
    
    # Test DOBIN - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    dob_knn[kk,i] <- roc_obj4$auc
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    # Test PCA2 - KNN
    knndist3 <- FNN::knn.dist(pca_x2, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca2_knn[kk,i] <- roc_obj$auc
    

    # Test COV4 - KNN
    knndist4 <- FNN::knn.dist(cov4_x, k=kkk)
    knndist4 <- knndist4[ ,kkk]
    roc_obj <- roc(labs, knndist4, direction="<")
    cov4_knn[kk,i] <- roc_obj$auc

    # -------------------------------------------------------
    
    # -------------------------------------------------------
    # iForest
  
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu, ntree=20) 
    #evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    u1_ifor[kk, i] <- roc_obj$auc
    
    
    # Test DOBIN - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc, ntree=20, colSamp = TRUE, nColSamp = 5, colWeight = 10:1 )
    #evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    dob_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x, ntree=20, colSamp = TRUE, nColSamp = 5, colWeight = 10:1 )
    # evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    pca_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca2 - isolationForest
    pca_x2 <- as.data.frame(pca_x2)
    tr <-IsolationTrees(pca_x2, ntree=20, colSamp = TRUE, nColSamp = 5, colWeight = 1:10 )
    #evaluate anomaly score
    as <-AnomalyScore(pca_x2,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    pca2_ifor[kk, i] <- roc_obj$auc
    
    
    # Test COV4 - isolationForest
    cov4_x <- as.data.frame(cov4_x)
    tr <-IsolationTrees(cov4_x, ntree=20, colSamp = TRUE, nColSamp = 5, colWeight = 10:1 ) 
    #evaluate anomaly score
    as <-AnomalyScore(cov4_x,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    cov4_ifor[kk, i] <- roc_obj$auc
    
    # -------------------------------------------------------
    
    
  }
}

lab <- c("Outlier Dimension", "AUC")
gg_title <- "Area under ROC - LOF - EXP3 "
mu <- seq(1, 20, by =1)


df_lof <- draw_five_curves(u1_lof, dob_lof, pca_lof, pca2_lof, cov4_lof,  mu, lab, gg_title, 10)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN - EXP3"
df_knn <- draw_five_curves(u1_knn, dob_knn,  pca_knn, pca2_knn, cov4_knn, mu, lab, gg_title, 10)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest - EXP3"
df_ifor <- draw_five_curves(u1_ifor,  dob_ifor, pca_ifor, pca2_ifor, cov4_ifor, mu, lab, gg_title, 10)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"


gg_title <- "Experiment 3"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab("Iteration") + ylab(lab[2]) + facet_grid(.~Method) + coord_fixed(ratio=62) + theme_bw()




# -------------------------------------------------------------------------------
# EXP 4 - ANNULUS EXAMPLE -  WITH ONE MOVING INTO THE CENTRE
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# PLOT THIS DATA
# -------------------------------------------------------------------------------
r1 <-runif(805)
r2 <-rnorm(805, mean=5)
theta = 2*pi*r1;
R1 <- 2
R2 <- 2
dist = r2+R2;
x =  dist * cos(theta) 
y =  dist * sin(theta) 

X <- data.frame(
  x1 = x,
  x2 = y,
  x3 = rnorm(805),
  x4 = rnorm(805),
  x5 = rnorm(805),
  x6 = rnorm(805),
  x7 = rnorm(805),
  x8 = rnorm(805),
  x9 = rnorm(805),
  x10 = rnorm(805),
  x11 = rnorm(805),
  x12 = rnorm(805),
  x13 = rnorm(805),
  x14 = rnorm(805)
)
labs <- as.factor(c(rep(0,800), rep(1,5)))
i <- 3
mu2 <-  0.9 - (i-1)*0.1
z <- cbind(rnorm(5,mu2, sd=0.1), rnorm(5,0, sd=0.1))

X[801:805, 1:2] <- z
ggplot(X[ ,1:2], aes(x1,x2)) + geom_point(aes(color=labs)) + theme_bw()


# -------------------------------------------------------------------------------
# EXPERIMENT
# -------------------------------------------------------------------------------
values <- rep(0, 10)
pp <- 10


u1_lof <- u1_knn <- u1_ifor <- matrix(0, nrow=pp, ncol=10)
pca_lof <- pca_knn <- pca_ifor <- matrix(0, nrow=pp, ncol=10)
pca2_lof <- pca2_knn <- pca2_ifor <- matrix(0, nrow=pp, ncol=10)
cov4_lof <- cov4_knn <- cov4_ifor <- matrix(0, nrow=pp, ncol=10)
dob_lof <- dob_knn <- dob_ifor <- matrix(0, nrow=pp, ncol=10)

hh <- 7
kkk <- min(max( ceiling(805/50),2), 25)
set.seed(123)
for(kk in 1:pp){
  
  r1 <-runif(805)
  r2 <-rnorm(805, mean=5)
  theta = 2*pi*r1;
  R1 <- 2
  R2 <- 2
  dist = r2+R2;
  x =  dist * cos(theta) 
  y =  dist * sin(theta) 
  
  X <- data.frame(
    x1 = x,
    x2 = y,
    x3 = rnorm(805),
    x4 = rnorm(805),
    x5 = rnorm(805),
    x6 = rnorm(805),
    x7 = rnorm(805),
    x8 = rnorm(805),
    x9 = rnorm(805),
    x10 = rnorm(805),
    x11 = rnorm(805),
    x12 = rnorm(805),
    x13 = rnorm(805),
    x14 = rnorm(805)
  )
  labs <- c(rep(0,800), rep(1,5))
  
  
  
  for(i in 1:10){
    # mu2 <-  0.9 - (i-1)*0.1
    mu2 <-  5 - (i-1)*0.5
    z <- cbind(rnorm(5,mu2, sd=0.1), rnorm(5,0, sd=0.1))
    
    X[801:805, 1:2] <- z
    #----------------------------------------------------
    # COORDINATES 
    
    # DOBIN
    out <- dobin(X, frac=0.95, norm=1)
    XSc <- out$coords[ ,1:hh]
    
    # PCA
    pca <- prcomp(X, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:hh]
    pca_x2 <- pca$x[ ,(hh+1):(2*hh)]
    
    # ICS package - cov4
    icscov4 <- cov4(X)
    cov4_x <- solve(icscov4)%*%t(X)
    cov4_x <- t(cov4_x[1:hh,])
    #----------------------------------------------------
    
    #----------------------------------------------------
    # LOF
    
    # Test Min-Max Scaling - LOF
    Xu <- apply(X, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    u1_lof[kk,i] <- roc_obj1$auc
    
    # Test DOBIN - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    dob_lof[kk,i] <- roc_obj2$auc
    
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    
    # Test PCA2 - LOF
    lof3 <- lofactor(pca_x2, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca2_lof[kk,i] <- roc_obj$auc
    
    
    # Test COV4 - LOF
    lof4 <- lofactor(cov4_x, k=kkk)
    roc_obj <- roc(labs, lof4, direction="<")
    cov4_lof[kk,i] <- roc_obj$auc
    #----------------------------------------------------
    
    
    #----------------------------------------------------
    # KNN
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    u1_knn[kk,i] <- roc_obj3$auc
    
    
    # Test DOBIN - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    dob_knn[kk,i] <- roc_obj4$auc
    
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    
    # Test PCA2 - KNN
    knndist3 <- FNN::knn.dist(pca_x2, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca2_knn[kk,i] <- roc_obj$auc
    
    
    # Test COV4 - KNN
    knndist4 <- FNN::knn.dist(cov4_x, k=kkk)
    knndist4 <- knndist4[ ,kkk]
    roc_obj <- roc(labs, knndist4, direction="<")
    cov4_knn[kk,i] <- roc_obj$auc
    #----------------------------------------------------
    
    #----------------------------------------------------
    # iForest
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu, ntree=20)
    # evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    u1_ifor[kk, i] <- roc_obj$auc
    
    
    # Test DOBIN - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc, ntree=20)
    # evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    dob_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x, ntree=20)
    # evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    pca_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca2 - isolationForest
    pca_x2 <- as.data.frame(pca_x2)
    tr <-IsolationTrees(pca_x2, ntree=20)
    #evaluate anomaly score
    as <-AnomalyScore(pca_x2,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    pca2_ifor[kk, i] <- roc_obj$auc
    
    
    
    # Test COV4 - isolationForest
    cov4_x <- as.data.frame(cov4_x)
    tr <-IsolationTrees(cov4_x, ntree=20)
    #evaluate anomaly score
    as <-AnomalyScore(cov4_x,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    cov4_ifor[kk, i] <- roc_obj$auc
    
  }
}

lab <- c("Iteration", "AUC")
gg_title <- "Area under ROC - LOF "
mu <- seq(1, 10, by =1)
df_lof <- draw_five_curves(u1_lof, dob_lof, pca_lof, pca2_lof, cov4_lof,  mu, lab, gg_title, 7)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"


gg_title <- "Area under ROC - KNN"
df_knn <- draw_five_curves(u1_knn, dob_knn,  pca_knn,  pca2_knn, cov4_knn, mu, lab, gg_title, 7)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"


gg_title <- "Area under ROC - iForest"
df_ifor <- draw_five_curves(u1_ifor,  dob_ifor, pca_ifor, pca2_ifor, cov4_ifor, mu, lab, gg_title, 7)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"


gg_title <- "Experiment 4"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
dfexp4 <- df
# df <- dfexp4
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab("Iteration") + ylab(lab[2]) + facet_grid(.~Method) + scale_x_discrete(limits=1:10)  + ylim(0.35,0.95)  + theme_bw()


# -------------------------------------------------------------------------------
# EXP 5 - PARABOLA EXAMPLE -  WITH ONE MOVING INTO THE CENTRE
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# EXP 5 - PLOT IT
# -------------------------------------------------------------------------------
set.seed(321)

x <-rnorm(805)
y <-x^2 + rnorm(805, sd=0.5)

X <- data.frame(
  x1 = x,
  x2 = y,
  x3 = rnorm(805),
  x4 = rnorm(805),
  x5 = rnorm(805),
  x6 = rnorm(805),
  x7 = rnorm(805),
  x8 = rnorm(805),
  x9 = rnorm(805),
  x10 = rnorm(805),
  x11 = rnorm(805),
  x12 = rnorm(805),
  x13 = rnorm(805),
  x14 = rnorm(805),
  x15 = rnorm(805),
  x16 = rnorm(805)
)
labs <- as.factor(c(rep(0,800), rep(1,5)))
mu2 <-  (8-1)*0.5
z <- cbind(rnorm(5,0, sd=0.1), rnorm(5,mu2, sd=0.1))
X[801:805, 1:2] <- z
ggplot(X[ ,1:2], aes(x1,x2)) + geom_point(aes(color=labs)) + theme_bw()




# -------------------------------------------------------------------------------
# EXP 5 - EXPERIMENT
# -------------------------------------------------------------------------------


u1_lof <- u1_knn <- u1_ifor <- matrix(0, nrow=pp, ncol=10)
pca_lof <- pca_knn <- pca_ifor <- matrix(0, nrow=pp, ncol=10)
pca2_lof <- pca2_knn <- pca2_ifor <- matrix(0, nrow=pp, ncol=10)
cov4_lof <- cov4_knn <- cov4_ifor <- matrix(0, nrow=pp, ncol=10)
dob_lof <- dob_knn <- dob_ifor <- matrix(0, nrow=pp, ncol=10)

hh <- 8
kkk <- min(max( ceiling(805/50),2), 25)
set.seed(123)
for(kk in 1:pp){
  
  x <-rnorm(805)
  y <-x^2 + rnorm(805, sd=0.5)
  
  X <- data.frame(
    x1 = x,
    x2 = y,
    x3 = rnorm(805),
    x4 = rnorm(805),
    x5 = rnorm(805),
    x6 = rnorm(805),
    x7 = rnorm(805),
    x8 = rnorm(805),
    x9 = rnorm(805),
    x10 = rnorm(805),
    x11 = rnorm(805),
    x12 = rnorm(805),
    x13 = rnorm(805),
    x14 = rnorm(805),
    x15 = rnorm(805),
    x16 = rnorm(805)
  )
  labs <- c(rep(0,800), rep(1,5))
  
  
  
  for(i in 1:10){
    mu2 <-  (i-1)*0.5
    z <- cbind(rnorm(5,0, sd=0.1), rnorm(5,mu2, sd=0.1))
    
    X[801:805, 1:2] <- z
    #----------------------------------------------------
    # COORDINATES 
    
    # DOBIN
    out <- dobin(X, frac=0.95, norm=1)
    XSc <- out$coords[ ,1:hh]
    
    # PCA
    pca <- prcomp(X, scale = TRUE, center = TRUE )
    pca_x <- pca$x[ ,1:hh]
    pca_x2 <- pca$x[ ,(hh+1):(2*hh)]
    
    
    # ICS package - cov4
    icscov4 <- cov4(X)
    cov4_x <- solve(icscov4)%*%t(X)
    cov4_x <- t(cov4_x[1:hh,])
    #----------------------------------------------------
    
    #----------------------------------------------------
    # LOF
    
    # Test Min-Max Scaling - LOF
    Xu <- apply(X, 2, unitize_1)
    lof1 <- lofactor(Xu, k=kkk)
    roc_obj1 <- roc(labs, lof1, direction="<")
    u1_lof[kk,i] <- roc_obj1$auc
    
    # Test DOBIN - LOF
    lof2 <- lofactor(XSc, k=kkk)
    roc_obj2 <- roc(labs, lof2, direction="<")
    dob_lof[kk,i] <- roc_obj2$auc
    
    
    # Test PCA - LOF
    lof3 <- lofactor(pca_x, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca_lof[kk,i] <- roc_obj$auc
    
    
    # Test PCA2 - LOF
    lof3 <- lofactor(pca_x2, k=kkk)
    roc_obj <- roc(labs, lof3, direction="<")
    pca2_lof[kk,i] <- roc_obj$auc
    
    
    # Test COV4 - LOF
    lof4 <- lofactor(cov4_x, k=kkk)
    roc_obj <- roc(labs, lof4, direction="<")
    cov4_lof[kk,i] <- roc_obj$auc
    #----------------------------------------------------
    
    
    #----------------------------------------------------
    # KNN
    
    # Test Min-Max Scaling - KNN
    knndist1 <- FNN::knn.dist(Xu, k=kkk)
    knndist1 <- knndist1[,kkk]
    roc_obj3 <- roc(labs, knndist1, direction="<")
    u1_knn[kk,i] <- roc_obj3$auc
    
    
    # Test DOBIN - KNN
    knndist2 <- FNN::knn.dist(XSc, k=kkk)
    knndist2 <- knndist2[,kkk]
    roc_obj4 <- roc(labs, knndist2, direction="<")
    dob_knn[kk,i] <- roc_obj4$auc
    
    
    # Test PCA - KNN
    knndist3 <- FNN::knn.dist(pca_x, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca_knn[kk,i] <- roc_obj$auc
    
    
    # Test PCA2 - KNN
    knndist3 <- FNN::knn.dist(pca_x2, k=kkk)
    knndist3 <- knndist3[ ,kkk]
    roc_obj <- roc(labs, knndist3, direction="<")
    pca2_knn[kk,i] <- roc_obj$auc
    
    
    # Test COV4 - KNN
    knndist4 <- FNN::knn.dist(cov4_x, k=kkk)
    knndist4 <- knndist4[ ,kkk]
    roc_obj <- roc(labs, knndist4, direction="<")
    cov4_knn[kk,i] <- roc_obj$auc
    #----------------------------------------------------
    
    #----------------------------------------------------
    # iForest
    
    # Test Min-Max Scaling - isolationForest
    Xu <- as.data.frame(Xu)
    tr <-IsolationTrees(Xu)
    # evaluate anomaly score
    as <-AnomalyScore(Xu,tr)
    roc_obj <- roc(labs, as$outF, direction="<")
    u1_ifor[kk, i] <- roc_obj$auc
    
    
    # Test DOBIN - isolationForest
    XSc <- as.data.frame(XSc)
    tr <-IsolationTrees(XSc)
    # evaluate anomaly score
    as <-AnomalyScore(XSc,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    dob_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca - isolationForest
    pca_x <- as.data.frame(pca_x)
    tr <-IsolationTrees(pca_x)
    # evaluate anomaly score
    as <-AnomalyScore(pca_x,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    pca_ifor[kk, i] <- roc_obj$auc
    
    
    # Test pca2 - isolationForest
    pca_x2 <- as.data.frame(pca_x2)
    tr <-IsolationTrees(pca_x2)
    # evaluate anomaly score
    as <-AnomalyScore(pca_x2,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    pca2_ifor[kk, i] <- roc_obj$auc
    
    
    # Test COV4 - isolationForest
    cov4_x <- as.data.frame(cov4_x)
    tr <-IsolationTrees(cov4_x)
    # evaluate anomaly score
    as <-AnomalyScore(cov4_x,tr)
    roc_obj <- roc(labs, as$outF, direction="<") 
    cov4_ifor[kk, i] <- roc_obj$auc
    
  }
}

lab <- c("mu", "AUC")
gg_title <- "Area under ROC - LOF "
mu <- seq(1, 10, by =1)
df_lof <- draw_five_curves(u1_lof, dob_lof, pca_lof, pca2_lof, cov4_lof,  mu, lab, gg_title, 8)
df_lof2 <- cbind.data.frame(df_lof, rep("LOF", dim(df_lof)[1]))
colnames(df_lof2)[4] <- "Method"

gg_title <- "Area under ROC - KNN"
df_knn <- draw_five_curves(u1_knn, dob_knn,  pca_knn, pca2_knn, cov4_knn, mu, lab, gg_title, 8)
df_knn2 <- cbind.data.frame(df_knn, rep("KNN", dim(df_lof)[1]))
colnames(df_knn2)[4] <- "Method"

gg_title <- "Area under ROC - iForest"
df_ifor <- draw_five_curves(u1_ifor,  dob_ifor, pca_ifor,  pca2_ifor, cov4_ifor, mu, lab, gg_title, 8)
df_ifor2 <- cbind.data.frame(df_ifor, rep("IFOREST", dim(df_lof)[1]))
colnames(df_ifor2)[4] <- "Method"


gg_title <- "Experiment 5"
df <- rbind.data.frame(df_lof2, df_knn2, df_ifor2)
dfexp5 <- df
ggplot(df, aes(xvals, AUC)) + geom_point(aes(color=Coordinates, shape=Coordinates), size=2.5) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab("Iteration") + ylab(lab[2]) + facet_grid(.~Method) + scale_x_discrete(limits=1:10)  + ylim(0.35, 0.95)  + theme_bw()






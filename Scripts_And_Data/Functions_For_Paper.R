draw_five_curves <- function(roc_u1, roc_sc1, roc_sc2, roc_sc3, roc_sc4,  xvals, lab, gg_title,  n){
  min_max <- apply(roc_u1, 2, mean)
  
  # dobin  ROC values
  dobin <- apply(roc_sc1, 2, mean)
  
  pca <- apply(roc_sc2, 2, mean)
  
  pca2 <- apply(roc_sc3, 2, mean)
  
  cov4 <- apply(roc_sc4, 2, mean)
  
  df <- cbind.data.frame(xvals, min_max, dobin, pca,  pca2,  cov4)
  
  df2 <- tidyr::gather(df, key, value, -xvals)
  
  colnames(df2)[3] <- "AUC"
  colnames(df2)[2] <- "Coordinates"
  df2[df2[ ,2]=="min_max" ,2] <- "All Vars"
  df2[df2[ ,2]=="dobin" ,2] <- paste(n, "DOBIN comps")
  df2[df2[ ,2]=="pca" ,2] <- paste(n, "First PC comps")
  df2[df2[ ,2]=="pca2" ,2] <- paste(n, "Last PC comps")
  df2[df2[ ,2]=="cov4" ,2] <- paste(n, "COV4 comps")
  
  if(lab[1]=="mu"){
    print(ggplot(df2, aes(xvals, AUC)) + geom_point(aes(color=Coordinates)) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab(expression(mu)) + ylab(lab[2])   + theme_bw()   )
  }else{
    print(ggplot(df2, aes(xvals, AUC)) + geom_point(aes(color=Coordinates)) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab(lab[1]) + ylab(lab[2])  + theme_bw()   )
    
  }
  return(df2)
}


draw_nemenyi_plots <- function(dat_all, method = 1, sd = 123){
  # dat_all is the merged data set
  
  
  na_sums <- apply(dat_all, 2, function(x)sum(is.na(x)))
  print(na_sums)
  print("Removing NA rows . . . ")
  na_sums <- apply(dat_all, 1, function(x)sum(is.na(x)))
  inds <- which(na_sums >0 )
  
  dat_all <- dat_all[-inds, ]
  
  
  if(method==1){ # LOF
    title <- "LOF"
    colset <- c(2:6)
  }else if(method==2){ # KNN
    title <- "KNN"
    colset <- c(7:11)
  }else if(method==3){ # iforest
    title <- "iForest"
    colset <- c(12:16)
  }
  
  file_source <- GetFileSources(dat_all[ ,1])
  
  dat2 <- dat_all[ ,colset]
  
  colnames(dat2) <- c("All Vars", "1/2 DOBIN  ", "First 1/2 PCA", "Last 1/2 PCA", "1/2 COV4")
  
  df <- cbind.data.frame(file_source, dat2)
  
  set.seed(sd)
  num_srcs <- length(unique(file_source))
  recs <- c()
  for(kk in 1:num_srcs){
    inds <- which(file_source ==unique(file_source)[kk])
    sam <- sample(inds, min(20, length(inds)))
    recs <- c(recs, sam)
  }  
  
  dfs <- df[recs, ]
  
  dfs[ ,-1] <- -1*dfs[ ,-1]
  
  nemenyi <- tsutils::nemenyi(as.matrix(dfs[ ,-1]), conf.level=0.95, sort=TRUE,  plottype="none", main="Nemenyi test average ranks") # sort=TRUE,
  order_methods <- names(nemenyi$means)
  df4 <- dfs[ ,-1]
  df5 <- df4[ ,order_methods]
  print("Doing the ordered graph!")
  nemenyi <- tsutils::nemenyi(as.matrix(df5), conf.level=0.95, sort=TRUE,  plottype="matrix", main=paste("Nemenyi test average ranks - " , title) )
  
  out <- list()
  out$friedman <- nemenyi$fpavl
  out$nemenyi <- nemenyi
  out$len2 <- dim(dfs)[1]
  return(out)
}


unitize_1 <- function(z) {
  # min-max normalization - 0 -1
  min.z <- min(z)
  max.z <- max(z)
  if ((max.z-min.z)==0)
    return(z)
  (z - min.z )/(max.z-min.z)
}


draw_box_plots <- function( dat_all, method=1 ){
  if(method==1){ # LOF
    title <- "LOF"
    colset <- 2:6
  }else if(method==2){ # KNN
    title <- "KNN"
    colset <- 7:11
  }else if(method==3){ # iforest
    title <- "iForest"
    colset <- 12:16
  }


  df <- dat_all[, colset]
  colnames(df) <- c("All Vars", "1/2 DOBIN", "1/2 PCA", "Last 1/2 PCA", "COV4")
  df2 <- cbind.data.frame((df[ ,2] - df[ ,1]), (df[ ,2] - df[ ,3]) , (df[ ,2] - df[ ,4]),(df[ ,2] - df[ ,5]))
  colnames(df2) <- c("1/2 DOBIN - All Vars", "1/2 DOBIN - First 1/2 PCA",  "1/2 DOBIN - Last 1/2 PCA", "1/2 DOBIN - 1/2 COV4" )
  df3 <- tidyr::gather(df2, key, value)
  colnames(df3)[1] <- "Coordinates"
  print( ggplot2::ggplot(df3, aes(x=Coordinates, y=value, fill=Coordinates)) + ggplot2::geom_boxplot() + ggplot2::ggtitle(title) + ggplot2::theme_bw() + xlab("Coordinates") + ylab("AUC Difference") ) # + coord_fixed(ratio=2)

}


GetFileSources <- function(filenames){
  file_source <-c()
  for(ll in 1:length(filenames)){
    fname <- filenames[ll]
    regobj1 <- regexpr("_C", fname)
    regobj2 <- regexpr("_withoutdupl", fname)
    if(regobj1[1]<0){
      regobj <- regobj2
    }else if(regobj2[1]<0){
      regobj <- regobj1
    }else{
      regobj <- regobj1
    }
    end.ind <- regobj[1]-1
    file_source <- c(file_source, substring(fname, 1, end.ind))
  }
  return(file_source)
}


# TAKEN FROM https://robjhyndman.com/hyndsight/crop-r-figures/
savepdf <- function(file, width=16, height=10)
{
  #folder <- paste( getwd(), "/paper/", sep="")
  folder <- "C:/Users/Sevvandi/Documents/repos/ssoutliers_private/paper/"
  fname <- paste(folder, file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}


# ------------------------------------------------------------------------
# COV3 COMPUTATIONS - BEGIN - TAKEN FROM ICS PACKAGE 
# ------------------------------------------------------------------------
cov3 <- function(X, location = "Mean", na.action = na.fail){
  X <- na.action(X)
  X.matrix <- as.matrix(X)
  p <- dim(X)[2]
  if (p < 2) 
    stop("'X' must be at least bivariate")
  if (is.numeric(location)) {
    if (length(location) != p) 
      stop("'location' is of wrong dimension")
    X.matrix <- sweep(X.matrix, 2, location)
    location = "Origin"
  }
  loc <- match.arg(location, c("Mean", "Origin"))
  if (loc == "Mean") 
    V <- .cov3moments.mean(X.matrix)
  if (loc == "Origin") 
    V <- .cov3moments.origin(X.matrix)
  return(V)
}


.cov3moments.mean<-function(X){
  p<-dim(X)[2]
  n<-dim(X)[1]
  data.centered<-sweep(X,2,colMeans(X),"-")
  Sigma.data.sqrt<-mat.sqrt(cov(X)) 
  radius<-sqrt(rowSums((data.centered %*% solve(Sigma.data.sqrt))^2))
  # INSERT SQRT - BEGIN
  radius <- sqrt(radius)
  # INSERT SQRT - END
  y<-radius*data.centered
  V<-(1/(n*(p+2)))*crossprod(y) 
  return(V) 
}

### covariance matrix based on 4th moments wrt to origin
### subroutine of cov4
###

.cov3moments.origin<-function(X){
  p<-dim(X)[2]
  n<-dim(X)[1]
  Sigma.data.sqrt<-mat.sqrt(ICS::covOrigin(X)) 
  radius<-sqrt(rowSums((X %*% solve(Sigma.data.sqrt))^2))
  V<-(1/(p+2))*ICS::covOrigin(radius*X)  
  return(V) 
}

mat.sqrt <- function(A){
  eigen.A<-eigen(A)
  sqrt.A<-eigen.A$vectors %*% tcrossprod(diag(eigen.A$values^0.5), eigen.A$vectors)
  return(sqrt.A)
}

# ------------------------------------------------------------------------
# COV3 COMPUTATIONS - END 
# ------------------------------------------------------------------------

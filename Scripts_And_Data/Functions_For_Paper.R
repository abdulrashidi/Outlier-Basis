draw_three_curves <- function(roc_u1, roc_sc1, roc_sc2,  xvals, lab, gg_title,  n){
  min_max <- apply(roc_u1, 2, mean)
  
  # dobin LOF ROC values
  dobin <- apply(roc_sc1, 2, mean)
  
  pca <- apply(roc_sc2, 2, mean)
  
  df <- cbind.data.frame(xvals, min_max, dobin, pca)
  
  df2 <- melt(df, id="xvals")
  colnames(df2)[3] <- "AUC"
  colnames(df2)[2] <- "Coordinates"
  levels(df2[ ,2])[levels(df2[ ,2])=="min_max"] <- "All Vars"
  levels(df2[ ,2])[levels(df2[ ,2])=="dobin"] <- paste(n, "DOBIN comps")
  levels(df2[ ,2])[levels(df2[ ,2])=="pca"] <- paste(n, "PC comps")
  if(lab[1]=="mu"){
    print(ggplot(df2, aes(xvals, AUC)) + geom_point(aes(color=Coordinates)) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab(expression(mu)) + ylab(lab[2])   + theme_bw()   )
  }else{
    print(ggplot(df2, aes(xvals, AUC)) + geom_point(aes(color=Coordinates)) + geom_line(aes(color=Coordinates)) + ggtitle(gg_title) + xlab(lab[1]) + ylab(lab[2])  + theme_bw()   )
    
  }
  return(df2)
}


unitize_1 <- function(z) {
  # min-max normalization - 0 -1
  min.z <- min(z)
  max.z <- max(z)
  if ((max.z-min.z)==0)
    return(z)
  (z - min.z )/(max.z-min.z)
}


draw_nemenyi_plots <- function(dat_all, method=2){
  # dat_all is the merged data set
  dat_all[is.na(dat_all)] <- 0
  
  if(method==1){ # LOF
    title <- "LOF"
    colset <- 2:4 
  }else if(method==2){ # KNN
    title <- "KNN"
    colset <- 5:7 
  }else if(method==3){ # iforest
    title <- "iForest"
    colset <- 8:10 
  }
  
  file_source <- GetFileSources(dat_all[ ,1])
  
  dat2 <- dat_all[ ,colset]
  colnames(dat2) <- c("All Vars", "1/2 DOBIN  ", "1/2 PCA")
  
  df <- cbind.data.frame(file_source, dat2)
  
  # --- Make a big data frame with source, outlier method and sensitivity to normalization
  for(j in 1:3){
    temp <- reshape::melt(df[ ,c(1,(j+1))] )
    colnames(temp) <- c("source", "method", "value")
    if(j==1){
      dat <- temp
    }else{
      dat <- rbind.data.frame(dat, temp)
    }
  }
  
  df2 <- stats::aggregate(dat$value, by=list(m=dat$method, s=dat$source), FUN=median)
  friedman_test <- stats::friedman.test(y=df2$x, groups=df2$m, blocks=df2$s)
  
  df[ ,-1] <- -1*df[ ,-1] 
  df3 <- stats::aggregate(df[,-1], by=list(file_source), FUN=median)
  nemenyi <- tsutils::nemenyi(as.matrix(df3[ ,-1]), conf.level=0.95, sort=TRUE,  plottype="none", main="Nemenyi test average ranks") # sort=TRUE,
  order_methods <- names(nemenyi$means)
  df4 <- df3[ ,-1]
  df5 <- df4[ ,order_methods]
  print("Doing the ordered graph!")
  nemenyi <- tsutils::nemenyi(as.matrix(df5), conf.level=0.95, sort=TRUE,  plottype="matrix", main=paste("Nemenyi test average ranks - " , title))
  
  out <- list()
  out$friedman <- friedman_test
  out$nemenyi <- nemenyi
  out$len2 <- dim(df3)[1]
  return(out)
}


draw_density_plots <- function( dat_all, method=1 ){
  if(method==1){ # LOF
    title <- "LOF"
    colset <- 2:4
  }else if(method==2){ # KNN
    title <- "KNN"
    colset <- 5:7
  }else if(method==3){ # iforest
    title <- "iForest"
    colset <- 8:10
  }
  
  
  df <- dat_all[, colset]
  colnames(df) <- c("All Variables", "1/2 DOBIN comps", "1/2 PCA comps")
  df3 <- reshape::melt(df)
  colnames(df3)[1] <- "Coordinates"
  print( ggplot2::ggplot(df3, aes(x=value, color=Coordinates)) + ggplot2::geom_density(alpha=0.1, size=1) + ggplot2::ggtitle(title) + ggplot2::theme_bw() + xlab("AUC") )
  
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

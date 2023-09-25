#---------------------------------------------------------------------------------------------------
#----------------------------------- function for PCA ----------------------------------------------
f.pca = function (x){
  ca <- match.call()
  if (ncol(x) > nrow(x)){
    u = princomp(t(x))
    u$call = ca
    return(u)
  }
  
  mu = rowMeans(x)
  xb <- x - mu
  xb.svd <- svd(xb)
  pc <- t(xb) %*% xb.svd$u
  dimnames(pc)[[2]] <- paste("PC", 1:ncol(pc), sep = "")
  loading <- xb.svd$u
  dimnames(loading) <- list(paste("V", 1:nrow(loading), sep = ""), 
                            paste("Comp.", 1:ncol(loading), sep = ""))
  class(loading) <- "loadings"
  sd = xb.svd$d/sqrt(ncol(x))
  names(sd) <- paste("Comp.", 1:length(sd), sep = "")
  pc <- list(sdev = sd, loadings = loading, center = mu, 
             scale = rep(1, length(mu)), n.obs = ncol(x), scores = pc, call = ca)
  class(pc) <- "princomp"
  return(pc)
}

#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
#----------------------------------- Intermediate PCA ----------------------------------------------

IntermediatePCA<-function(Micro.mat,
                          Prognostic,
                          Survival,
                          Censor,
                          index){
  if (is.matrix(Micro.mat)) {
    pc1 <- f.pca(as.matrix(Micro.mat[ ,index]))[[6]][,1]
  } else {
    pc1<-Micro.mat[index]
  }

  if (is.null(Prognostic)) {

    cdata <- data.frame(Survival=Survival[index], Censor=Censor[index], pc1)
    m0 <- survival::coxph(Surv(Survival, Censor==1) ~ pc1, data=cdata)
  }

  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nPrgFac<-ncol(Prognostic)
      NameProg<-colnames(Prognostic)
      cdata <- data.frame(Survival[index], Censor[index], pc1, Prognostic[index,])
      colnames(cdata) = c("Survival", "Censor", "pc1", NameProg)
      eval(parse(text=paste( "m0 <-survival::coxph(Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {
      stop(" Argument 'Prognostic' is NOT a data frame ")
    }

  }

  return(list(m0=m0,pc1=pc1,cdata=cdata))
}


#---------------------------------------------------------------------------------------------------
#----------------------------------- Intermediate PLS ----------------------------------------------


IntermediatePLS<-function(Micro.mat,
                          Prognostic,
                          Survival,
                          Censor,
                          index){
  if (is.matrix(Micro.mat)) {

    DataPLS<-data.frame(1:length(index))
    DataPLS$g<-as.matrix(t(Micro.mat[ ,index]))
    colnames(DataPLS)[1]<-c("Survival")
    DataPLS[,1]<-Survival[index]
    plsr.1 <- pls::plsr(Survival ~ g, method="simpls", ncomp = 2, scale =TRUE,data = DataPLS, validation =  "CV")
    pc1<-pls::scores(plsr.1)[,1] # extract the first com
  } else {
    pc1<-Micro.mat[index]
  }

  if (is.null(Prognostic)) {
    cdata <- data.frame(Survival=Survival[index], Censor=Censor[index], pc1)
    m0 <- survival::coxph(Surv(Survival, Censor==1) ~ pc1, data=cdata)
  }
  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nPrgFac<-ncol(Prognostic)
      NameProg<-colnames(Prognostic)
      cdata <- data.frame(Survival[index], Censor[index], pc1, Prognostic[index,])
      colnames(cdata) = c("Survival", "Censor", "pc1", NameProg)
      eval(parse(text=paste( "m0 <-survival::coxph(Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {
      stop(" Argument 'Prognostic' is NOT a data frame ")
    }
  }
  return(list(m0=m0,pc1=pc1,cdata=cdata))
}

#'Frequency of Selected Taxa from the LASSO, Elastic-net Cross-Validation
#'
#' The function selects the frequency of selection from the shrinkage method (LASSO, Elastic-net) based on cross validation,
#' that is the number of times each taxon occur during the cross-validation process.

#' This function outputs the mostly selected taxa during the LASSO and Elastic-net cross validation.
#' Selected top taxa are ranked based on frequency of selection and also a particular frequency can be selected.
#' In addition, it visualizes the selected top taxa based on the minimum frequency specified.
#' @param Object An object of class \code{\link[MicrobiomeSurv]{cvle}} returned from the function \code{\link[MicrobiomeSurv]{CVLasoelacox}}.
#' @param TopK The number of Top K taxa (5 by default) to be displayed in the frequency of selection graph.
#' @return A vector of taxa and their frequency of selection. Also, a graphical representation is displayed.
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{cvmm}}, \code{\link[survival]{coxph}},
#' \code{\link[MicrobiomeSurv]{EstimateHR}},\code{\link[MicrobiomeSurv]{CVLasoelacox}}
#' @examples
#' \donttest{
#' # Cross-Validation for LASSO and ELASTIC-NET
#' CV_lasso_fam_shan_w3 = CVLasoelascox(Survival = survival_data_w3$Survival,
#'                                      Censor = survival_data_w3$Censor,
#'                                      Micro.mat = fam_shan_trim_w3,
#'                                      Prognostic = prog_fam_w3,
#'                                      Standardize = TRUE,
#'                                      Alpha = 1,
#'                                      Fold = 4,
#'                                      Ncv = 100,
#'                                      nlambda = 100)
#'
#'
#' # Using the function
#' MiFreq_fam_shan_w3 = MiFreq(Object = CV_lasso_fam_shan_w3, TopK=5)
#' }
#' @export MiFreq

MiFreq<-function(Object,
                 TopK=20){

  #Decrease=FALSE
  if (class(Object)!="cvle") stop("Invalid object class.")

  MFreq <- Object@mi.mat
  fr<-colSums(MFreq)
  names(fr)<-rownames(Object@Micro.mat)

  if(!is.null(TopK)){
    top <- sort(fr,decreasing=TRUE)
    topn <- top[1:TopK]
    barplot(topn,las=2,col=rainbow(length(topn)), ylim = c(0,max(topn) + 1),  ylab="",cex.names=0.6,main=paste( "Top ", TopK, " most selected taxa ",sep=""),cex.lab=1,cex.main=1.5  )
    return(top)
  } else
  {
    #top.fr<-fr[fr==N]
    top.fr <- sort(fr,decreasing=TRUE)
    topnn <- top.fr[1:TopK]
    barplot(top.fr,las=2,col=rainbow(length(topnn)), ylim = c(0,max(topnn)),  ylab="",cex.names=0.6,main=paste( "Taxa ", "with selection frequency=", N ,sep=""),cex.lab=1,cex.main=1.5  )
    return(fr)
  }
}

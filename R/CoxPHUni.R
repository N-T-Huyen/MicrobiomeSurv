
#' This function will fit the full and reduced models and calculate LRT raw p-value and adjusted p-value based on BH Method

#'@param Survival The time to event outcome.
#'@param Censor An indicator variable indicate the subject is censored or not.
#'@param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#'@param Micro.mat a microbiome matrix, can be at otu, family or any level of the ecosystem. Rows are taxa, columns are subjectsc.
#'@param Method A multiplicity adjustment Method that user can choose. The default is BH Method.

#'@return A relative abundance matrix of OTUs
#' \item{coef}{coefficient of one microbiome (OTU or family, ...)}
#' \item{exp.coef}{exponential of the coefficient}
#' \item{p.value.LRT}{raw LRT p-value}
#' \item{p.value}{adjusted p-value based on chosen Method}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' \code{\link[MicrobiomeSurv]{CoxPHUni}}

#' @import utils
#' @import stats
#' @import Biobase


#' @export CoxPHUni



CoxPHUni = function(Survival, Censor, Prognostic, Micro.mat, Method = "BH"){
  n.mi = dim(Micro.mat)[1]
  n.obs = dim(Micro.mat)[2]
  coef = exp.coef = p.value.LRT = c(1 : n.mi)
  if (is.data.frame(Prognostic)) {
    nPrgFac<-ncol(Prognostic)
    NameProg<-colnames(Prognostic)
  }

  if (dim(Prognostic)[2] == 1){
    cox.prog = eval(parse(text = paste("survival::coxph(survival::Surv(Survival, Censor)", paste(" ~ Prognostic[ ,1])"))))
  } else{
    cox.prog = eval(parse(text = paste("survival::coxph(survival::Surv(Survival, Censor) ~ NameProg[1]", paste("+", NameProg[2:nPrgFac], sep="", collapse =""), paste(")"))))
  }

  for(i in 1 : n.mi){
    xi = Micro.mat[i, ]
    datai <- data.frame(Survival, Censor, xi, Prognostic)
    modeli = eval(parse(text = paste("survival::coxph(survival::Surv(Survival, Censor) ~ xi", paste("+", NameProg[1:nPrgFac], sep="", collapse =""), ",data=datai)" , sep="" )))
    coef[i] = round(summary(modeli)$coefficients[1,1], 4)
    exp.coef[i] = round(summary(modeli)$coefficients[1,2], 4)
    p.value.LRT[i] = round(lmtest::lrtest(cox.prog, modeli)[2,5], 4)
  }

  p.value = round(p.adjust(p.value.LRT, method = Method, n = length(p.value.LRT)), 4)
  summary = cbind(coef, exp.coef, p.value.LRT, p.value)
  rownames(summary) = rownames(Micro.mat)
  colnames(summary) = c("coef", "exp.coef", "p.value.LRT", "p.value")
  return(summary)
}

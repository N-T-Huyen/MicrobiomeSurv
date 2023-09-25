#'Cross validation for sequentially increases taxa
#'
#' This function does cross validation for the taxon by taxon analysis while sequentially increasing the number of taxa as specified.
#'
#' The function is a cross validation version of the function \code{\link[MicrobiomeSurv]{SITaxa}}.
#' This function firstly processes the cross validation for the taxon by taxon analysis results, and then sequentially considers top k taxa.
#' The function recompute first PCA or PLS on train data and estimate risk scores on both test and train data only on the microbiome matrix with top k taxa.
#' Patients are then classified as having low or high risk based on the test data where the cutoff used is mean of the risk score.
#' The process is repeated for each top K taxa sets.
#' @param Object An object of class \code{\link[MicrobiomeSurv]{cvmm}}.
#' @param Top The Top k number of taxa to be used.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Censor A vector of censoring indicator.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @return A object of class \code{\link[MicrobiomeSurv]{cvsit}} is returned with the following values
#'    \item{HRpca}{A 3-way array in which first, second, and third dimensions correspond to number of taxa, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PCA.}
#'    \item{HRpls}{A 3-way array in which first, second, and third dimensions correspond to number of taxa, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PLS.}
#'    \item{Ntaxa}{The number of taxa in the reduced matrix.}
#'   \item{Ncv}{The number of cross validation done.}
#'   \item{Top}{A sequence of top k taxa considered. Default is Top = seq(5, 100, by=5)}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{MSpecificCoxPh}}, \code{\link[MicrobiomeSurv]{SITaxa}}
#' @examples
#' \donttest{
#' # Prepare data
#' Week3_response = read_excel("Week3_response.xlsx")
#' Week3_response = data.frame(Week3_response)
#' Week3_response = Week3_response[order(Week3_response$SampleID), ]
#' Week3_response$Treatment_new = ifelse(Week3_response$Treatment=="3PATCON",0,1)
#' surv_fam_shan_w3 = data.frame(cbind(as.numeric(Week3_response$T1Dweek),
#' as.numeric(Week3_response$T1D)))
#' colnames(surv_fam_shan_w3) = c("Survival", "Censor")
#' prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
#' colnames(prog_fam_shan_w3) = c("Treatment")
#' fam_shan_trim_w3 = read_excel("fam_shan_trim_w3.xlsx")
#' names_fam_shan_trim_w3 = c(fam_shan_trim_w3[ ,1])$X.
#' fam_shan_trim_w3 = data.matrix(fam_shan_trim_w3[ ,2:82])
#' rownames(fam_shan_trim_w3) = names_fam_shan_trim_w3

#' # Getting the cvmm object
#' CVCox_taxon_fam_shan_w3 = CVMSpecificCoxPh(Fold=3,
#'                                            Survival = survival_data_w3$Survival,
#'                                            Micro.mat = fam_shan_trim_w3,
#'                                            Censor = survival_data_w3$Censor,
#'                                            Reduce=TRUE,
#'                                            Select=5,
#'                                            Prognostic=prog_fam_w3,
#'                                            Mean = TRUE,
#'                                            Ncv=100)
#'
#' # Using the function
#'  CVSITaxa_fam_shan_w3 = CVSITaxa(Object = CVCox_taxon_fam_shan_w3,
#'                                  Top=seq(1, 6, by=2),
#'                                  Survival = survival_data_w3$Survival,
#'                                  Censor = survival_data_w3$Censor,
#'                                  Prognostic=prog_fam_w3)
#'
#' # Get the class of the object
#' class(CVSITaxa_fam_shan_w3)     # An "cvsit" Class
#' }
#'
#' @export CVSITaxa

CVSITaxa = function(Object,
                     Top=seq(5, 100, by=5),
                     Survival,
                     Censor,
                     Prognostic=NULL){

   Decrease=FALSE
  if (class(Object)!="cvmm") stop("Invalid class object.")
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  mi.mat=matrix(0, Object@Ncv,Object@n.mi)

  Micro.mat=Object@Rdata
  mi.names = rownames(Micro.mat)

  n.mi=Object@n.mi
  n.patients=ncol(Micro.mat)


  if (Object@n.mi<max(Top)) stop("The max(Top) should be less than or equal to total number of taxa")


  HRPC=array(NA,dim=c(Object@Ncv,3,length(Top)))
  HRPL=array(NA,dim=c(Object@Ncv,3,length(Top)))

  #i=1
  options( warn = -1)
  for (j in 1: length(Top)){
    message('Analysis for ',paste0("Top",Top[j]))
    for (i in 1:Object@Ncv){
      message('Cross validation loop ',j)
      hrsetTE=Object@HRTest[,c(1,3,4),i]
      hrsetTr=Object@HRTrain[,c(1,3,4),i]

      Names.KTaxaT=rownames(Micro.mat)[order(hrsetTE[,1],decreasing=Decrease)]
      index.Top.KTaxaT = order(hrsetTE[,1],decreasing=Decrease)
      index.Top.KTaxa=index.Top.KTaxaT[1:Top[j]]

      taxoni =Micro.mat[intersect(rownames(Micro.mat),Names.KTaxaT[1:Top[j]]),]

      #PCA ----------------

      TrainTemp=IntermediatePCA(taxoni, Prognostic, Survival, Censor, Object@train[i, ])
      TestTemp= IntermediatePCA(taxoni, Prognostic, Survival, Censor, Object@test[i, ])
      m1 = TrainTemp$m0
      TrtandMi1=summary(m1)[[7]][c("pc1"),1]
      rm(m1)
      p1.train     = TrtandMi1[1]*TrainTemp$pc1
      p1.test     =  TrtandMi1[1]*TestTemp$pc1

      TemptaxoniTE =EstimateHR(p1.test,
                                Data.Survival=data.frame(Survival=Survival[Object@test[i,]], Censor=Censor[Object@test[i,]]),
                                Prognostic=data.frame(Prognostic[Object@test[i,],]),
                                Plots = FALSE,
                                Mean = TRUE,
                                Quantile = 0.5)

      HRPC[i, ,j]= (summary(TemptaxoniTE$SurvResult)[[8]][1,])[-2]

      #PLS -------------------
      TrainTemp1=IntermediatePLS(taxoni, Prognostic, Survival, Censor, Object@train[i,])
      TestTemp1 =IntermediatePLS(taxoni, Prognostic, Survival, Censor, Object@test[i,])
      m1 = TrainTemp1$m0
      TrtandMi2=summary(m1)[[7]][c("pc1"),1]
      rm(m1)
      p1.train     = TrtandMi2[1]*TrainTemp1$pc1
      p1.test     = TrtandMi2[1]*TestTemp1$pc1
      TemptaxoniTE =EstimateHR(p1.test,
                                Data.Survival=data.frame(Survival=Survival[Object@test[i,]],Censor=Censor[Object@test[i,]]),
                                Prognostic=data.frame(Prognostic[Object@test[i,],]),
                                Plots = FALSE,
                                Mean = TRUE,
                                Quantile = 0.5)

    HRPL[i, ,j]= (summary( TemptaxoniTE$SurvResult)[[8]][1,])[-2]

    }

  }

  Ncv=Object@Ncv
  n.mi=Object@n.mi
  return(new("cvsit",HRpca=HRPC,HRpls=HRPL,Ntaxa=n.mi,Ncv=Ncv,Top=Top))
  }

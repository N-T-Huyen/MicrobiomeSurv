% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CVPcaPls.R
\name{CVPcaPls}
\alias{CVPcaPls}
\title{Cross Validations for PCA and PLS based methods}
\usage{
CVPcaPls(
  Fold = 3,
  Survival,
  Micro.mat,
  Censor,
  Reduce = TRUE,
  Select = 15,
  Prognostic = NULL,
  Ncv = 5,
  DR = "PCA"
)
}
\arguments{
\item{Fold}{Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to test the results.}

\item{Survival}{A vector of survival time with length equals to number of subjects.}

\item{Micro.mat}{A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.}

\item{Censor}{A vector of censoring indicator.}

\item{Reduce}{A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.}

\item{Select}{Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.}

\item{Prognostic}{A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.}

\item{Ncv}{The Number of cross validation loop. Default is 100.}

\item{DR}{The dimension reduction method. It can be either "PCA" for Principle components analysis or "PLS" for Partial least squares.}
}
\value{
A object of class \code{\link[MicrobiomeSurv]{cvpp}} is returned with the following values
  \item{Result}{A dataframe containg the estimated Hazard ratio of the test dataset and the training dataset.}
  \item{Ncv}{The number of cross validation performed.}
  \item{Method}{The dimesion reduction method used.}
  \item{CVtrain}{The training dataset indices matrix used for the cross validation.}
  \item{CVtest}{The test dataset indices matrix used for the cross validation.}
  \item{Select}{The number of taxa used for the dimesion reduction method used.}
}
\description{
This function does cross validation for the analysis performs by \code{\link[MicrobiomeSurv]{SurvPcaClass}}
and \code{\link[MicrobiomeSurv]{SurvPlsClass}} functions where the dimension reduction methods can either be PCA and PLS.
}
\details{
This function does cross validation for the analysis using two reduction method. The reduction method can be PCA or PLS.
If it is PCA then the \code{\link[MicrobiomeSurv]{SurvPcaClass}} is internally used for the cross validation
and \code{\link[MicrobiomeSurv]{SurvPlsClass}} otherwise.
}
\examples{
\donttest{
# Prepare data
data(Week3_response)
Week3_response = data.frame(Week3_response)
surv_fam_shan_w3 = data.frame(cbind(as.numeric(Week3_response$T1Dweek),
as.numeric(Week3_response$T1D)))
colnames(surv_fam_shan_w3) = c("Survival", "Censor")
prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
colnames(prog_fam_shan_w3) = c("Treatment")
data(fam_shan_trim_w3)
names_fam_shan_trim_w3 =
c("Unknown", "Lachnospiraceae", "S24.7", "Lactobacillaceae", "Enterobacteriaceae", "Rikenellaceae")
fam_shan_trim_w3 = data.matrix(fam_shan_trim_w3[ ,2:82])
rownames(fam_shan_trim_w3) = names_fam_shan_trim_w3
# Using the function
CVPls_fam_shan_w3 = CVPcaPls(Fold = 3,
                            Survival = surv_fam_shan_w3$Survival,
                            Micro.mat = fam_shan_trim_w3,
                            Censor = surv_fam_shan_w3$Censor,
                            Reduce=TRUE,
                            Select=5,
                            Prognostic = prog_fam_shan_w3,
                            Ncv=10,
                            DR = "PLS")

# Get the class of the object
class(CVPls_fam_shan_w3)     # An "cvpp" Class

# Method that can be used for the result
show(CVPls_fam_shan_w3)
summary(CVPls_fam_shan_w3)
plot(CVPls_fam_shan_w3)
}
}
\seealso{
\code{\link[MicrobiomeSurv]{SurvPlsClass}},
\code{\link[MicrobiomeSurv]{SurvPcaClass}}
}
\author{
Thi Huyen Nguyen, \email{thihuyen.nguyen@uhasselt.be}

Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@gsk.com}

Ziv Shkedy
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSpecificCoxPh.R
\name{MSpecificCoxPh}
\alias{MSpecificCoxPh}
\title{Taxon by taxon Cox proportional analysis}
\usage{
MSpecificCoxPh(
  Survival,
  Micro.mat,
  Censor,
  Reduce = FALSE,
  Select = 5,
  Prognostic = NULL,
  Mean = TRUE,
  Quantile = 0.5,
  Method = "BH"
)
}
\arguments{
\item{Survival}{A vector of survival time with length equals to number of subjects}

\item{Micro.mat}{A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of subjects.}

\item{Censor}{A vector of censoring indicator.}

\item{Reduce}{A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach.}

\item{Select}{Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.}

\item{Prognostic}{A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.}

\item{Mean}{The cut off value for the classifier, default is the mean cutoff.}

\item{Quantile}{If users want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they wish to use. The default is the median cutoff.}

\item{Method}{Multiplicity adjustment methods.}
}
\value{
A object of class \code{\link[MicrobiomeSurv]{ms}} is returned with the following values
  \item{Result}{The cox proportional regression result for each taxon}
  \item{HRRG}{The hazard ratio statistics (Hazard-ratio, Lower confidence interval and upper confidence interval) of the riskgroup based on the riskscore and the cut off value for each taxon}
  \item{Group}{The classification of the subjects based on each taxon analysis}
  \item{Mi.names}{The names of the taxa for the analysis}
}
\description{
The Function fits cox proportional hazard model and does classification for each taxon separately
}
\details{
This function fits  taxon by taxon Cox proportional hazard model and perform the classification based on a microbiome risk score which has been estimated using a single taxon.
Function is useful for majority vote classification method and taxon by taxon analysis and also for top K taxa.
}
\examples{
\donttest{
# Prepare data
data(Week3_response)
Week3_response = data.frame(Week3_response)
surv_fam_shan_w3 =
data.frame(cbind(as.numeric(Week3_response$T1Dweek), as.numeric(Week3_response$T1D)))
colnames(surv_fam_shan_w3) = c("Survival", "Censor")
prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
colnames(prog_fam_shan_w3) = c("Treatment")
data(fam_shan_trim_w3)
names_fam_shan_trim_w3 =
c("Unknown", "Lachnospiraceae", "S24.7", "Lactobacillaceae", "Enterobacteriaceae", "Rikenellaceae")
fam_shan_trim_w3 = data.matrix(fam_shan_trim_w3[ ,2:82])
rownames(fam_shan_trim_w3) = names_fam_shan_trim_w3
# Using the function
Cox_taxon_fam_shan_w3 = MSpecificCoxPh(Survival = surv_fam_shan_w3$Survival,
                                      Micro.mat = fam_shan_trim_w3,
                                      Censor = surv_fam_shan_w3$Censor,
                                      Reduce=FALSE,
                                      Select=5,
                                      Prognostic = prog_fam_shan_w3,
                                      Mean = TRUE,
                                      Method = "BH")

# Results
show(Cox_taxon_fam_shan_w3)
summary(Cox_taxon_fam_shan_w3)
}
}
\seealso{
\code{\link[survival]{coxph}},  \code{\link[MicrobiomeSurv]{EstimateHR}}
}
\author{
Thi Huyen Nguyen, \email{thihuyen.nguyen@uhasselt.be}

Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@gsk.com}

Ziv Shkedy
}

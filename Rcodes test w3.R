rm(list = ls())


library(survival)
library(survminer)
library(lubridate)
library("Biobase")
library("grid")
library(ggplot2)
library("nlme")
library(stringr)
library(dplyr)
library(phyloseq)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(IntegratedJM)
library(lmtest)
library(kableExtra)
library(readxl)
library(writexl)
library(loo)
library(reshape2)
library(glmnet)
library(MicrobiomeSurv)



# **************** #
# Data preparation #
# **************** #

set.seed(1012387)

# Taxonomic information of otus
metadata_taxonomy <- read.csv("metadata_taxonomy.csv")

otus_data <- lapply(str_split(metadata_taxonomy$Taxon, ";"), FUN = function(k) {

  y <- sapply(str_split(k, "__"), FUN = "[", 2) # take only column 2
  y[y == ""] <- NA
  return(y)

}) %>% do.call(rbind.data.frame, .)

names(otus_data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# ********* #
# OTU level #
# ********* #

# otu table week 3

otus_data <- otus_data %>% dplyr::mutate(OTUID = metadata_taxonomy$OTUID)
otus_data$Family = replace_na(as.vector(otus_data$Family), "Unknown")

Week3_otu = read_excel("Week3_otu.xlsx")
Week3_otu = data.frame(Week3_otu)
dim(Week3_otu)

# Preparing data for analysis at OTU level
otu_mat_w3 = t(data.matrix(Week3_otu[ , 1:2720]))
dim(otu_mat_w3)
colnames(otu_mat_w3) = Week3_otu$SampleID

## Filtering first step

otu_w3 = FirstFilter(Micro.mat = otu_mat_w3)
dim(otu_w3) #420  81

## 2nd stage filtering of OTUs based on Percentage of zero counts
Week3_otu$Treatment_new = ifelse(Week3_otu$Treatment == "3PATCON", 0, 1)
table(Week3_otu$Treatment_new)

n_obs = dim(otu_w3)[2]
n_control = table(Week3_otu$Treatment_new)[1]
n_treated = table(Week3_otu$Treatment_new)[2]
n_otu = dim(otu_mat_w3)[1]

### Calculate zeros per groups
zero_per_group_otu_w3 = ZerosPerGroup(Micro.mat = otu_mat_w3,
                                      groups = Week3_otu$Treatment_new, week = 3,
                                      n.obs = n_obs, n.control = n_control,
                                      n.treated = n_treated, n.mi = n_otu,
                                      plot = TRUE)


### Removing OTUs having at least 70% of zeros in at least one treatment group
otu_trim_w3 = SecondFilter(zero.per.group = zero_per_group_otu_w3[[1]],
                           Micro.mat = otu_w3, threshold = 0.7, week = 3)

dim(otu_trim_w3) # 7 81

## Convert absolute abundance to relative abundance
ra_otu_trim_w3 = GetRA(Micro.mat = otu_trim_w3)
dim(ra_otu_trim_w3)

## Plot of percentage of zeros per groups
windows()
zero_per_group_otu_w3[[2]]
savePlot(filename = "Percentage of zeros otu week 3", type = "png")



# ************ #
# Family level #
# ************ #

# Preparing data for analysis at family level
fam = otus_data[c("Family", "OTUID")]
otus$OTUID = row.names(otus)
otufam = inner_join(otus, fam)
otu_fam_w3 = otufam[,which(colnames(otufam) %in% c(Week3_response$SampleID, "OTUID", "Family"))]
fam_w3 = otu_fam_w3[ , c("OTUID", "Family")]
otu_data_w3 = otu_fam_w3[ , 1:n_obs]
rownames(otu_data_w3) = otu_fam_w3$OTUID

fam_shan_w3 = SummaryData(Micro.mat = otu_data_w3, info = fam_w3, measure = "shannon")

fam_shan_filter_w3 = FirstFilter(t(fam_shan_w3[, -c(dim(fam_shan_w3)[2])]))
n_fam_w3 = dim(fam_shan_filter_w3)[1]

### Calculate zeros per groups
zero_per_group_fam_shan_w3 = ZerosPerGroup(Micro.mat = fam_shan_filter_w3,
                                           groups = Week3_otu$Treatment_new, week = 3,
                                           n.obs = n_obs, n.control = n_control,
                                           n.treated = n_treated, n.mi = n_fam_w3,
                                           plot = TRUE)

### Removing OTUs having at least 90% of zeros in at least one treatment group
#### 70% -- 0 families
#### 80% -- 2 families
#### 90% -- 6 families

fam_shan_trim_w3 = SecondFilter(zero.per.group = zero_per_group_fam_shan_w3[[1]],
                                Micro.mat = fam_shan_filter_w3, threshold = 0.9, week = 3)

dim(fam_shan_trim_w3) # 6 81

# Analysis at family level

n_fam_shan_w3 = dim(fam_shan_trim_w3)[1]
Week3_response = read_excel("Week3_response.xlsx")
Week3_response = data.frame(Week3_response)
Week3_response = Week3_response[order(Week3_response$SampleID), ]
identical(colnames(fam_shan_trim_w3), Week3_response$SampleID)
Week3_response$Treatment_new = ifelse(Week3_response$Treatment=="3PATCON",0,1)


## The univariate analysis
surv_fam_shan_w3 = data.frame(cbind(as.numeric(Week3_response$T1Dweek),
                                    as.numeric(Week3_response$T1D)))
colnames(surv_fam_shan_w3) = c("Survival", "Censor")

prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
colnames(prog_fam_shan_w3) = c("Treatment")
table(prog_fam_shan_w3[, 1])

summary_fam_shan_w3 = CoxPHUni(Survival = surv_fam_shan_w3$Survival, Censor = surv_fam_shan_w3$Censor,
                               Prognostic = prog_fam_shan_w3, Micro.mat = fam_shan_trim_w3, Method = "BH")
summary_fam_shan_w3

## Analysis of the taxon having smallest p-value

top1_fam_shan_w3 = Top1Uni(Result = summary_fam_shan_w3, Micro.mat = fam_shan_trim_w3,
                           Survival = surv_fam_shan_w3$Survival, Censor = surv_fam_shan_w3$Censor,
                           Plots = TRUE)

### Plot the observed Kaplan-Meier curves per group
windows()
top1_fam_shan_w3[[3]]
savePlot(filename = "Kaplan-Meier curves per group fam shan top 1 week 3", type = "png")

### log-rank test
top1_fam_shan_w3[[4]]



## LASSO without MCCV

set.seed(1234)
windows()
lasso_fam_shan_w3 = Lasoelascox(Survival = surv_fam_shan_w3$Survival,
                                Censor = surv_fam_shan_w3$Censor,
                                Micro.mat = fam_shan_trim_w3,
                                Prognostic = prog_fam_shan_w3,
                                Plots = TRUE,
                                Standardize = TRUE,
                                Alpha = 1,
                                Fold = 4,
                                nlambda = 100,
                                Mean = TRUE)

savePlot(filename = "LASSO fam shan week 3", type = "png")


### Non zeros coefficients
lasso_fam_shan_w3$Coefficients.NonZero

### Selected microbiome
selected_mi_fam_shan_w3 = lasso_fam_shan_w3$Selected.mi

### summary

summary(lasso_fam_shan_w3$SurvFit)

### Microbiome risk score
mi_score_fam_shan_w3 = cbind(Score = lasso_fam_shan_w3$Risk.Scores, prog_fam_shan_w3)
mean(lasso_fam_shan_w3$Risk.Scores)

windows()
ggplot(mi_score_fam_shan_w3, aes(x = Score, color = as.factor(Treatment))) +
  geom_density()+ylab("Density")+
  geom_vline(xintercept = mean(lasso_fam_shan_w3$Risk.Scores), linetype="dotted")+
  guides(color = guide_legend(title = "Treatmet"))
savePlot(filename = "Density of microbiome score shannon w3", type = "png")


windows()
den_lasso = density(mi_score_fam_shan_w3[ ,1])
plot(den_lasso, frame = FALSE, col = "blue",main = " ",
     xlab = "Score")
abline(v=mean(mi_score_fam_shan_w3[ ,1]), lty = "dashed")
savePlot(filename = "Density of microbiome score LASSO shannon w3", type = "png")


### Plot the observed Kaplan-Meier curves per group and a boxplot per group

windows()
est_HR_fam_shan_w3 = EstimateHR(Risk.Scores = lasso_fam_shan_w3$Risk.Scores,
                                Data.Survival = lasso_fam_shan_w3$Data.Survival,
                                Prognostic = prog_fam_shan_w3, Plots = TRUE,
                                Mean = TRUE)

### Estimate the HR
est_HR_fam_shan_w3$SurvResult

### Kaplan Meier plot

est_HR_fam_shan_w3$KMplot
savePlot(filename = "Kaplan Meier LASSO fam shan w3", type = "png")

### Boxplot of low and high risks
windows()
est_HR_fam_shan_w3$SurvBPlot
savePlot(filename = "Boxplot LASSO fam shan w3", type = "png")


## MCCV for LASSO model

set.seed(7890)
CV_lasso_fam_shan_w3 = CVLasoelascox(Survival = surv_fam_shan_w3$Survival,
                                     Censor = surv_fam_shan_w3$Censor,
                                     Micro.mat = fam_shan_trim_w3,
                                     Prognostic = prog_fam_shan_w3,
                                     Standardize = TRUE,
                                     Alpha = 1,
                                     Fold = 4,
                                     Ncv = 100,
                                     nlambda = 100)

windows()
plot(CV_lasso_fam_shan_w3, ylim = c(0, 4))
savePlot(filename = "HR LASSO fam shan w3", type = "png")


show(CV_lasso_fam_shan_w3)
summary(CV_lasso_fam_shan_w3)

## MCCV for LASSO model with a grid of alpha
grid.alpha <- seq(0.1, 0.9, by = 0.1)

alpha_fun = function(Alpha){CVLasoelascox(Survival = surv_fam_shan_w3$Survival,
                                          Censor = surv_fam_shan_w3$Censor,
                                          Micro.mat = fam_shan_trim_w3,
                                          Prognostic = prog_fam_shan_w3,
                                          Standardize = TRUE,
                                          Alpha = Alpha,
                                          Fold = 4,
                                          Ncv = 10,
                                          nlambda = 100,
                                          Mean = TRUE)}



CV_lasso_fam_shan_w3_grid = sapply(grid.alpha,
                                   FUN = function(Alpha){CVLasoelascox(Survival = surv_fam_shan_w3$Survival,
                                                                       Censor = surv_fam_shan_w3$Censor,
                                                                       Micro.mat = fam_shan_trim_w3,
                                                                       Prognostic = prog_fam_shan_w3,
                                                                       Standardize = TRUE,
                                                                       Alpha = Alpha,
                                                                       Fold = 4,
                                                                       Ncv = 100,
                                                                       nlambda = 100,
                                                                       Mean = TRUE)})

str(CV_lasso_fam_shan_w3_grid)
CV_lasso_fam_shan_w3_grid[1]
str(CV_lasso_fam_shan_w3_grid[1])
slot(CV_lasso_fam_shan_w3_grid[[1]], "HRTest")[ ,1]




HRTest = sapply(1:length(grid.alpha),
                function(k) slot(CV_lasso_fam_shan_w3_grid[[k]], "HRTest")[, 1])


Ncv = 100
windows()
boxplot(HRTest,
        names = grid.alpha[1:length(grid.alpha)],
        xlab = expression(alpha), ylab = "HR estimate",
        main = paste("HR on test set", Ncv, "runs"),
        ylim = c(0, 5), #max(HRTest)
        col = 2:(length(grid.alpha)+1),
        cex.main = 1.5,
        cex.lab = 1.3)
savePlot(filename = paste("HR on test set", ncv, "runs LASSO fam shan w3"), type = "png")




windows()
boxplot(HRTest,
        names = grid.alpha[1:length(grid.alpha)],
        xlab = expression(alpha), ylab = "HR estimate",
        main = "HR on Test Set 100 runs",
        ylim = c(0, 5), #max(HRTest)
        col = 2:(length(grid.alpha)+1),
        cex.main = 1.5,
        cex.lab = 1.3)
savePlot(filename = "HR on Test Set 100 runs LASSO fam shan w3", type = "png")




## Supervised principal component analysis without MCCV
set.seed(7890)

SPCA_fam_shan_w3 = SurvPcaClass(Survival = surv_fam_shan_w3$Survival,
                                Micro.mat = fam_shan_trim_w3,
                                Censor = surv_fam_shan_w3$Censor,
                                Reduce=TRUE,
                                Select=5,
                                Prognostic = prog_fam_shan_w3,
                                Plots = TRUE,
                                Mean = TRUE)

rownames(SPCA_fam_shan_w3$ReduMicro.mat)


windows()
SPCA_fam_shan_w3$KMplot
savePlot(filename = "Kaplan-Meier PcaClass fam shan week 3", type = "png")

windows()
SPCA_fam_shan_w3$SurvBPlot
savePlot(filename = "Boxplot PcaClass fam shan week 3", type = "png")

summary(SPCA_fam_shan_w3$SurvFit)[[8]][,-2]
SPCA_fam_shan_w3$SurvFit

### microbiome risk score
mi_score_SPCA_fam_shan_w3 = round(SPCA_fam_shan_w3$Riskscores, 4)
mean(mi_score_SPCA_fam_shan_w3)
sum(mi_score_SPCA_fam_shan_w3)

windows()
den = density(mi_score_SPCA_fam_shan_w3)
plot(den, frame = FALSE, col = "blue",main = " ",
     xlab = "Score")
abline(v=mean(mi_score_SPCA_fam_shan_w3), lty = "dashed")
savePlot(filename = "Density of microbiome score SPCA shannon w3", type = "png")



## Supervised partial least squares without MCCV

set.seed(7890)
SPLS_fam_shan_w3 = SurvPlsClass(Survival = surv_fam_shan_w3$Survival,
                                Micro.mat = fam_shan_trim_w3,
                                Censor = surv_fam_shan_w3$Censor,
                                Reduce=TRUE,
                                Select=5,
                                Prognostic = prog_fam_shan_w3,
                                Plots = TRUE,
                                Mean = TRUE)

windows()
SPLS_fam_shan_w3$KMplot
SPLS_fam_shan_w3$SurvBPlot
summary(SPLS_fam_shan_w3$SurvFit)


## MCCV for PCA and PLS
set.seed(7890)
CVPca_fam_shan_w3 = CVPcaPls(Fold = 3,
                             Survival = surv_fam_shan_w3$Survival,
                             Micro.mat = fam_shan_trim_w3,
                             Censor = surv_fam_shan_w3$Censor,
                             Reduce=TRUE,
                             Select=5,
                             Prognostic = prog_fam_shan_w3,
                             Ncv=100,
                             DR = "PCA")

windows()
plot(CVPca_fam_shan_w3)
summary(CVPca_fam_shan_w3)
show(CVPca_fam_shan_w3)

CVPls_fam_shan_w3 = CVPcaPls(Fold = 3,
                             Survival = surv_fam_shan_w3$Survival,
                             Micro.mat = fam_shan_trim_w3,
                             Censor = surv_fam_shan_w3$Censor,
                             Reduce=TRUE,
                             Select=5,
                             Prognostic = prog_fam_shan_w3,
                             Ncv=100,
                             DR = "PLS")

windows()
plot(CVPls_fam_shan_w3)
savePlot(filename = "HR PLS fam shan week 3", type = "png")

summary(CVPls_fam_shan_w3)
show(CVPls_fam_shan_w3)



## Majority voting without MCCV

set.seed(7890)
Cox_taxon_fam_shan_w3 = MSpecificCoxPh(Survival = surv_fam_shan_w3$Survival,
                                       Micro.mat = fam_shan_trim_w3,
                                       Censor = surv_fam_shan_w3$Censor,
                                       Reduce=TRUE,
                                       Select=5,
                                       Prognostic = prog_fam_shan_w3,
                                       Mean = TRUE,
                                       Method = "BH")

show(Cox_taxon_fam_shan_w3)
summary(Cox_taxon_fam_shan_w3)

windows()
plot(Cox_taxon_fam_shan_w3)
savePlot(filename = "Feature specific fam shan week 3", type = "png")

CVCox_taxon_fam_shan_w3 = CVMSpecificCoxPh(Fold=3,
                                           Survival = surv_fam_shan_w3$Survival,
                                           Micro.mat = fam_shan_trim_w3,
                                           Censor = surv_fam_shan_w3$Censor,
                                           Reduce=TRUE,
                                           Select=5,
                                           Prognostic=prog_fam_shan_w3,
                                           Mean = TRUE,
                                           Ncv=100)


show(CVCox_taxon_fam_shan_w3)

windows()
plot(CVCox_taxon_fam_shan_w3)
savePlot(filename = "HR training and test sets first taxon fam shan week 3", type = "png")

summary(CVCox_taxon_fam_shan_w3)

windows()
plot(CVCox_taxon_fam_shan_w3, which = 5)
savePlot(filename = "HR training and test sets taxon 5 fam shan week 3", type = "png")

windows()
Majority_fam_shan_w3 = Majorityvotes(Result = Cox_taxon_fam_shan_w3,
                                     Prognostic = prog_fam_shan_w3,
                                     Survival = surv_fam_shan_w3$Survival,
                                     Censor = surv_fam_shan_w3$Censor,
                                     J=1)

savePlot(filename = "Majority votes fam shan week 3", type = "png")


## MCCV for Majority voting
set.seed(7890)
CVMajority_fam_shan_w3 = CVMajorityvotes(Survival = surv_fam_shan_w3$Survival,
                                         Micro.mat = fam_shan_trim_w3,
                                         Censor = surv_fam_shan_w3$Censor,
                                         Reduce=TRUE,
                                         Select=5,
                                         Mean = TRUE,
                                         Prognostic = prog_fam_shan_w3,
                                         Fold=3,
                                         Ncv=100)
show(CVMajority_fam_shan_w3)
summary(CVMajority_fam_shan_w3)
windows()
plot(CVMajority_fam_shan_w3)
savePlot(filename = "CVMajorityvotes fam shan week 3", type = "png")



## Sequential Increase in Taxa for the PCA or PLS classifier
set.seed(7890)
windows()
SITaxa_fam_shan_w3 = SITaxa(TopK=5,
                            Survival = surv_fam_shan_w3$Survival,
                            Micro.mat = fam_shan_trim_w3,
                            Censor = surv_fam_shan_w3$Censor,
                            Reduce=TRUE,
                            Select=5,
                            Prognostic=prog_fam_shan_w3,
                            Plot = TRUE,
                            DM="PLS")
SITaxa_fam_shan_w3
savePlot(filename = "HR with CI PLS for top 5 fam shan week 3", type = "png")



## Cross validation for sequentially increases taxa
set.seed(7890)
CVSITaxa_fam_shan_w3 = CVSITaxa(Object = CVCox_taxon_fam_shan_w3,
                                Top=seq(1, 6, by=2),
                                Survival = surv_fam_shan_w3$Survival,
                                Censor = surv_fam_shan_w3$Censor,
                                Prognostic=prog_fam_shan_w3)
show(CVSITaxa_fam_shan_w3)
summary(CVSITaxa_fam_shan_w3)

windows()
plot(CVSITaxa_fam_shan_w3, type = 2)
savePlot(filename = "PLS fam shan week 3", type = "png")


windows()
plot(CVSITaxa_fam_shan_w3, type = 1)
savePlot(filename = "PCA fam shan week 3", type = "png")



## Frequency of Selected Taxa from the LASSO, Elastic-net Cross-Validation
set.seed(7890)
windows()
MiFreq_fam_shan_w3 = MiFreq(Object = CV_lasso_fam_shan_w3,
                            TopK=5)
savePlot(filename = "MiFreq fam shan week 3", type = "png")


windows()
MiFreq_fam_shan_w3 = MiFreq(Object = CV_lasso_fam_shan_w3,
                            TopK = 6)
savePlot(filename = "MiFreq fam shan week 3", type = "png")

## Quantile sensitivity analysis


QuantileAnalysis_fam_shan_w3 = QuantileAnalysis(Survival = surv_fam_shan_w3$Survival,
                                                Micro.mat = fam_shan_trim_w3,
                                                Censor = surv_fam_shan_w3$Censor,
                                                Reduce=TRUE,
                                                Select= 5,
                                                Prognostic=prog_fam_shan_w3,
                                                Plots = FALSE,
                                                DM="SM",
                                                Alpha =1)
QuantileAnalysis_fam_shan_w3


set.seed(NULL)
windows()
QuantileAnalysis_PCA_fam_shan_w3 = QuantileAnalysis(Survival = surv_fam_shan_w3$Survival,
                                                    Micro.mat = fam_shan_trim_w3,
                                                    Censor = surv_fam_shan_w3$Censor,
                                                    Reduce=TRUE,
                                                    Select= 5,
                                                    Prognostic=prog_fam_shan_w3,
                                                    Plots = TRUE,
                                                    DM="PCA",
                                                    Alpha =1)

savePlot(filename = "HR and CI for the sensitivity analysis of the cutoff value fam shan week 3", type = "png")

QuantileAnalysis_PCA_fam_shan_w3




## Null Distribution of the Estimated HR

DistHR_fam_shan_w3 = DistHR(Survival = surv_fam_shan_w3$Survival,
                            Micro.mat = fam_shan_trim_w3,
                            Censor = surv_fam_shan_w3$Censor,
                            Prognostic=prog_fam_shan_w3,
                            Mean = TRUE,
                            Quantile=0.5,
                            Reduce= FALSE,
                            Select = 5,
                            nperm=100,
                            case=4,
                            Validation="PCAbased")
show(DistHR_fam_shan_w3)
summary(DistHR_fam_shan_w3)

windows()
plot(DistHR_fam_shan_w3)
savePlot(filename = "An approximation of the null distribution of the HR PCAbased fam shan week 3", type = "png")




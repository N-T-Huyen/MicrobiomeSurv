---
output:
  pdf_document: default
  html_document: default
---

# MicrobiomeSurv

<!-- badges: start -->
<!-- badges: end -->

R package : A biomarker validation approach for predicting time to event using microbiome data, this package develope biomarker detection for microbiome data. It contains a set of functions and cross validation methods to validate and select biomarkers when the outcome of interest is time to event. The package can handle prognostic factors and mainly microbiome matrix as input, the package can served as biomarker validation tool.

## Installation

You can install the development version of MicrobiomeSurv from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("N-T-Huyen/MicrobiomeSurv")
```


## Why use the package
* It can be used with any form of high dimensional/omics data such as: Microbiome data, Metabolic data, Gene expression matrix.
* It develops any level of the ecosystem from the high dimensional data to be used for other purpose.
* It also employs data reduction techniques such as PCA, PLS, Lasso and elastic net, and majority votes. 
* It classifies subjects into low and high risk group based on the microbiome scoro.
* It incorporates the use of subject prognostic information for the to enhance the biomarker for classification.
* It gives information about the time to develop a desease of subjects depending on the classification.



## A quick Demostration to solve a problem

``` r
"Problem of interest"

"Given a set of subjects with known risks cores and prognostic features how can we use this information to obtain their risk of develop T1D and what group does each respective subject belongs to?"

```

``` r
##  Loading the package
library("MicrobiomeSurv")

## Dataset
Dataset and the full codes example is published on github 

##  This function does classification, survival estimation and visualization
est_HR_fam_shan_w3 = EstimateHR(Risk.Scores = lasso_fam_shan_w3$Risk.Scores,
                                Data.Survival = lasso_fam_shan_w3$Data.Survival,
                                Prognostic = prog_fam_shan_w3, Plots = TRUE,
                                Mean = TRUE)

## Survival information
est_HR_fam_shan_w3$SurvResult


## Group information
est_HR_fam_shan_w3$Riskgroup
```

## Functions in the package

| Category	|	Functions	|	Description	                                |
| --------- | --------- | ------------------------------------------- |
| Basic	|	FirstFilter	|	Remove taxa with all zeros count.	
|				|	ZerosPerGroup	|	Obtain number and proportion of zeros per group.	
|				|	SecondFilter	|	Filter taxa using a pre-specified threshold.	
|				|	GetRA	|	Convert OTU count to relative abundance.	
|				|	SummaryData	|	Obtain microbiome measurements at higher level such as family, order, etc. 
|				|	CoxPHUni	|	Summary statistics for all taxa using univariate CoxPh model. BH method is used for multiplicity adjustment.	
|				|	Top1Uni	|	Do analysis for the taxon that has the smallest adjusted p-value from the CoxPh model. 
|				| MSpecificCoxPh	| Feature specific analysis. 
|				|	SurvPcaClass	|	Classifier based on first PCA.	
|				|	SurvPlsClass	|	Classifier based on first PLS. 	
|				|	Majorityvotes	|	Classifiction for Majority Votes.	
|				|	Lasoelascox	|	Wrapper function for glmnet.	
|				| MSData | Generate Artificial Microbiome Survival Data.
|				| EstimateHR | Classify subjects into low or high risk groups. Visualize survival fit along with HR estimates.
|				Advance	|	CVLasoelascox	|	Cross-Validations for Lasso Elastic Net predictive models and Classification.
|				|	CVSITaxa	|	Cross-Validation for Top $K_{1}, \ldots, K_{n}$ taxa.	
|				|	CVPcaPls	|	Cross-Validations for PCA and PLS based methods.	
|				|	CVMajorityvotes	|	Cross-Validation for majority votes.	
|				|	MiFreq	|	Frequency of selected taxa from Cross-Validation.	
|				|	QuantileAnalysis	|	Sensitivity of the quantile used for classification.	
|				|	DistHR	|	Null distribution of the estimated HR.	
|				|	SITaxa	|	Sequentially increase the number of top $K$ taxa.


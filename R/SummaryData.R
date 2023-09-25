
#' This function gives indices such as Observed richness, Shannon index, Inverse Simpson, ... of higher level such as levelily, order, phylum, ...

#'@param Micro.mat an OTU matrix with OTUs in rows and subjects in columns.
#'@param info A n x 2 matrix containing a column of OTU's names and a column of the corresponding information of the chosen level.
#'@param measure The indices at chosen level that user wishes to use. It can be observed richness, Shannon index, inverse Simpson, ...
#'@return A matrix of the selected measurement of the chosen level.
#' \item{level.measure}{A matrix of measurements at levelily level of patients}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{SummaryData}}
#' @examples
#' \donttest{
#' # Preparing data for analysis at family level
#' fam = otus_data[c("Family", "OTUID")]
#' otus$OTUID = row.names(otus)
#' otufam = inner_join(otus, fam)
#' otu_fam_w3 = otufam[ ,which(colnames(otufam) %in% c(Week3_response$SampleID, "OTUID", "Family"))]
#' fam_w3 = otu_fam_w3[ ,c("OTUID", "Family")]
#' otu_data_w3 = otu_fam_w3[ , 1:n_obs]
#' rownames(otu_data_w3) = otu_fam_w3$OTUID
#'
#' # USing the function
#' fam_shan_w3 = SummaryData(Micro.mat = otu_data_w3, info = fam_w3, measure = "shannon")
#' }
#' @export SummaryData

SummaryData = function(Micro.mat, info, measure = "observed"){
  Micro.mat = data.frame(Micro.mat)
  Micro.mat$OTUID = rownames(Micro.mat)
  otu.level = merge(Micro.mat, info, by = "OTUID")
  nest = otu.level[, -c(1)]%>% group_by(across(dim(otu.level)[2] - 1)) %>% nest() # group by levelily or other level
  level = function(nest){
    lapply(1:length(nest$data), function(i) {
    level.measure = microbiome::alpha(data.matrix(nest$data[[i]]), index = measure)
    names(level.measure) = nest[[1]][i]
  return(level.measure)
    })
  }
  level.measure = data.frame(level(nest))
  level.measure$SampleID = rownames(level.measure)
  return(level.measure)
}


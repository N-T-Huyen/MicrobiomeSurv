
#' This function convert OTU matrix to RA matrix.

#'@param Micro.mat an OTU matrix with OTUs in rows and subjects in columns.
#'@return A relative abundance matrix of OTUs
#' \item{ra}{Relative abundance matrixs}

#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{GetRA}}

#' @import utils
#' @import stats
#' @import Biobase

#' @export GetRA



GetRA <- function(Micro.mat) {

  ra <- matrix(0,nrow(Micro.mat),ncol(Micro.mat))
  rownames(ra) <- rownames(Micro.mat)
  colnames(ra) <- colnames(Micro.mat)

  for(i in 1:ncol(Micro.mat)){

    if(sum(Micro.mat[,i]) == 0){
      ra[,i] <- 0
    } else {
      ra[,i] <- Micro.mat[,i]/sum(Micro.mat[,i])
    }

  }

  return(ra)
}

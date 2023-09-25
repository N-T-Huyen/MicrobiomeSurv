

ElasticNetVaryAlpha = function(Survival,
                               Censor,
                               Micro.mat,
                               Prognostic,
                               Standardize,
                               grid.alpha,
                               Fold,
                               Ncv,
                               nlambda,
                               Mean,
                               Quantile,
                               Plot = TRUE){
  CV_lasso_fam_shan_w3_grid = sapply(grid.alpha,
                                     FUN = function(Alpha){CVLasoelascox(Survival = Survival,
                                                                         Censor = Censor,
                                                                         Micro.mat = Micro.mat,
                                                                         Prognostic = Prognostic,
                                                                         Standardize = Standardize,
                                                                         Alpha = Alpha,
                                                                         Fold = Fold,
                                                                         Ncv = Ncv,
                                                                         nlambda = nlambda,
                                                                         Mean = Mean,
                                                                         Quantile = Quantile)})

  str(CV_lasso_fam_shan_w3_grid)
  CV_lasso_fam_shan_w3_grid[1]
  str(CV_lasso_fam_shan_w3_grid[1])
  slot(CV_lasso_fam_shan_w3_grid[[1]], "HRTest")[ ,1]




  HRTest = sapply(1:length(grid.alpha),
                  function(k) slot(CV_lasso_fam_shan_w3_grid[[k]], "HRTest")[, 1])

  if (Plot){
    windows()
    boxplot(HRTest,
            names = grid.alpha[1:length(grid.alpha)],
            xlab = expression(alpha), ylab = "HR estimate",
            main = paste("HR on test set", Ncv, "runs"),
            ylim = c(0, 5), #max(HRTest)
            col = 2:(length(grid.alpha)+1),
            cex.main = 1.5,
            cex.lab = 1.3)
  }

}

#' Performs the Hotelling T2 test in Fourier space.
#' Modified based on the 'graph.T2.test' function from R package 'DEGraph'.
#' Updates include: (1) updated the tolerance for detecting linear dependencies;
#' (2) added an argument "nmin", the minimum of sample size between the two groups.
#' @param X1 A n1 x p numeric matrix, observed data for class 1: p variables, n1 observations.
#' @param X2 A n2 x p numeric matrix, observed data for class 2: p variables, n2 observations.
#' @param G An object of class graphAM or graphNEL, the graph to be used in the two-sample test.
#' @param alpha0 significance level of the hypothesis test.
#' @param lfA A list returned by laplacianFromA(), containing the Laplacian eigen vectors and eigen values.
#' @param k A numeric value, number of Fourier components retained for the test.
#' @param nmin the minimum of sample size between the two groups.
#' @return A list with class "htest".
#' @export
graph.T2.test = function (X1, X2, G = NULL, lfA = NULL, k = ncol(X1),nmin){
  tol <- 1e-80 # update the the tolerance for detecting linear dependencies in the columns of a.
  p <- ncol(X1)
  U <- lfA$U
  egVal <- lfA$l
  kIdx <- (egVal <= max(egVal[k], tol))
  rk <- max(which(kIdx))
  rk = min(rk,nmin)
  ut <- T2.test(X1 %*% U[, 1:rk, drop = FALSE], X2 %*% U[, 1:rk, drop = FALSE])
}

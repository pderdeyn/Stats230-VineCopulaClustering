#' label points by vcmm
#'
#' @param x data
#' @param margin_pars marginal parameters
#' @param M1 vine structure for component 1
#' @param M2 vine structure for component 2
#' @param F1 copula families for component 1
#' @param F2 copula families for component 2
#' @param cop1 copula params for component 1
#' @param cop2 copula params for component 2
#' @param mix mixing parameters for components
#' @return labels
#' @export
fit_points_to_vcmm <- function(x,margin_pars,M1,M2,F1,F2,cop1,cop2,mix) {
  gam1<-margin_pars[,,1]
  gam2<-margin_pars[,,2]

  n<-dim(x)[1]

  u1 <- fGarch::psstd(x[,1], mean=gam1[1,1], sd=gam1[2,1], nu=gam1[3,1], xi=gam1[4,1])
  u2 <- pnorm(x[,2], mean=gam1[1,2], sd=gam1[2,2])
  u3 <- plogis(x[,3], location=gam1[1,3], scale=gam1[2,3])
  udata1<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly

  u1 <- fGarch::psstd(x[,1], mean=gam2[1,1], sd=gam2[2,1], nu=gam2[3,1], xi=gam2[4,1])
  u2 <- pnorm(x[,2], mean=gam2[1,2], sd=gam2[2,2])
  u3 <- plogis(x[,3], location=gam2[1,3], scale=gam2[2,3])
  udata2<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly
  r<-matrix(rep(0,2*n),ncol=2)

  RVM1 <- VineCopula::RVineMatrix(M1,F1,cop1)
  rvine_density1 <- VineCopula::RVinePDF(udata1,RVM1)
  m1 <- fGarch::dsstd(x[,1], mean=gam1[1,1], sd=gam1[2,1], nu=gam1[3,1], xi=gam1[4,1])
  m2 <- dnorm(x[,2], mean=gam1[1,2], sd=gam1[2,2])
  m3 <-  dlogis(x[,3], location=gam1[1,3], scale=gam1[2,3])
  margin_density1 <- m1 * m2 * m3

  RVM2 <- VineCopula::RVineMatrix(M2,F2,cop2)
  rvine_density2 <- VineCopula::RVinePDF(udata2,RVM2)
  m1 <- fGarch::dsstd(x[,1], mean=gam2[1,1], sd=gam2[2,1], nu=gam2[3,1], xi=gam2[4,1])
  m2 <- dnorm(x[,2], mean=gam2[1,2], sd=gam2[2,2])
  m3 <-  dlogis(x[,3], location=gam2[1,3], scale=gam2[2,3])
  margin_density2 <- m1 * m2 * m3

  denom <- mix[1]*margin_density1*rvine_density1+mix[2]*margin_density2*rvine_density2
  r[,1]<-mix[1]*margin_density1*rvine_density1
  r[,2]<-mix[2]*margin_density2*rvine_density2


  r[,1]<-r[,1]/denom
  r[,2]<-r[,2]/denom
  return(r)


}

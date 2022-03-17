#' internal function
#' @noRd
vcmm1_margins <- function(pars,x,gam,p,r,M,fam,cop) {
  if (p==1) { # skew t => 4 parameters
    gam[1,p]<-pars[1]
    gam[2,p]<-pars[2]
    gam[3,p]<-pars[3]
    gam[4,p]<-pars[4]
  }
  else if (p==2) { # normal => 2 params
    gam[1,p]<-pars[1]
    gam[2,p]<-exp(pars[2])
  }
  else { # skew normal => 3 params
    gam[1,p]<-pars[1]
    gam[2,p]<-pars[2]
  }
  RVM <- VineCopula::RVineMatrix(M,fam,cop)
  # TODO make sure x is indexed correctly
  u1 <- fGarch::psstd(x[,1], mean=gam[1,1], sd=gam[2,1], nu=gam[3,1], xi=gam[4,1])
  u2 <- pnorm(x[,2], mean=gam[1,2], sd=gam[2,2])
  u3 <- plogis(x[,3], location=gam[1,3], scale=gam[2,3])
  u<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly
  rvinedensity <- VineCopula::RVinePDF(u,RVM)

  m1 <- fGarch::dsstd(x[,1], mean=gam[1,1], sd=gam[2,1], nu=gam[3,1], xi=gam[4,1])
  m2 <- dnorm(x[,2], mean=gam[1,2], sd=gam[2,2])
  m3 <- dlogis(x[,3], location=gam[1,3], scale=gam[2,3])

  margin_density <- m1 * m2 * m3
  density <- rvinedensity * margin_density
  density[which(density == 0)] <- 1e-100
  -sum(r*log(density))
}

#' internal function
#' @noRd
vcmm2_margins <- function(pars,x,gam,p,r,M,fam,cop) {
  if (p==2) { # skew t => 4 parameters
    gam[1,p]<-pars[1]
    gam[2,p]<-pars[2]
    gam[3,p]<-pars[3]
    gam[4,p]<-pars[4]
  }
  else if (p==1) { # normal => 2 params
    gam[1,p]<-pars[1]
    gam[2,p]<-exp(pars[2])
  }
  else { # skew normal => 3 params
    gam[1,p]<-pars[1]
    gam[2,p]<-pars[2]
  }
  RVM <- VineCopula::RVineMatrix(M,fam,cop)
  # TODO make sure x is indexed correctly
  u2 <- fGarch::psstd(x[,2], mean=gam[1,2], sd=gam[2,2], nu=gam[3,2], xi=gam[4,2])
  u1 <- pnorm(x[,1], mean=gam[1,1], sd=gam[2,1])
  u3 <- plogis(x[,3], location=gam[1,3], scale=gam[2,3])
  u<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly
  rvinedensity <- VineCopula::RVinePDF(u,RVM)

  m2 <- fGarch::dsstd(x[,2], mean=gam[1,2], sd=gam[2,2], nu=gam[3,2], xi=gam[4,2])
  m1 <- dnorm(x[,1], mean=gam[1,1], sd=gam[2,1])
  m3 <- dlogis(x[,3], location=gam[1,3], scale=gam[2,3])

  margin_density <- m1 * m2 * m3
  density <- rvinedensity * margin_density
  density[which(density == 0)] <- 1e-100
  -sum(r*log(density))
}



#' fit model to toy example
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
#' @param iter iterations of ECM
#' @param maxit iterations for each CM step
#' @param tol tolerance for stopping by log likelihood
#' @return model
#' @export
test_vineclust <- function(x,margin_pars,M1,M2,F1,F2,cop1,cop2,mix,
                           iter=100,maxit=10,tol=0.01) {
  gam1<-margin_pars[,,1]
  gam2<-margin_pars[,,2]
  lls<-numeric(iter)


  n<-dim(x)[1]

  u1 <- fGarch::psstd(x[,1], mean=gam1[1,1], sd=gam1[2,1], nu=gam1[3,1], xi=gam1[4,1])
  u2 <- pnorm(x[,2], mean=gam1[1,2], sd=gam1[2,2])
  u3 <- plogis(x[,3], location=gam1[1,3], scale=gam1[2,3])
  udata1<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly

  u2 <- fGarch::psstd(x[,2], mean=gam2[1,2], sd=gam2[2,2], nu=gam2[3,2], xi=gam2[4,2])
  u1 <- pnorm(x[,1], mean=gam2[1,1], sd=gam2[2,1])
  u3 <- plogis(x[,3], location=gam2[1,3], scale=gam2[2,3])
  udata2<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly
  r<-matrix(rep(0,2*n),ncol=2)


  for (t in seq(iter)) {
    RVM1 <- VineCopula::RVineMatrix(M1,F1,cop1)
    rvine_density1 <- VineCopula::RVinePDF(udata1,RVM1)
    m1 <- fGarch::dsstd(x[,1], mean=gam1[1,1], sd=gam1[2,1], nu=gam1[3,1], xi=gam1[4,1])
    m2 <- dnorm(x[,2], mean=gam1[1,2], sd=gam1[2,2])
    m3 <-  dlogis(x[,3], location=gam1[1,3], scale=gam1[2,3])
    margin_density1 <- m1 * m2 * m3

    RVM2 <- VineCopula::RVineMatrix(M2,F2,cop2)
    rvine_density2 <- VineCopula::RVinePDF(udata2,RVM2)
    m2 <- fGarch::dsstd(x[,2], mean=gam2[1,2], sd=gam2[2,2], nu=gam2[3,2], xi=gam2[4,2])
    m1 <- dnorm(x[,1], mean=gam2[1,1], sd=gam2[2,1])
    m3 <-  dlogis(x[,3], location=gam2[1,3], scale=gam2[2,3])
    margin_density2 <- m1 * m2 * m3

    denom <- mix[1]*margin_density1*rvine_density1+mix[2]*margin_density2*rvine_density2
    r[,1]<-mix[1]*margin_density1*rvine_density1
    r[,2]<-mix[2]*margin_density2*rvine_density2

    ls<-r[,1]+r[,2]
    lls[t]<-sum(log(ls))
    if (t>1) {
      if (lls[t]-lls[t-1] < tol) {
        break
      }
    }


    r[,1]<-r[,1]/denom
    r[,2]<-r[,2]/denom

    mix[1]<-sum(r[,1])/n
    mix[2]<-sum(r[,2])/n

    # CM-margins component 1 node 1
    pars<-gam1[,1]
    vcmm1_margins(par=pars,x=x,gam=gam1,p=1,r=r[,1],M=M1,fam=F1,cop=cop1)
    opt_margins <- optim(par=pars, vcmm1_margins,
                         lower = c(min(x[,1]), 0.01*sd(x[,1]), 2.0001, 0.0001),
                         upper = c(max(x[,1]), 100*sd(x[,1]), 100, 100),
                         x=x,gam=gam1,p=1,r=r[,1],M=M1,fam=F1,cop=cop1,
                         method="L-BFGS-B",control=list(maxit=maxit))
    gam1[1,1]<-opt_margins$par[1]
    gam1[2,1]<-opt_margins$par[2]
    gam1[3,1]<-opt_margins$par[3]
    gam1[4,1]<-opt_margins$par[4]

    # CM-margins component 1 node 2
    pars<-gam1[,2]
    opt_margins <- optim(par=pars, vcmm1_margins,
                         x=x,gam=gam1,p=2,r=r[,1],M=M1,fam=F1,cop=cop1,
                         method="BFGS", control=list(maxit=maxit))
    gam1[1,2]<-opt_margins$par[1]
    gam1[2,2]<-exp(opt_margins$par[2])

    # CM-margins component 1 node 3
    pars<-gam1[,3]
    opt_margins <- optim(par=pars, vcmm1_margins,
                         lower = c(min(x[,3]), 0.01*sd(x[,3]), 0.0001),
                         upper = c(max(x[,3]), 100*sd(x[,3]), 100),
                         x=x,gam=gam1,p=3,r=r[,1],M=M1,fam=F1,cop=cop1,
                         method="L-BFGS-B",control=list(maxit=maxit))
    gam1[1,3]<-opt_margins$par[1]
    gam1[2,3]<-opt_margins$par[2]
    gam1[3,3]<-opt_margins$par[3]

    # component 2 margins


    vcmm2_margins(par=pars,x=x,gam=gam2,p=1,r=r[,2],M=M2,fam=F2,cop=cop2)
    # CM-margins component 2 node 1
    pars<-gam2[,2]
    opt_margins <- optim(par=pars, vcmm2_margins,
                         lower = c(min(x[,2]), 0.01*sd(x[,2]), 2.0001, 0.0001),
                         upper = c(max(x[,2]), 100*sd(x[,2]), 100, 100),
                         x=x,gam=gam2,p=2,r=r[,2],M=M2,fam=F2,cop=cop2,
                         method="L-BFGS-B",control=list(maxit=maxit))
    gam2[1,2]<-opt_margins$par[1]
    gam2[2,2]<-opt_margins$par[2]
    gam2[3,2]<-opt_margins$par[3]
    gam2[4,2]<-opt_margins$par[4]

    # CM-margins component 2 node 2
    pars<-gam2[,1]
    opt_margins <- optim(par=pars, vcmm2_margins,
                         x=x,gam=gam2,p=1,r=r[,2],M=M2,fam=F2,cop=cop2,
                         method="BFGS", control=list(maxit=maxit))
    gam2[1,1]<-opt_margins$par[1]
    gam2[2,1]<-exp(opt_margins$par[2])

    # CM-margins component 2 node 3
    pars<-gam2[,3]
    opt_margins <- optim(par=pars, vcmm2_margins,
                         lower = c(min(x[,3]), 0.01*sd(x[,3]), 0.0001),
                         upper = c(max(x[,3]), 100*sd(x[,3]), 100),
                         x=x,gam=gam2,p=3,r=r[,2],M=M2,fam=F2,cop=cop2,
                         method="L-BFGS-B",control=list(maxit=maxit))
    gam2[1,3]<-opt_margins$par[1]
    gam2[2,3]<-opt_margins$par[2]
    gam2[3,3]<-opt_margins$par[3]

    # component 1 copula pair parameters

    u1 <- fGarch::psstd(x[,1], mean=gam1[1,1], sd=gam1[2,1], nu=gam1[3,1], xi=gam1[4,1])
    u2 <- pnorm(x[,2], mean=gam1[1,2], sd=gam1[2,2])
    u3 <- plogis(x[,3], location=gam1[1,3], scale=gam1[2,3])
    udata1<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly
    RVM <- VineCopula::RVineMatrix(M1,F1,cop1)
    seq_RVM <-VineCopula::RVineSeqEst(udata1, RVM, weights=r[,1],  progress=FALSE)
    cop1<-seq_RVM$par

    # component 2 copula pair parameters

    u2 <- fGarch::psstd(x[,2], mean=gam2[1,2], sd=gam2[2,2], nu=gam2[3,2], xi=gam2[4,2])
    u1 <- pnorm(x[,1], mean=gam2[1,1], sd=gam2[2,1])
    u3 <- plogis(x[,3], location=gam2[1,3], scale=gam2[2,3])
    udata2<-matrix(c(u1,u2,u3),ncol=3) # TODO make sure it's wrapped correctly
    RVM <- VineCopula::RVineMatrix(M2,F2,cop2)
    seq_RVM <-VineCopula::RVineSeqEst(udata2, RVM, weights=r[,2],  progress=FALSE)
    cop2<-seq_RVM$par

  }
  return(list(lls[1:t],gam1,gam2,cop1,cop2,mix))
}


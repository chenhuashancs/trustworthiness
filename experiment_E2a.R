dat <- read.csv("fpfn_realdat.csv")
dat <- signif(dat,3)
nseq <- c(seq(5,35,by=10),47)
nsample <- 1000

m <- 100000
pi1 <- 0.58579
m1 <- m*pi1
m0 <- m-m1
datdir <- "C:\\where_outputs_from_experiment_E2.R_is_stored\\"

# List of bias, RB, and estimates of naive, adjusted and Kant methods
nbiaslist <- biaslist <- kbiaslist <- nrbiaslist <- rbiaslist <- krbiaslist <- nestlist <- estlist <- kestlist <- as.list(NULL)
for (n in nseq) {
  nbiasmat <- biasmat <- ncovmat <- covmat <- as.list(NULL)
  datfile <- paste(datdir,"outapp_n",n,"_new.R",sep="")
  datfile1 <- paste(datdir,"outapp_n",n,"_new_Kant.R",sep="")
  ## True FPRs and FNRs
  pp <- dat$pp[1:n] 
  pn <- dat$pn[1:n]
  ## True positive and negative predictive values
  qp <- pi1*(1-pn)/(pi1*(1-pn)+(1-pi1)*pp)
  qn <- (1-pi1)*(1-pp)/((1-pi1)*(1-pp)+pi1*pn)
  source(datfile)
  # bias and estimate from naive and adjusted methods
  pi1bias <- apply(epi1mat-pi1,2,mean)
  ppbias <- apply(t(eppmat)-c(pp,pp),1,mean)
  pnbias <- apply(t(epnmat)-c(pn,pn),1,mean)
  qpbias <- apply(t(eqpmat)-c(qp,qp),1,mean)
  qnbias <- apply(t(eqnmat)-c(qn,qn),1,mean)
  nbiaslist[[n]] <- c(pi1bias[1],ppbias[1:n],pnbias[1:n],qpbias[1:n],qnbias[1:n])
  biaslist[[n]] <- c(pi1bias[2],ppbias[(1:n)+n],pnbias[(1:n)+n],qpbias[(1:n)+n],qnbias[(1:n)+n])
  nrbiaslist[[n]] <- signif(abs(nbiaslist[[n]]/c(pi1,pp,pn,qp,qn))*100,5)
  rbiaslist[[n]] <- signif(abs(biaslist[[n]]/c(pi1,pp,pn,qp,qn))*100,5)
  nestlist[[n]] <- c(pi1bias[1]+pi1,ppbias[1:n]+pp,pnbias[1:n]+pn,qpbias[1:n]+qp,qnbias[1:n]+qn)
  estlist[[n]] <- c(pi1bias[2]+pi1,ppbias[(1:n)+n]+pp,pnbias[(1:n)+n]+pn,qpbias[(1:n)+n]+qp,qnbias[(1:n)+n]+qn)
  # bias and estimate from Kant method
  source(datfile1)
  pi1bias <- mean(epi1vec-pi1)
  ppbias <- apply(t(eppmat)-pp,1,mean)
  pnbias <- apply(t(epnmat)-pn,1,mean)
  kbiaslist[[n]] <- round(c(pi1bias,ppbias,pnbias),5)
  krbiaslist[[n]] <- signif(abs(kbiaslist[[n]]/c(pi1,pp,pn))*100,5)
  kestlist[[n]] <- round(kbiaslist[[n]]+c(pi1,pp,pn),5)
}
save(dat,nseq,nbiaslist,biaslist,kbiaslist,nrbiaslist,rbiaslist,krbiaslist,nestlist,estlist,kestlist,file="out_E2.R")




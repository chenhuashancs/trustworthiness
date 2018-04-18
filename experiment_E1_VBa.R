options(scipen=999)
### For comparison with Kantchelian use these values for plowseq and nseq
plowseq <- c(0.1,0.4)
nseq <- c(5,15,35)
### For repeating Experiment III.E. 1), use these values for plowseq and nseq and comment out 
### the parts below that are relevant to processing of the outputs from Kantchelian method.
#plowseq <- seq(0,0.7,by=0.1)
#nseq <- c(5,10,15,25,35)
m <- 50000
pi1 <- 0.2
m1 <- m*pi1
m0 <- m-m1
datdir <- "C:\\where_outputs_from_experiment_E1.R_is_stored\\"

# lists of naive estimate bias, adjusted estimator bias, Kant estimator bias, 
#    naive estimate, adjusted estimate, Kant estimate
nbiaslist <- biaslist <- kbiaslist <- nestlist <- estlist <- kestlist <- as.list(NULL)
for (n in nseq) {
  nbiasmat <- biasmat <- kbiasmat <- nestmat <- estmat <- kestmat <- as.list(NULL)
  for (plow in plowseq)
  {
    datfile <- paste(datdir,"out_plow",10*plow,"n",n,".R",sep="")
    datfile1 <- paste(datdir,"outKant_plow",10*plow,"n",n,".R",sep="")
    # bias and estimate from naive and adjusted methods
    source(datfile)
    pi1bias <- apply(epi1mat-pi1,2,mean)
    ppbias <- apply(t(eppmat)-c(pp,pp),1,mean)
    pnbias <- apply(t(epnmat)-c(pn,pn),1,mean)
    qpbias <- apply(t(eqpmat)-c(qp,qp),1,mean)
    qnbias <- apply(t(eqnmat)-c(qn,qn),1,mean)
    nbiasmat[[10*plow+1]] <- round(c(pi1bias[1],ppbias[1:n],pnbias[1:n],qpbias[1:n],qnbias[1:n]),5)
    biasmat[[10*plow+1]] <- round(c(pi1bias[2],ppbias[(1:n)+n],pnbias[(1:n)+n],qpbias[(1:n)+n],qnbias[(1:n)+n]),5)
    nestmat [[10*plow+1]] <- round(nbiasmat[[10*plow+1]]+c(pi1,pp,pn,qp,qn),5)
    estmat [[10*plow+1]] <- round(biasmat[[10*plow+1]]+c(pi1,pp,pn,qp,qn),5)
    # bias and estimate from Kant method
    source(datfile1)
    pi1bias <- mean(epi1vec-pi1)
    ppbias <- apply(t(eppmat)-pp,1,mean)
    pnbias <- apply(t(epnmat)-pn,1,mean)
    kbiasmat[[10*plow+1]] <- round(c(pi1bias,ppbias,pnbias),5)
    kestmat[[10*plow+1]] <- round(kbiasmat[[10*plow+1]]+c(pi1,pp,pn),5)
  }
  nbiaslist[[n]] <- nbiasmat
  biaslist[[n]] <- biasmat
  kbiaslist[[n]] <- kbiasmat
  nestlist[[n]] <- nestmat
  estlist[[n]] <- estmat
  kestlist[[n]] <- kestmat
}
save(plowseq,nseq,nbiaslist,biaslist,kbiaslist,nestlist,estlist,kestlist,file="out_E1_VB.R")


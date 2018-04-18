plowseq <- seq(0,0.4,by=0.1)
nseq <- c(5,10,15,25,35)
# sequence of half-range for random pp,pn
dseq <- c(0.01,0.02,0.03,0.04,seq(0.05, 0.2, by=0.05))
m <- 50000
pi1 <- 0.2
m1 <- m*pi1
m0 <- m-m1
datdir <- "C:\\where_outputs_of_experiment_VC.R_is_stored\\"

nbiaslistdelta <- biaslistdelta <- as.list(NULL)
for (delta in dseq) {
# lists of naive estimate bias, adjusted estimator bias
nbiaslist <- biaslist <- as.list(NULL)
for (n in nseq) {
  nbiasmat <- biasmat <- as.list(NULL)
  for (plow in plowseq)
  {
    datfile <- paste(datdir,"outapp_plow", plow, "_n", n, "_delta",100*delta,".R",sep="")
    # bias from naive and adjusted methods
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
  }
  nbiaslist[[n]] <- nbiasmat
  biaslist[[n]] <- biasmat
}
nbiaslistdelta[[100*delta]] <- nbiaslist
biaslistdelta[[100*delta]] <- biaslist
}
save(dseq,plowseq,nseq,nbiaslistdelta,biaslistdelta,file="out_VC.R")


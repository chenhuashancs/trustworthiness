plowseq <- seq(0,0.4,by=0.1)
nseq <- c(5,10,15,25,35)
nsample <- 1000

# sequence of half-range for random pp,pn
dseq <- c(0.01,0.02,0.03,0.04,seq(0.05, 0.2, by=0.05))

m <- 50000
pi1 <- 0.2
m1 <- m*pi1
m0 <- m-m1

for (plow in plowseq)
for (n in nseq)
for (delta in dseq)
  {
set.seed(100*(n+delta)+10000*plow)
  
# output data file for each plow
datfile <- paste("outapp_plow", plow, "_n", n, "_delta",100*delta,".R",sep="")

## True FPRs and FNRs
pp <- runif(n,plow,plow+0.1) 
pn <- runif(n,plow,plow+0.1)
## True positive and negative predictive values
qp <- pi1*(1-pn)/(pi1*(1-pn)+(1-pi1)*pp)
qn <- (1-pi1)*(1-pp)/((1-pi1)*(1-pp)+pi1*pn)

# matrices of estimates, each row is a sample, 
# 1st column or 1st n columns is naive estimator, 2nd column or 2nd n columns is adjusted estimator.
epi1mat <- eppmat <- epnmat <- eqpmat <- eqnmat <- NULL

for (isample in 1:nsample) {
cat("plow=", plow, "n=",n,"delta=", delta, "isample=",isample,"...\n")

## Ground truth labels
aa <- rep(0,m)
id1 <- sample(1:m,m*pi1)
id0 <- (1:m)[-id1]
aa[id1] <- 1
## data matrix
xmat <- matrix(aa,m,n)
for (j in 1:n) {
  # pp[j] or pn[j] are from random uniform distribution on [0,2p] if p<delta, [2p-1,1] if p>1-delta, or [p-delta,p+delta] otherwise
  xmat[id0,j] <- rbinom(m0,1,ifelse(pp[j]<delta,runif(m0,0,2*pp[j]),
                                    ifelse(pp[j]>1-delta,runif(m0,2*pp[j]-1,1), runif(m0,pp[j]-delta,pp[j]+delta))))
  xmat[id1,j] <- 1-rbinom(m1,1,ifelse(pn[j]<delta,runif(m1,0,2*pn[j]),
                                      ifelse(pn[j]>1-delta,runif(m1,2*pn[j]-1,1), runif(m1,pn[j]-delta,pn[j]+delta))))
}

yy <- ifelse(apply(xmat,1,sum)>=n/2,1,0)

## Monte carlo computation of a probability involving sum
## of multiple independent Bernoulli random variables.
## Inputs: 
##   pvec - vector of Bernoulli probabilities
##   cval - value in probability
##   ge - TRUE or FALSE, P(>=cval) or P(<cval)
mcp <- function(pvec,cval,ge=TRUE) {
  N <- 10000
  npvec <- length(pvec)
  datwk <- NULL
  for (jj in 1:npvec) datwk <- cbind(datwk,rbinom(N,1,pvec[jj]))
  if(ge) val <- mean(apply(datwk,1,sum)>=cval)
  else val <- mean(apply(datwk,1,sum)<cval)
  val
}

## True probabilities
p11 <- mcp(1-pn,n/2)
p01 <- mcp(pp,n/2)
p101 <- p100 <- p110 <- p111 <- rep(0,n)
p001 <- p000 <- p010 <- p011 <- rep(0,n)
for (j in 1:n) {
  p101[j] <- (1-pn[j])*mcp((1-pn)[-j],(n-2)/2,FALSE)
  p100[j] <- pn[j]*mcp((1-pn)[-j],n/2,FALSE)
  p110[j] <- pn[j]*mcp((1-pn)[-j],n/2)
  p111[j] <- (1-pn[j])*mcp((1-pn)[-j],(n-2)/2)
  p001[j] <- pp[j]*mcp(pp[-j],(n-2)/2,FALSE)
  p000[j] <- (1-pp[j])*mcp(pp[-j],n/2,FALSE)
  p010[j] <- (1-pp[j])*mcp(pp[-j],n/2)
  p011[j] <- pp[j]*mcp(pp[-j],(n-2)/2)
}

## Naive estimators
nepi1 <- sum(yy)/m
nepp <- apply(xmat & (1-yy),2,sum)/sum(1-yy)
nepn <- apply((1-xmat)&yy,2,sum)/sum(yy)
neqp <- apply(xmat&yy,2,sum)/apply(xmat,2,sum)
neqn <- apply((1-xmat)&(1-yy),2,sum)/apply(1-xmat,2,sum)

# Estimated probabilities using naive estimators
ep11 <- mcp(1-nepn,n/2)
ep01 <- mcp(nepp,n/2)
ep101 <- ep100 <- ep110 <- ep111 <- rep(0,n)
ep001 <- ep000 <- ep010 <- ep011 <- rep(0,n)
alpha1 <- beta1 <- alpha2 <- beta2 <- rep(0,n)
for (j in 1:n) {
  alpha1[j] <- mcp((1-nepn)[-j],(n-2)/2,FALSE)
  alpha2[j] <- mcp((1-nepn)[-j],n/2,FALSE)
  beta1[j] <- mcp(nepp[-j],(n-2)/2,FALSE)
  beta2[j] <- mcp(nepp[-j],n/2,FALSE)
  ep101[j] <- (1-nepn[j])*alpha1[j]
  ep100[j] <- nepn[j]*alpha2[j]
  ep110[j] <- nepn[j]*(1-alpha2[j])
  ep111[j] <- (1-nepn[j])*(1-alpha1[j])
  ep001[j] <- nepp[j]*beta1[j]
  ep000[j] <- (1-nepp[j])*beta2[j]
  ep010[j] <- (1-nepp[j])*(1-beta2[j])
  ep011[j] <- nepp[j]*(1-beta1[j])
}


# Adjust estimate for pi1
epi1 <- max((nepi1-ep01)/(ep11-ep01),0)
aem1 <- round(m*epi1)
aem0 <- m - aem1
# Adjust estimates for p+, p-, q+, q-
epp <- epn <- eqp <- eqn <- rep(0,n)
for (j in 1:n) {
  a1wk <- aem0*beta1[j]*(1-nepp[j])+aem0*beta2[j]*nepp[j]
  b1wk <- aem1*alpha1[j]*(1-nepp[j])+aem1*alpha2[j]*nepp[j]
  c1wk <- aem0*beta2[j]*nepp[j]-aem1*alpha1[j]*(1-nepp[j])
  a2wk <- aem0*(1-beta2[j])*(1-nepn[j])+aem0*(1-beta1[j])*nepn[j]
  b2wk <- aem1*(1-alpha2[j])*(1-nepn[j])+aem1*(1-alpha1[j])*nepn[j]
  c2wk <- aem0*(1-beta2[j])*(1-nepn[j])-aem1*(1-alpha1[j])*nepn[j]
  wk <- solve(matrix(c(a1wk,a2wk,-b1wk,-b2wk),2),c(c1wk,c2wk))
  epp[j] <- wk[1]
  epn[j] <- wk[2]
  eqp[j] <- epi1*(1-epn[j])/(epi1*(1-epn[j])+(1-epi1)*epp[j])
  eqn[j] <- (1-epi1)*(1-epp[j])/((1-epi1)*(1-epp[j])+epi1*epn[j])
}

epi1mat <- rbind(epi1mat,c(nepi1,epi1))
eppmat <- rbind(eppmat,c(nepp,epp))
epnmat <- rbind(epnmat,c(nepn,epn))
eqpmat <- rbind(eqpmat,c(neqp,eqp))
eqnmat <- rbind(eqnmat,c(neqn,eqn))
dump(c("pi1","pp","pn","qp","qn","epi1mat","eppmat","epnmat","eqpmat","eqnmat"),datfile)


  } # end for (isample in 1:nsample)
} # end for (plow in plowseq)








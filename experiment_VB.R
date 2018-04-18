plowseq <- c(0.1,0.4)
nseq <- c(5,15,35)
nsample <- 1000

# stopping criterion constants for the Kantchelian paper
delta <- 1e-12
epsilon <- 1e-3

for (plow in plowseq)
  for (n in nseq)
  {
# random seed for p+,p-,q+,q- generation    
set.seed(10*(plow+n))
# output data file for each plow
datfile <- paste("outKant_plow",10*plow,"n",n,".R",sep="")

# hyperparameters in theh Kantchelian paper
theta <- 0
phi <- 1.7*n
psi <- 0.1*n

# m is N and n is M in the Kantchelian paper
m <- 50000
#n <- 20
## True malicious proportion
pi1 <- 0.2
m1 <- m*pi1
m0 <- m-m1
## True FPRs and FNRs
pp <- runif(n,plow,plow+0.1) 
pn <- runif(n,plow,plow+0.1)
## True positive and negative predictive values
qp <- pi1*(1-pn)/(pi1*(1-pn)+(1-pi1)*pp)
qn <- (1-pi1)*(1-pp)/((1-pi1)*(1-pp)+pi1*pn)

# vector/matrices of estimates, each entry/row is a sample, pp=beta, pn=1-alpha
epi1vec <- eppmat <- epnmat <- NULL

for (isample in 1:nsample) {
# random seed for file labels generation    
set.seed(10000*(plow+n)+isample-1)
cat("plow=",plow,"n=",n,"isample=",isample,"...\n")
## Ground truth labels
aa <- rep(0,m)
id1 <- sample(1:m,m*pi1)
id0 <- (1:m)[-id1]
aa[id1] <- 1
## data matrix
xmat <- matrix(aa,m,n)
for (j in 1:n) {
  xmat[id0,j] <- rbinom(m0,1,pp[j])
  xmat[id1,j] <- 1-rbinom(m1,1,pn[j])
}

# initialize z vector for the E-M algorithm as in 5.1 of the Kantchelian paper
zmean <- ifelse(apply(xmat,1,sum)>4,1,0)
epi1 <- 0.5; ealpha <- ebeta <- rep(0.5,n)
itn <- 1; cntdiv <- 0
# loop for the E-M algorithm
repeat{
  if(itn%%100==0) cat("EM iteration", itn, "...\n")
  epi1new <- (theta+sum(zmean))/(2*theta+m)
  ealphanew <- apply(zmean*xmat,2,sum)/(phi+sum(zmean))
  ebetanew <- apply((1-zmean)*xmat,2,sum)/(psi+sum(1-zmean))
  disc <- sum(abs(epi1new-epi1),abs(ealphanew-ealpha),abs(ebetanew-ebeta))
  # clip the new parameter values to the interval [epsilon,1-epsilon] if necessary
  epi1 <- ifelse(epi1new<=epsilon, epsilon, epi1new); epi1 <- ifelse(epi1>=1-epsilon,1-epsilon,epi1)
  ealpha <- ifelse(ealphanew<=epsilon, epsilon, ealphanew); ealpha <- ifelse(ealpha>=1-epsilon,1-epsilon,ealpha) 
  ebeta <- ifelse(ebetanew<=epsilon, epsilon, ebetanew); ebeta <- ifelse(ebeta>=1-epsilon,1-epsilon,ebeta) 
  amatwk <- matrix(rep(ebeta/ealpha,rep(m,n)),nrow=m)
  bmatwk <- matrix(rep((1-ebeta)/(1-ealpha),rep(m,n)),nrow=m)
  zmean <- 1/(1+(1-epi1)/epi1*apply(ifelse(xmat,amatwk,bmatwk),1,prod))
  if (disc<=delta) break
  itn <- itn + 1
  if (itn > 500) {
    cntdiv <- cntdiv + 1
    break
  }
}

epi1vec <- c(epi1vec,epi1)
eppmat <- rbind(eppmat,ebeta)
epnmat <- rbind(epnmat,1-ealpha)
dump(c("pi1","pp","pn","qp","qn","epi1vec","eppmat","epnmat"),datfile)

  } # end for (isample in 1:nsample)
} # end for (plow in plowseq)



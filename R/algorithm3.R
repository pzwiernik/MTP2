#library(MASS)
### simulate data
#n <- 100
#p <- 10
#Sigma <- diag(p)
#X <- mvrnorm(n, rep(0,p), Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#S <- (t(X)%*%X)/n

##
library(gRbase)
data("carcass")
dat <- carcass
S <- cov(dat)[1:6,1:6]
p <- 6

##### Algorithm 3
K <- solve(diag(diag(S)))
# start with G, the complete graph
K0 <- 100*diag(p)
it <- 0

while (sum(abs(K-K0))>1e-16){
K0 <- K
it <- it+1
for (i in 1:(p-1)){
  for (j in (i+1):p){
    A <- union(i,j)
    B <- setdiff(1:6,A)
    L <- K[A,B]%*%solve(K[B,B])%*%K[B,A]
    if (L[1,2]<=S[i,j]/(S[i,i]*S[j,j]-S[i,j]^2)) K[A,A] <- solve(S[A,A])+L
    else {
      au <- (1+sqrt(1+4*L[1,2]^2*S[i,i]*S[j,j]))/2
      K[A,A] <- diag(c(L[1,1]+au/S[i,i],L[2,2]+au/S[j,j]))
     }
  }
}
}
K
Sig <- solve(K)
cov2cor(Sig)

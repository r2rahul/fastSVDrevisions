library(Matrix)

## Build the intial Matrix and factorization
m <- 17 # m > n, this condition must hold 
n <- 10
X <- Matrix(runif(m*n,min=0,max=5),ncol =n,nrow=m, sparse = TRUE)
test <- svd(X,nu=n)

### Let us append a row to X equivalent to real time update.
new.row <- Matrix(runif(n,min=0,max=2),ncol=1, sparse = TRUE)
w <- m+1
A <- Matrix(rep(0,w),ncol=1, sparse = TRUE)
B <- new.row
A[w,1] <- 1 
U <- Matrix(test$u, sparse = TRUE)
V <- Matrix(test$v, sparse = TRUE)
S <- Matrix(diag(test$d), sparse = TRUE)
temp <- Matrix(rep(0,n),nrow=1,sparse = TRUE)
U0 <- rBind(U,temp)
X0 <- rBind(X, temp)

## Perform the Rank-1 Update
update.svd <- increment_svd(U0, S, V, A, B)

## Final check for the correction of Algorithm
U2 <- update.svd$U
V2 <- update.svd$V
S2 <- update.svd$S
temp <- U2 %*% S2
error <- base::norm((X0 + A %*% t(B)) -  temp %*% t(V2),type="2")
print(error)


increment_svd <- function(U, S, V, A, B){
  d_u <- dim(U)
  d_a <- dim(A)
  rankU <- d_u[2]
  qrP <- base::qr(cBind(U,A))
  QP <- Matrix(qr.Q(qrP, complete = TRUE), sparse = TRUE)
  RP <- Matrix(qr.R(qrP, complete = TRUE), sparse = TRUE)
  w <- rankU+1
  d_qp <- dim(QP)
  d_rp <- dim(RP)
  
  
  if (d_qp[2] == rankU) {
    P <- Matrix(numeric(0),ncol=0,nrow=0,sparse = TRUE)
    Ra <- Matrix(numeric(0),ncol=0,nrow=0,sparse = TRUE)
    dim.Ra <- 0; 
  }
  else {
    P <- QP[ ,w:d_qp[2]]
    Ra <- RP[w:d_rp[1], w:d_rp[2]]
    dim.Ra <- d_rp[1] - rankU
  }
    
  
  v <- rankU+d_a[2]
  M <- Matrix(RP[1:rankU, w:v], sparse=TRUE)
  
  
  qrQ <- base::qr(cBind(V,B))
  QQ <- Matrix(qr.Q(qrQ, complete = TRUE), sparse=TRUE)
  RQ <- Matrix(qr.R(qrQ, complete = TRUE), sparse=TRUE)
  d_qq <- dim(QQ)
  d_rq<- dim(RQ)
  
  if (d_qq[2] == rankU) {
    Q <- Matrix(numeric(0), ncol=0,nrow=0, sparse = TRUE)
    Rb <- Matrix(numeric(0), ncol=0,nrow=0, sparse = TRUE)
    dim.Rb <- 0
  }
  else {
    Q <- QQ[,w:d_qq[2]] 
    Rb <- RQ[w:d_rq[1], w:d_rq[2]]
    dim.Rb <- d_rq[1] - rankU
  }
    
  N <-  t(V) %*% B 
  ##########################################
  ####### Creating Augmented Sigma
  ##########################################
  #l <- d_u[1]-nrow(M)-1
  #Ra <- rBind(Ra,Matrix(rep(0,l),ncol=1))
  temp0 <- rBind(M,Matrix(Ra,ncol=1, sparse = TRUE))
  if (ncol(Rb) == 0){
    temp1 <- N
  }
  else {
    temp1 <- rBind(N,Rb)
  }
  
  S.augmented <- temp0 %*% t(temp1)
  
  if (dim.Rb == 0) {
    
    temp2 <- cBind(S, Matrix(numeric(0),
                             nrow=rankU, ncol=dim.Rb,
                             sparse = TRUE))
    
    temp3 <- cBind(Matrix(rep(0, rankU*dim.Ra),
                          nrow=dim.Ra, ncol=rankU,
                          sparse = TRUE), 
                   Matrix(numeric(0),
                          nrow=dim.Ra, ncol=dim.Rb,
                          sparse = TRUE))
  }
  else {
    i <- rankU* dim.Rb 
    temp2 <- cBind(S, Matrix(rep(0,i),
                             nrow=rankU, ncol=dim.Rb,
                             sparse = TRUE))
    
    j <- dim.Ra*dim.Rb
    temp3 <- cBind(Matrix(rep(0, rankU*dim.Ra),
                          nrow=dim.Ra, ncol=rankU,
                          sparse = TRUE), 
                   Matrix(rep(0, j),
                          nrow=dim.Ra, ncol=dim.Rb,
                          sparse = TRUE))
  }
  
  
 
  S.new <- rBind(temp2, temp3) + S.augmented
  
  update.svd <- svd(S.new, nu =rankU)
  
  U.new <- cBind(U,P) %*% update.svd$u
  S.new <- update.svd$d
  
  if (nrow(Q) ==0){
    V.new <- V %*% update.svd$v
  }
  else {
    V.new <- cBind(V,Q) %*% update.svd$v
  }
  
  
  U.new <- Matrix(U.new, sparse = TRUE)
  S.new <- Matrix(diag(S.new), sparse = TRUE)
  V.new <- Matrix(V.new, sparse = TRUE)
  
  return(list(U=U.new, S= S.new, V = V.new))
}
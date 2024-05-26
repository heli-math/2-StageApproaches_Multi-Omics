# function 'generate_data_inclZ'
# based on the function 'generate_data' in package PO2PLS
# aim: generate X, Y, Z at same time in normal distr

generate_data_inclZ <- function(N, params, distr = rnorm){
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  SigT = params$SigT
  SigTo = params$SigTo + 1e-6*SigT[1]*(params$SigTo[1]==0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT%*%B^2 + SigH
  SigUo = params$SigUo + 1e-6*SigU[1]*(params$SigUo[1]==0)
  a = params$a # added, different from PO2PLS
  b = params$b # added
  sig2G = params$sig2G # added
  
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  
  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co),
                cbind(a, matrix(0,1,r+rx+ry))) # added
  VarM = blockm( # change from VarZ
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)

  M <- scale(matrix(distr(N*(2*r+rx+ry)), N)) # change Z into M
  M <- M %*% chol(VarM)
  M[,2*r+1:rx] <- sign(ssq(Wo))*M[,2*r+1:rx]
  M[,2*r+rx+1:ry] <- sign(ssq(Co))*M[,2*r+rx+1:ry]
  
  EFG <- cbind(scale(matrix(distr(N*p), N))*sqrt(sig2E), 
               scale(matrix(distr(N*q), N))*sqrt(sig2F),
               scale(matrix(rnorm(N), N))*sqrt(sig2G)) # added
  
  dat <- M %*% t(Gamma) + EFG # added
  return(list(X = dat[,1:p], Y = dat[,(p+1):(p+q)], Z = dat[,-(1:(p+q))]))
}

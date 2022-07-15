hadamardF <- function(n){
  if(n <= 1){
    res <- matrix(0:1, ncol=2)
    return(res)
  } else {
    top <- matrix(c(numeric(2^(n-1)), rep(1,(2^(n-1)))),
                            nrow=1, ncol=(2^n), byrow=TRUE)
    bot <- cbind(hadamardF((n-1)), hadamardF((n-1)))
    res <- rbind(top, bot)
    return(res)
  }
}
hadamard <- function(n) {
  F  <- hadamardF(floor(log2(n)))
  FF <- t(F)%*%F
  FF[FF%%2==1] <- -1
  FF[FF%%2==0] <- 1
  return(FF)
}

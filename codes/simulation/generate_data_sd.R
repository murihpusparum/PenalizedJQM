# function to obtain seed that can generate datasets without warning in LQMM

source("./codes/main/LQMM_Function.R")

d<-function(alpha, N, n, beta, a, b, seed, NSim, FUN){
  seed.d<-matrix(nrow=NSim,ncol=10)
  seed.d<-as.data.frame(seed.d)
  names(seed.d)<-c("iter", "N", "n", "alpha", "a", "b","new","min.sd","max.sd","sd.sd")
  set.seed(seed)
  for (si in 1:NSim) {
    seed.n<-runif(1, 2904713, 16397462)
    res<-lqmm.fun(alpha=alpha, N=N, n=n, beta=beta, a=a, b=b, seed=seed.n, si=si, FUN=FUN) 
    df<-res$df
    seed.d[si,1:10]<-c(si,N,n,alpha,a,b,res$seed,min(res$sd),max(res$sd),sd(res$sd))
  }
  return(list(seed=seed.d,
              df=df))
}
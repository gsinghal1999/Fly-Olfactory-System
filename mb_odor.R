odor <-function(no,params=c(.1,100)){
  size <- no
  signal <- rexp(size,params)
  signal <- signal*(1+rnorm(no,0,0.1)) #noise
  signal
}

conn_matrix <- function(r,c){#at six random locations synapses are present randomly
    matrix <- as.vector(do.call("cbind", lapply(1:r, function(x) sample(c(rgamma(6,4,4), rep(0, 44))))))
    matrix <- matrix*(1+sapply(1:c,function(x) rnorm(r,mean = 0,sd = 0.1))) #noise
}

mb_odor <-function(signal, matrix,op=1){
  mb <- signal %*% matrix 
  if (is.null(signal)) return(as.vector(mb)) 
  mb
}

mb_new<- function(mb){ #negative feedback
  new_odor<- 0.1*mb
  new_odor_len<-length(new_odor) #length of final odor signal
  new_odor <- as.vector(new_odor*(1+rnorm(new_odor_len,0,0.1))) #noise added
}

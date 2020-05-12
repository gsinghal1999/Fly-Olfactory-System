#this function contains code to implement the model from the OB to the PCx



#ths function geenerates the OB code

#for binary the glomeruli are randomly chosen unless specified by by bulb vec

#no - # glomeruli

#type, 1 - binary on/off, 2 - exponential distribution, 3- nearly binary - the others are

# 4 - gamma distribution

#not 0 but within 10 % of 0, no negative firing rates

#params - if type ==1, (fraction active,firing rate), if type = 2 (exp. parameter)

#noise: the noise to be added, it is gaussian mean 0 with stdev given by noise

#op: 1: generate active glomeruli randomly, 2- use the obvec vector

genOBCode <-function(no,type=1,params=c(.1,100),obvec=c(),noise=0,op=1){

  #cat('genOBcode',no,type)

  size <- no

  if(type==1){#binary code

    tmpn <- params[1]*no #no of active glom

    #cat('genOBcode binary',params[1],',',no,',',tmpn)

    if (length(obvec)==0) glom <- sample(no)[1:tmpn] #the glomeruli that are active

    else glom <- obvec

    glomvec <- sapply(1:no,function(x) 0) #rep(0,times=size) # initialize glomvec

    #cat('first glm',length(glomvec),'no',no)

    if(length(glomvec)<no) cat('glomvec',glomvec,'\t')

    glomvec[glom] <- params[2] #set the active glom to the firing rate

  }

  if(type==2){#exponential code, simpler, just generate the code

    glomvec <- rexp(size,params)

  }

  if(type==4){#gamma distribution

    glomvec <- rgamma(size,shape = params[1],scale = params[2])

  }

  #add noise, noise % is added to the vectors firing rate 

  #tst <- sapply(1:no,function(x) rnorm(10,0,noise))

  #if (length(tst) < no ) cat(' lesser ')

  #cat('obcode',no,'glomvec noise length',length(tst))

  glomvec <- glomvec*(1+rnorm(no,0,noise))

  #cat('2nd obcode',no,'glomvec',length(glomvec))

  glomvec

}



#ths function geenerates a vector of OB code, returns a matrix where each odor is a col.

#for binary the glomeruli are randomly chosen unless specified by by bulb vec

#no - # glomeruli

#type, 1 - binary on/off, 2 - exponential distribution, 3- nearly binary - the others are 

#not 0 but within 10 % of 0, no negative firing rates

#params - if type ==1, (fraction active,firing rate), if type = 2 (exp. parameter)

#noise: the noise to be added, it is gaussian mean 0 with stdev given by noise

#n: no of odors

#op: 1: generate active glomeruli randomly, 2- use the obvec vector

genOBCodeVec <-function(n,no,type=1,params=c(.1,100),obvec=c(),noise=0,op=1){

   res <- lapply(1:n,function(x) genOBCode(no,type = type,params = params,obvec = obvec,

                                           noise = noise,op = op))

   combineListMat(res)

}



#given a bindary odor code adds some noise to the zero and nonzero posns

#dist: distribution to be used

#params : tje sd of this distribution

genOBCodeBinary <- function(ob,dist=1,params=.2,op=1){

  maxr <- max(ob)

  minr <- min(ob)

  #they both are around 10% of their firing rates

  zeropos <- which(ob == minr)

  nonzeropos <- which(ob != minr)

  res <- ob

  res[zeropos] <- rnorm(length(zeropos),mean = minr,sd = maxr*params)

  res[nonzeropos] <- rnorm(length(nonzeropos),mean = maxr,sd = maxr*params)

  scaleMat(res) #get rid of non-zero entries

}



#this function generates a second ob vector that overlaps with the first one to n percent

#caution: curently if the number of overlap posns comes to be less than 1 returns 0

#type, 1 - binary on/off, 2 - exponential distribution, 3- nearly binary - the others are 

#not 0 but within 10 % of 0, no negative firing rates

# 4 - a gamma distribution

#noise: the noise to be added, it is gaussian mean 0 with stdev given by noise

#overlap: percentage of overlap

genOverlapOBCode <-function(obvec,overlap=50,type=1,params=100,noise=0,op=1){

  glomlen <- length(obvec) #length of the ob vector

  glompos <- 1:glomlen #vec of glomeruluar posns

  if(type==1){

    actpos <- which(obvec>0) #get all the active positions

    posns <- sample(actpos,length(actpos)*overlap/100) #sample the overlapping fraction

    inactpos <- setdiff(glompos,actpos) #get all non-active posns

    otherposns <- sample(inactpos,length(actpos)*(1-overlap/100)) # sample the rest from here

    glomvec <- rep(0,glomlen)

    #cat(' posns',posns,',',otherposns)

    glomvec[c(posns,otherposns)] <- params

  }

  if(type==2){

    noposns <- (1-(overlap/100)) * glomlen #noposns to replace, i.e. overlap fraction

    replaceposns <- sample(glompos,noposns) #choose the posns

    #cat(glomlen,',',noposns,'m',replaceposns)  

    glomvec <- obvec

    glomvec[replaceposns] <- rexp(noposns,params)

  }

  if(type==4){#gamma distribution

    noposns <- (1-(overlap/100)) * glomlen #noposns to replace, i.e. overlap fraction

    replaceposns <- sample(glompos,noposns) #choose the posns

    #cat(glomlen,',',noposns,'m',replaceposns)  

    glomvec <- obvec

    glomvec[replaceposns] <- rgamma(noposns,shape = params[1],scale = params[2])

  }

  #add noise, noise % is added to the vectors firing rate 

  glomvec <- glomvec*(1+rnorm(glomlen,0,noise))

  glomvec

}





#adds noise to the glomerular vector

#noise to be added to the glomerular vector

addNoiseGlom <-function(glomvec,noise=0.1){

  #get the dimension of the vector

  glomlen <- length(glomvec)

  #add noise, noise % is added to the vectors firing rate 

  glomvec <- glomvec*(1+rnorm(glomlen,0,noise))

  glomvec

  

}



#generates connection matrix, default gaussiam, whose params are specified by params

#the matrix is size rows x cols

#rows - the rows, for eg each one might specify a pcx neuron, the target neurons

#cols - the cols e.g each one might specify a glomerulus, the source neurons

#type - 1- agaussian matrix, 2 -gamma distribution, 3 - compound poisson Gamma distribution

#nosyn - not needed: average number of synapses between a glomerulus and neuron 

#noise: the noise to be added, it is gaussian mean 0 with stdev given by noise

#noise, e.g. noise = .1 is 10 % noise

#if compound poisson gamma params = (poisson mean,alpha,beta)

#params - , no of synapses, and paramss of the distribution

genConnMat <- function(rows,cols,params=c(1,1.15,.16),type=3,noise=0,op=1){

  if(type==1){#gaussian matrix

    mat <- matrix(rnorm(rows*cols,mean = params[1],sd = params[2]),nrow = rows)

  }

  if(type==3){#compund poisson gamma matrix, could be made faster by breaking all 

    #generate 1.2 times total synapse no of gamma distributions, and then use them

    #although since the no of synapses is from the mean, we shouldnt need them

    strvec <- rgamma(rows*cols*params[1]*1.2,shape=params[2],scale=params[3])

    synno <- rpois(rows*cols,params[1]) #no of synapses for all connections

    index <- cumsum(synno) #the cumulative sum which acts as an index

    zero <- synno #multiplier to turn 0 synapses to 0

    zero[zero > 0] <- 1

    #cat(index,'\n',synno,'\nzero',zero,'\n')

    mat <- sapply(1:(rows*cols),function(i) {

      #strength <- sapply(0:syn,function(x) (x>0)*rgamma(1,shape=params[2],scale=params[3]))

      #start and end give the indices of the strvec that need to be used

      start <- (index[i]-synno[i]+1)*zero[i] 

      end <- index[i]*zero[i]

      #cat(i,'st ',start,',',' end',end,'\t',sum(start,end),' ;')

      strength <- sum(strvec[start:end])

      strength

    })

    mat <- matrix(mat,nrow = rows)

  }

  if(type==2){#gamma matrix

    mat <- matrix(rnorm(rows*cols,mean = params[1],sd = params[2]),nrow = rows)

  }

  #adding noise if needed

  mat*(1+sapply(1:cols,function(x) rnorm(rows,mean = 0,sd = noise)))

}



#add noise to the synaptic conn matrix

#noise, e.g. noise = .1 is 10 % noise

addConnMatNoise <-function(mat,noise=0.1){

  rows <- dim(mat)[1]

  cols <- dim(mat)[2]

  #cat('rows',rows,cols,str(mat))

  mat*(1+sapply(1:cols,function(x) rnorm(rows,mean = 0,sd = noise)))

}





#given the conn matrix and OB code, this generates the PCx vector

#mat: matrix of pcx rows vs ob cols

genPcxVec <-function(obvec,mat,op=1){

  #pcx <- obvec %*% t(mat) #generate the pcx vector or matrix

  pcx <- mat %*% obvec #generate the pcx vector

  if (is.null(obvec)) return(as.vector(pcx)) #return a vector if obvec is a vector

  pcx

}







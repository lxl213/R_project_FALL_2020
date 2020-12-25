##########################
# Dec.2020 by Yan Ruoru and Luo Xiaolong adopted by author's code
# Convex sliced inverse regression code 
# For three part:
# Part1:Estimating the conditional covariance matrix 
# Part2:Convex formulation for SSIR
# Part3:Cross validation to choosing the best rank K and tuning parameter 
###########################

# Part1:
# Estimating the conditional covariance
# Using cov(E[X|Y])=cov(x)-E[cov(X|Y)] by estimating E[cov(X|Y)]
# according to equation (7) by computing weighted average of the sample covariance

calsigmafit <- function(y,X,nslice){
  #input: y:a vector of response; X:a n*d matrix of observations
  #       nslice:the number of slices; 
  #output:a matrix of estimator of X
  
  n <- length(y)
  X <- scale(X,T,F) #standardize

  
  # Test if this is integer or continuous
  if(sum(y-round(y))==0){
    nslice <- length(unique(y))
    quany <- sort(unique(y))
  }
  
  # For continuous response y
  else if(sum(y-round(y))!=0){
    quany <- quantile(y,seq(1/nslice,1,length.out=nslice))
  }
  
  indexy <- vector("list",nslice) #index vector initialization
  
  # Assign indices into indexy
  indexy[[1]] <- which(y<=quany[1])
  for(k in 2:nslice){
    indexy[[k]] <- which(y >quany[k-1] & y<=quany[k])
  }
  
  #the number of observations in each slice
  nindexy <- lapply(indexy,length)
  
  f <- matrix(0,n,nslice)
  for(k1 in 1:(nslice-1)){
    for(k2 in 1:nslice){
      if(k1==k2){
        f[indexy[[k1]],k2] <- 1 - nindexy[[k2]]/n
      }
      if(k1!=k2){
        f[indexy[[k1]],k2] <- -nindexy[[k2]]/n			
      }
    }
  }
  for(k in 1:nslice){
    f[indexy[[nslice]],k] <- -nindexy[[k]]/n
  }
  
  bigF <- f%*%solve(t(f)%*%f)%*%t(f)
  Sigmafit <- t(X)%*%bigF%*%X/(n)
  return(Sigmafit)
}

##########################################################
# Part2
# Algorithm 1:LADMM to solve convex optimazation problem 9
########################################################## 

ssir <- function(covxy,covx,lambda,K,nu=1,epsilon=1e-3,maxiter=1000,trace=FALSE,init=FALSE,initPi=NULL,initH=NULL,initGamma=NULL){
  #input:cov_xy:the estimator of conditional covariance in part 1
  #      cov_x:the estimator of covariance of observation x
  #      lambda:the tuning parameter; K:rank constraint; nu:L-ADMM parameter
  #      epsilon:tolerance; maxiter:maxium rounds of iteration trace:to tract the performing
  #output:a list including object variable pi,H and dual variable gamma
  
  p <- nrow(covx) 
  eigencovx <- eigen(covx) 
  #computing the cov_x decomposition of cov_x^(1/2)
  sqcovx <- eigencovx$vectors%*%sqrt(diag(pmax(eigencovx$values,0)))%*%t(eigencovx$vectors)	
 
  tau <- 4*nu*eigencovx$values[1]^2	
  criteria <- 1e10
  i <- 1


# Initialize 
   H <- Pi <- oldPi <-  diag(1,p,p)
   Gamma <- matrix(0,p,p)

   if(init==TRUE){
      H <- initH #H when having input in cross-validation
      Pi <- initPi 
      Gamma <- initGamma 
    }

# While loop for the iterations
  while(criteria > epsilon && i <= maxiter){
    Pi <- newPi(covx,sqcovx,covxy,H,Gamma,nu,lambda,Pi,tau)
    H <- newH(sqcovx,Gamma,nu,Pi,K)
    Gamma <- Gamma + sqcovx%*%Pi%*%sqcovx-H	
    criteria <- sqrt(sum((Pi-oldPi)^2)) 
    oldPi <- Pi
    i <- i+1
  
  #the tract decision
    if(trace==TRUE)
    {
      print(i)
      print(criteria)
    }
  
  }  

  return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
}


# Soft-thresholding Operator
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))
}


# Update Pi
newPi <- function(covx,sqcovx,covxy,H,Gamma,nu,lambda,Pi,tau){
  A <- Pi + 1/tau*covxy-nu/tau*covx%*%Pi%*%covx+nu/tau*sqcovx%*%(H-Gamma)%*%sqcovx
  B <- lambda/tau
  return(Soft(A,B))
}


# Update H
newH <- function(sqcovx,Gamma,nu,Pi,K){
  
  temp <- Gamma + sqcovx%*%Pi%*%sqcovx
  temp <- (temp+t(temp))/2 
  svdtemp <- eigen(temp)
  d <- svdtemp$values
  p <- length(d)
  
  if(sum(pmin(1,pmax(d,0)))<=K){
    dfinal <- pmin(1,pmax(d,0))
    return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
  }
  
  fr <- function(x){ # a simple operator with input of gamma 
    sum(pmin(1,pmax(d-x,0)))
  }
  
  knots <- sort(unique(c((d-1),d)),decreasing = TRUE)
  temp <- which(sapply(knots,fr)<=K)
  
  lentemp <- tail(temp,1) #The length of singular values meeting command
  a=knots[lentemp]
  b=knots[lentemp+1]
  fa <- sum(pmin(pmax(d-a,0),1))
  fb <- sum(pmin(pmax(d-b,0),1))
  theta <- a+ (b-a)*(K-fa)/(fb-fa) #Lagrange inserting method
  dfinal <- pmin(1,pmax(d-theta,0)) 
  res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
  return(res)
}	

################################
# Part3:
# cross-validation to select tuning parameter K and lambda
################################
predictpfc <- function(pfcobject,K,y,X,Xnew){
  #input: pfcobject:the object by training set; K:constraint rank; 
  # y:training response; X:training set Xnew:testing set of observation 
  #output:a vector of prediction y of ytest  
  
  Pi <- pfcobject$Pi 
  Pi <- (t(Pi)+Pi)/2 #transformed into symmetric matrix 
  temp <- eigen(Pi)
  if(max(temp$values)<0.01){
    return(rep(mean(y),nrow(Xnew)))
  }
  temp <- temp$vectors[,1:K]
  RhatX <- X%*%temp
  
  Xnew <- as.list(data.frame(t(Xnew)))
  
  predicty <- function(x){ 		
    temp2 <- x%*%temp	
    residual <- t(t(RhatX)-as.vector(temp2)) #residual in each row.
    weights <- exp(-0.5*apply((residual)^2,1,sum))
    weights <- weights/sum(weights)
    return(sum(weights*y))
  }
  yhat <- unlist(lapply(Xnew,predicty))
  return(yhat)	
} 


ssir.cv <- function(X,y,Ks,lambdas,nfold,nslice){
  #input: X:the observations; y:the response vector; 
  # Ks:the possible parameter of rank constraint 
  # lambdas:the possible parameter of sparsity tuning parameter
  # nslice:the number of slices; nfold:the number of partition sets
  #output:a list with elements of a matrix including
  
  p <- ncol(X) #The dimension of primal obeservation
  fold <- rep(NA,length(y))
  fold <- sample(rep(1:nfold,length(y)/nfold)) #The partition
 
   cv.error <- vector("list",length(Ks)) #output 
  for(K in Ks){
    cv.error[[K]] <- matrix(NA,nrow=nfold,ncol=length(lambdas))	
  }
   
  for(j in 1:nfold){
    print(paste("cross validation for dataset ",j,sep=""))
    tmp <- 1
    ytrain <- y[which(fold!=j)]
    ytest <- y[which(fold==j)]		
    Xtrain <- X[which(fold!=j),]		
    Xtest <- X[which(fold==j),]	
    
    Sigmax <- cov(Xtrain)
    #estimator of conditional covariance 
    Sigmafit <- calsigmafit(ytrain,Xtrain,nslice=nslice) 
    
    for(K in Ks){
      initH = initPi = diag(1,p,p)
      initGamma = diag(0,p,p)
      
      tmp <- 1
      for(lambda in lambdas){
        res <- ssir(Sigmafit,Sigmax,lambda,K,epsilon=5e-04,maxiter=1000,init=TRUE,initPi=initPi,initH=initH,initGamma=initGamma,trace=FALSE)
        initPi = res$Pi
        initH = res$H
        initGamma = res$Gamma
        yhat <- predictpfc(res,K,ytrain,Xtrain,Xtest)
        # Test if this is integer or continuous
        #if(sum(y-round(y))==0){
        #	yhat <- sign(yhat)
        # cv.error[[K]][j,tmp] <- 	1-sum(diag(table(yhat,ytest)))/sum(table(yhat,ytest))
        #}
        #else if(sum(y-round(y))!=0){
        cv.error[[K]][j,tmp] <- sum((ytest-yhat)^2)
        #}
        tmp <- tmp + 1
      }
    }
  }
  return(cv.error)
}

##################################################
#part 4:Numerical studies
##################################################

obseX = function(p,n){
  # Generate AR-1 type covariance for X
  # input: p:the number(dimension) of covariates X; n: the number of samples
  # Output: generate a matrix with n rows and p columns with observation X
  # X has the distribution from N(0,sigma_x)??where sigma_x[i,j] = 0.5^abs(i-j)
  Sigma <- matrix(0.5,p,p)
  tmpmat <- matrix(0,p,p)
  for(i in 1:(p-1)){
    tmpmat[i,i:p] <- c(0:(length(i:p)-1))
  }
  tmpmat = tmpmat+t(tmpmat)
  Sigma <- Sigma^tmpmat	
  
  X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
  return(X)
}


simul <- function(X,y,n,p,lambdas,seed,beta){
  #a function for simulation
  #input:X,y,beta:given by specific model indicating observations,responses and direction
  #     n:the number of samples; p:the number of dimension; lambdas:possible tuning parameter
  #Output: a list including TPR,FPR,chosen parameters K and lambda 
  
  
  set.seed(seed)
  
  # Cv to select lambda #
  Ks <- 1:3
  nfold <- 5
  nslice <- 5
  lambdas <- lambdas
  a <- ssir.cv(X,y,Ks,lambdas,nfold,nslice)
  print("a")
  print(a)
  calmean <- function(mat){ apply(mat,2,mean)}
  calse <- function(mat){ apply(mat,2,sd)/nfold}
  temp1 <- lapply(a,calmean)
  print("temp1")
  print(temp1)
  temp1se <- lapply(a,calse)
  print("temp1se")
  print(temp1se)
  cverror <- NULL
  
  for(k in Ks){
    
    tempcv <- min(temp1[[k]])#+temp1se[[k]][which(temp1[[k]]==min(temp1[[k]]))][1]
    
    tempcv2 <- temp1[[k]][tail(which(temp1[[k]]<=tempcv),1)]
    
    cverror <- c(cverror,tempcv2)
  }
  print("cverror")
  print(cverror)
  	temp2 <- unlist(lapply(temp1,min))
  	print("temp2")
  	print(temp2)
  chosenK <- which(cverror==min(cverror))
  chosenK <- chosenK[1]
  print(chosenK)
  
  templambda <- min(temp1[[chosenK]])#+temp1se[[chosenK]][which(temp1[[chosenK]]==min(temp1[[chosenK]]))][1]
  
  chosenlambda <- lambdas[tail(which(temp1[[chosenK]]<=templambda),1)]
  chosenlambda <- chosenlambda[1]
  print(chosenlambda)
  
  # Fit it using CV selected tuning parameters
  Sigmax <- cov(X)
  Sigmafit <- calsigmafit(y,X,nslice=5)
  res <- ssir(Sigmafit,Sigmax,chosenlambda,chosenK,epsilon=1e-04,maxiter=1000,trace=FALSE)
  temp <- eigen(round((res$Pi+t(res$Pi))/2,2))
  print(temp$vectors%*%t(temp$vectors))
  est <- temp$vectors[,1:chosenK]%*%t(temp$vectors[,1:chosenK])
  print(est)
  cvspec <- sum(abs(diag(est))<1e-5 & beta< 1e-5)/sum(beta==0)
  cvsens <- sum(abs(diag(est))>1e-5 & abs(beta)>1e-5)/sum(beta!=0)
  
  return(list(cvsens=cvsens,cvspec=cvspec,K=chosenK,lambda=chosenlambda,temp=temp))
}
###################
#Setting 1:
#y=(x1+x2+x3)/sqrt(3)+2*noise,where noise has distribution N(0,1)
###################
n = 100
p = 150
X = obseX(n=n,p=p)
beta <- rep(0,p)
beta[1:3] <- 1/sqrt(3)
y <- (X[,1]+X[,2]+X[,3])+2*rnorm(n,0,1)
lambdas=0.15

t1=proc.time()
z = simul(X=X,y=y,n=n,p=p,lambdas=lambdas,beta=beta,seed=1)
print(z)
t2=proc.time()
cat("time consuming:",t2-t1,"s")

if(z$K==1){
  abscor <-abs(cor(z$temp$vectors[,1],beta))
  
}
if(z$K>1){
  abscor <- apply(abs(cor(z$temp$vectors[,1:z$K],beta)),2,max)
}
cat("abscor",abscor)
####################
#Setting 2:
#y = 1 + exp((x1+x2+x3)/3^(1/2))+epsilon
####################

n = 100
p = 150
X = obseX(n=n,p=p)
beta <- rep(0,p)
beta[1:3] <- 1/sqrt(3)
#beta = as.matrix(beta)
y <- X%*%beta+2*rnorm(n,0,1)

t1 = proc.time()
lambdas=seq(from=0.02,to=0.5,by = 0.02)
t2 = proc.time()
z=simul(X=X,y=y,n=n,p=p,lambdas=lambdas,beta=beta,seed=2)
cat("time consuming:",t2-t1,"s")
print(z)

if(z$K==1){
  abscor <-abs(cor(z$temp$vectors[,1],beta))
  
}
if(z$K>1){
  abscor <- apply(abs(cor(z$temp$vectors[,1:z$K],beta)),2,max)
}
cat("abscor:",abscor)
############################################
# Setting 3:
# y = (x1+x2+x3)/(0.5+(x4+x5+1.5)^2) + 0.1*epislon
############################################
  n = 200
  p = 150
  X = obseX(n=n,p=p)
  
  beta1 <- rep(0,p)
  beta2 <- rep(0,p)
  beta1[1:3] <- 1
  beta2[4:5] <- 1
  beta <- beta1+beta2	
  y <- (X[,1]+X[,2]+X[,3])/(0.5+(X[,4]+X[,5]+1.5)^2)+0.1*rnorm(n,0,0.1)
  lambdas=seq(from=0.02,to=0.5,by = 0.02)
  
  t1 = proc.time()
  z = simul(X=X,y=y,n=n,p=p,lambdas=lambdas,beta=beta,seed=2)
  t2 = proc.time()
  cat("time consuming",t2-t1,"s")
  print(z)
  if(z$K==1){
    abscor <-abs(cor((z$temp$vectors[,1:z$K]),cbind(beta1,beta2)))
    abscor <- mean(abscor)
  }
  if(z$K>1){
    abscor <- apply(abs(cor(z$temp$vectors[,1:z$K],cbind(beta1,beta2))),2,max)
    abscor <- mean(abscor)
  }
  cat("abscor:",abscor)



############################################
# Plot_fuction: Plot the subspace distance
# Input: s, n1,n2 is a vector of length 8;
# V1,V2 is a matrix contain the predict beta (dim:8*d)
############################################
Plot_dist = function(s,d1 = 100, d2 =200,n1,n2 ,V1_0, V1_hat, V2_0,V2_hat){
    x1 = s* (log(d1)/n1)^(1/2)
    y1 = c(1:8)
    for(i in 1:8){
      y1[i] = sum(sqrt(abs(V1_0 - V1_hat[i,])^2))  # F-norm
    }
    
    par(pin = c(3.2,2.7))
    plot(x1,y1, xlab = "s × (log d /n)^(1/2)",ylab="Subspace distance",
         xlim = c(0,0.20),ylim = c(0.15,0.50),type = 'l',lwd = 2)
    points(x1,y1,pch= 21)
    
    x2 = s* (log(d2)/n2)^(1/2)
    y2 = c(1:8)
    for(i in 1:8){
      y2[i] = sum(sqrt(abs(V2_0 - V2_hat[i,])^2))  # F-norm
    }
    lines(x2,y2, type = 'l',lwd = 2, col ='gray')
    points(x2,y2,pch= 21, col='gray')
  }
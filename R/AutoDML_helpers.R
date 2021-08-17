two.norm <- function(x){
  return(sqrt(x %*% x))
} 

m<-function(y,d,z,gamma){ #all data arguments to make interchangeable with m2
  return(gamma(1,z)-gamma(0,z))
}

m2<-function(y,d,z,gamma){
  return(y*gamma(d,z))
}

psi_tilde<-function(y,d,z,m,alpha,gamma){
  
  return(m(y,d,z,gamma)+alpha(d,z)*(y-gamma(d,z)))
}

psi_tilde_bias<-function(y,d,z,m,alpha,gamma){
  return(m(y,d,z,gamma))
}

get_MNG<-function(Y,D,X,b){
  
  p=length(b(D[1],X[1,]))
  n.nl=length(D)
  
  B=matrix(0,n.nl,p)
  M=matrix(0,p,n.nl)
  N=matrix(0,p,n.nl)
  
  for (i in 1:n.nl){
    B[i,]=b(D[i],X[i,])
    M[,i]=m(Y[i],D[i],X[i,],b)
    N[,i]=m2(Y[i],D[i],X[i,],b)  # this is a more general formulation for N
  }
  
  M_hat=rowMeans(M)
  N_hat=rowMeans(N)
  G_hat=t(B)%*%B/n.nl
  
  return(list(M_hat,N_hat,G_hat,B))
}


get_D <- function(Y,D,X,m,rho_hat,b){
  n=nrow(X)
  p=length(b(D[1],X[1,]))
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=b(D[i],X[i,])*as.vector(rho_hat %*% b(D[i],X[i,]))-m(Y[i],D[i],X[i,],b)
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

RMD_stable<-function(Y,D,X,p0,D_LB,D_add,max_iter,b,is_alpha,is_lasso){
  
  k=1
  l=0.1
  
  p=length(b(D[1],X[1,]))
  n=length(D)
  
  # low-dimensional moments
  X0=X[,1:p0]
  MNG0<-get_MNG(Y,D,X0,b)
  M_hat0=MNG0[[1]]
  N_hat0=MNG0[[2]]
  G_hat0=MNG0[[3]]
  
  # initial estimate
  rho_hat0=solve(G_hat0,M_hat0)
  rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
  beta_hat0=solve(G_hat0,N_hat0)
  beta_hat=c(beta_hat0,rep(0,p-ncol(G_hat0)))
  
  # moments
  MNG<-get_MNG(Y,D,X,b)
  M_hat=MNG[[1]]
  N_hat=MNG[[2]]
  G_hat=MNG[[3]]
  
  c=0.5
  alpha=0.1
  tol=1e-6
  
  # penalty
  lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) # snippet
  
  ###########
  # alpha_hat
  ###########
  diff_rho=1
  while(diff_rho>tol & k<=max_iter){
    
    # previous values
    rho_hat_old=rho_hat+0
    
    # normalization
    D_hat_rho=get_D(Y,D,X,m,rho_hat_old,b)
    D_hat_rho=pmax(D_LB,D_hat_rho)
    D_hat_rho=D_hat_rho+D_add
    
    # RMD estimate
    #rho_hat=RMD_lasso(M_hat, G_hat, D_hat_rho, lambda)$coefficients
    
    L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,...)
    lambda_vec=lambda*L*D_hat_rho #v3: insert D here
    rho_hat = LassoShooting.fit(G_hat, M_hat, lambda_vec, XX = -G_hat/2, Xy = -M_hat/2,beta.start = rep(0,p))$coefficients
    # difference
    diff_rho=two.norm(rho_hat-rho_hat_old)
    k=k+1
    
  }
  
  return(rho_hat)
}



printer<-function(spec1){
  print(paste(" treated: ",spec1[1], " untreated: ", spec1[2], "   ATE:    ",round(spec1[3],2), "   SE:   ", round(spec1[4],2), sep=""))
}

for_tex<-function(spec1){
  print(paste(" & ",spec1[1], " & ", spec1[2], "   &    ",round(spec1[3],2), "   &   ", round(spec1[4],2), sep=""))
}






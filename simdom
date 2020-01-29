#Comuputing Dominance effects using QMSim outputs
simdom <- function(QTL,qtlmap, var_a, var_d,var_eps, var_res,qtlef,domef, p1, nDQTL, npairs, eps_ef){
  nRow <-dim (QTL)[1]
  nCol <-dim (QTL)[2]
  
  
  G <- matrix (0, nRow, nCol)
  D <- matrix (0, nRow, nCol)
  
  
  arg1 = as.numeric (qtlef[2])  
  arg2 <- as.numeric(qtlef[3])
  arg3 = as.numeric (domef[2])  
  arg4 <- as.numeric(domef[3])
  
  QTL <- QTL-1
  
  qtlmap$a <- NA
  qtlmap$d <- NA
  qtlmap$m <- NA
  if (qtlef[1] == "gamma"){
    for (i in 1:dim(QTL)[2]){
      qtlmap$a[i] <- rgamma (1,shape = arg1, rate=arg2)*sample(c(-1,1),1)
    }
  }else{
    if (qtlef[1] == "normal"){
      for (i in 1:dim(QTL)[2]){
        qtlmap$a[i] <- rnorm (1, arg1, arg2)
      }
      
    }else{
      if (qtlef[1] == "unif"){
        for (i in 1:dim(QTL)[2]){
          qtlmap$a[i] <- runif (1,min = arg1, max=arg2)}
        
      }else{
        cat('\n',"ERROR: Check your distribution parameters")
      }  
    }
  }
  
  if (domef[1] == "gamma"){
    for (i in 1:dim(QTL)[2]){
      qtlmap$d[i] <- rgamma (1,shape = arg3, rate=arg4)*abs(qtlmap$a[i])
    }
  }else{
    if (domef[1] == "normal"){
      for (i in 1:dim(QTL)[2]){
        qtlmap$d[i] <- rnorm (1, arg3, arg4)*abs(qtlmap$a[i])}
      
    }else{
      if (domef[1] == "unif"){
        for (i in 1:dim(QTL)[2]){
          qtlmap$d[i] <- runif (1,min = arg3, max=arg4)*abs(qtlmap$a[i])}
        
      }else{
        cat('\n',"ERROR: Check your Dominance distribution parameters")
      }  
    }
  }
  
  if (var_d==0){
    qtlmap$d <- 0
  }else{ 
    qtlmap$d[as.integer(sample(1:dim(QTL)[2],dim(QTL)[2]-nDQTL))] <-0
  }
  
  
  for (i in 1:nCol){
    p <- p1[i]
    q <- 1-p
    if(var_d==0){
      qtlmap$m[i] <- qtlmap$a[i]
    }else{
      qtlmap$m[i] <- qtlmap$a[i] + (q-p)*qtlmap$d[i]
    }
  }
  
  #Falconer & Mackay
  
  for (j in 1:nCol){
    for (i in 1:nRow){  
      p <- p1[j]
      q <- 1-p
      m <- qtlmap$m[j]
      d <- qtlmap$d[j]
      if (QTL[i,j]==1){
        G[i,j] <- 2*q*m
        D[i,j] <- -2*q^2*d
      }else{
        if (QTL[i,j]==0){
          G[i,j] <- (q-p)*m
          D[i,j] <- 2*p*q*d
        }else{
          if(QTL[i,j]==-1){
            G[i,j] <- -2*p*m
            D[i,j] <- -2*p^2*d
          }
        }
      }
    }
    print (j)
  }
  
  u <- as.matrix (rep(1, nCol),nrow=nCol, ncol=1)
  
  ebv <- G%*%u
  dom <- D%*%u
  
  
  cat('\n',"Breeding values and Dominance deviations: Done...")
  cat('\n',".")
  cat('\n',".")
  cat('\n',".")
  cat('\n',".")
  cat('\n',"Initializing standartizations...")
  
  if (var_d!=0){
    qtlmap$d <- sqrt(var_d)*qtlmap$d/sqrt((4*sum((p1*(1-p1)*qtlmap$d)^2)))
    qtlmap$a <- (sqrt(var_a)*qtlmap$m - (sqrt(2*sum(p1*(1-p1)*qtlmap$m^2))*(1-2*p1)*qtlmap$d))/sqrt(2*sum(p1*(1-p1)*qtlmap$m^2))
  }else{
    qtlmap$a <- (sqrt(var_a)*qtlmap$m) /sqrt(2*sum(p1*(1-p1)*qtlmap$m^2))  
    qtlmap$m <- qtlmap$a
  }
  
  for (i in 1:nCol){
    p <- p1[i]
    q <- 1-p
    qtlmap$m[i] <- qtlmap$a[i] + (q-p)*qtlmap$d[i]
  }                                                                                  
  
  
  
  #Standartized breeding values
  for (j in 1:nCol){
    for (i in 1:nRow){  
      p <- p1[j]
      q <- 1-p
      m <- qtlmap$m[j]
      d <- qtlmap$d[j]
      if (QTL[i,j]==1){
        G[i,j] <- 2*q*m
        D[i,j] <- -2*q^2*d
      }else{
        if (QTL[i,j]==0){
          G[i,j] <- (q-p)*m
          D[i,j] <- 2*p*q*d
        }else{
          if(QTL[i,j]==-1){
            G[i,j] <- -2*p*m
            D[i,j] <- -2*p^2*d
          }
        }
      }
    }
    print(j)
  }
  
  u <- as.matrix (rep(1, nCol),nrow=nCol, ncol=1)
  ebv <- G%*%u
  dom <- D%*%u
  cat('\n',"Done...")
  
  
  #Epistasis effects using a multiplicative model  
  if (eps_ef==T){
    
    E <- matrix (0, nRow, npairs)
    
    lociep <- matrix (0,2,npairs)
    lociep[1,]<- sample(1:nCol,npairs,replace = T)
    lociep[2,]<- sample(1:nCol,npairs,replace = T)
    
    eps_eff =NULL
    
    for(i in 1:npairs){
      eps_eff[i] <- rnorm(1,0,1)#qtlmap$m[lociep[1,i]]*qtlmap$m[lociep[2,i]]
    }
    
    for (i in 1:npairs){
      E[,i] <- G[,lociep[1,i]]*G[,lociep[2,i]]*eps_eff[i]
    }
    
    uz <-  as.matrix (rep(1, npairs),nrow=npairs, ncol=1)
    
    eps <- E%*%uz
    ve <- NULL
    for (i in 1:npairs){
      ve[i] <- 4*p1[lociep[1,i]]*(1-p1[lociep[1,i]])*p1[lociep[2,i]]*(1-p1[lociep[2,i]])*eps_eff[i]^2
    }
    
    vep<- sum(ve)
    
    eps_eff <- (sqrt(var_eps)*eps_eff)/sqrt(vep)
    
    for (i in 1:npairs){
      E[,i] <- G[,lociep[1,i]]*G[,lociep[2,i]]*eps_eff[i]
    }
    
    eps <- E%*%uz
    
  }
  
  res <- rnorm (nRow, 0, sqrt(var_res))
  
  
  if (eps_ef==T){
    y <- ebv + dom + eps + res 
  }else{
    y <- ebv + dom + res  
  }
  
  if(eps_ef==T){
    pheno <- data.frame (ebv=ebv, dom=dom, eps=eps,res=res, y=y)
  }else{pheno <- data.frame (ebv=ebv, dom=dom, res=res, y=y)}
  if(eps_ef==T){  
    out <- list (pheno, qtlmap, lociep)
  }else{out <- list (pheno, qtlmap)}
  #pheno <- data.frame (ebv=ebv, dom=dom, res=res, y=y)
  #out <- list (pheno, qtlmap)
  return (out)
}

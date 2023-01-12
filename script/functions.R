require(betafunctions)


pa<-function(datas){
  
  data<-2-datas
  
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  
  l <- rep(list(factor(1:2,levels=c(1,2),labels=c(1,2))), nrat)
  pi.label<-expand.grid(l)
  
  data.f_<-as.data.frame.table(table(data)/n)
  names(pi.label)<-names(data.f_)[1:nrat]
  all.pi<-as.data.frame(pi.label)
  
  data.f<-merge(data.f_,as.data.frame(pi.label),all.y=TRUE)
  data.f[is.na(data.f)] <- 0
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  l <- rep(list(1:2), nrat)
  pi.label<-expand.grid(l)
  
  d<-matrix(0,nrow=3*npair,ncol=nrow(pi.label))
  
  for (i in 1:npair){
    d[(i-1)*3+1,pi.label[,rat.pair[i,1]]==1 & pi.label[,rat.pair[i,2]]==1]<-1
    d[(i-1)*3+2,pi.label[,rat.pair[i,1]]==1 & pi.label[,rat.pair[i,2]]==2]<-1
    d[(i-1)*3+3,pi.label[,rat.pair[i,1]]==2 & pi.label[,rat.pair[i,2]]==1]<-1
  }
  
  
  
  var_cov_pi<-diag(data.f[,nrat+1])-matrix(data.f[,nrat+1],ncol=1)%*%t(data.f[,nrat+1])
  
  var_cov_pair<-d%*%var_cov_pi%*%t(d)
  
  d2<-matrix(NA,nrow=2,ncol=3*npair)
  d2[1,]<-rep(c(2,0,0),npair)
  d2[2,]<-rep(c(2,1,1),npair)
  
  var_cov_two<-d2%*%var_cov_pair%*%t(d2)
  
  a<-rep(NA,npair)
  b<-rep(NA,npair)
  
  for (i in 1:npair){
    a[i]<-sum(data.f[data.f[,rat.pair[i,1]]==1 & data.f[,rat.pair[i,2]]==1,(nrat+1)])
    b[i]<-sum(2*data.f[data.f[,rat.pair[i,1]]==1 & data.f[,rat.pair[i,2]]==1,(nrat+1)]+data.f[data.f[,rat.pair[i,1]]==1 & data.f[,rat.pair[i,2]]==2,(nrat+1)]
              +data.f[data.f[,rat.pair[i,1]]==2 & data.f[,rat.pair[i,2]]==1,(nrat+1)])
    
  }
  
  A<-2*sum(a)
  
  B<-sum(b)
  
  ppos<-A/B
  
  d3<-matrix(c(1/B,-A/B^2),nrow=1)
  
  var.ppos<-d3%*%var_cov_two%*%t(d3)/n
  
  result<-matrix(c(ppos,sqrt(var.ppos)),nrow=1)
  colnames(result) <- c("PPOS", "SE(PPOS)")
  rownames(result) <- ""
  return(result)
}





pca<-function(datas){
  
  data<-2-datas
  
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  
  # nij<-cbind(rowSums(data+1),nrat-rowSums(data+1))
  # ppos<-sum(nij[,2]*(nij[,2]-1))/(sum(nij[,2]*(nij[,2]-1))/2+n*npair-sum(nij[,1]*(nij[,1]-1))/2)
  
  l <- rep(list(factor(1:2,levels=c(1,2),labels=c(1,2))), nrat)
  pi.label<-expand.grid(l)
  
  data.f_<-as.data.frame.table(table(data)/n)
  names(pi.label)<-names(data.f_)[1:nrat]
  all.pi<-as.data.frame(pi.label)
  
  data.f<-merge(data.f_,as.data.frame(pi.label),all.y=TRUE)
  data.f[is.na(data.f)] <- 0
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  l <- rep(list(1:2), nrat)
  pi.label<-expand.grid(l)
  
  d<-matrix(0,nrow=3*npair,ncol=nrow(pi.label))
  
  for (i in 1:npair){
    d[(i-1)*3+1,pi.label[,rat.pair[i,1]]==1 & pi.label[,rat.pair[i,2]]==1]<-1
    d[(i-1)*3+2,pi.label[,rat.pair[i,1]]==1 & pi.label[,rat.pair[i,2]]==2]<-1
    d[(i-1)*3+3,pi.label[,rat.pair[i,1]]==2 & pi.label[,rat.pair[i,2]]==1]<-1
  }
  
  
  
  var_cov_pi<-diag(data.f[,nrat+1])-matrix(data.f[,nrat+1],ncol=1)%*%t(data.f[,nrat+1])
  
  var_cov_pair<-d%*%var_cov_pi%*%t(d)
  
  d2<-matrix(NA,nrow=2,ncol=3*npair)
  d2[1,]<-rep(c(2,0,0),npair)
  d2[2,]<-rep(c(2,1,1),npair)
  
  var_cov_two<-d2%*%var_cov_pair%*%t(d2)
  
  a<-rep(NA,npair)
  b<-rep(NA,npair)
  
  for (i in 1:npair){
    a[i]<-sum(data.f[data.f[,rat.pair[i,1]]==1 & data.f[,rat.pair[i,2]]==1,(nrat+1)])
    b[i]<-sum(2*data.f[data.f[,rat.pair[i,1]]==1 & data.f[,rat.pair[i,2]]==1,(nrat+1)]+data.f[data.f[,rat.pair[i,1]]==1 & data.f[,rat.pair[i,2]]==2,(nrat+1)]
              +data.f[data.f[,rat.pair[i,1]]==2 & data.f[,rat.pair[i,2]]==1,(nrat+1)])
    
  }
  
  A<-2*sum(a)
  
  B<-sum(b)
  
  ppos<-A/B
  
  pca<-ppos/(2-ppos)
  
  
  d3<-matrix(c(1/B,-A/B^2),nrow=1)
  
  var.ppos<-d3%*%var_cov_two%*%t(d3)/n
  
  var.pca<-(2/(2-ppos)^2)^2*var.ppos
  
  result<-matrix(c(pca,sqrt(var.pca)),nrow=1)
  colnames(result) <- c("PCA", "SE(PCA)")
  rownames(result) <- ""
  return(result)
}



po<-function(data){
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  Poip<- mapply(agree <- function(x, y) {ifelse(data[, x] == data[, y], 1, 0)},rat.pair[, 2],rat.pair[, 1])
  
  Poi<-rowMeans(Poip)
  
  Po<-mean(Poi)
  
  var.Po<-(sum(Poi^2)-n*Po^2)/n^2
  
  result<-matrix(c(Po,sqrt(var.Po)),nrow=1)
  colnames(result) <- c("Po", "SE(Po)")
  rownames(result) <- ""
  return(result)
  
}

  
po.bayes<-function(datas,quantiles){
  
  n<-nrow(datas)
  nrat<-ncol(datas)
  npair<-nrat*(nrat-1)/2
  number<-1/2^nrat
  
  #Step 1: obtain the frequencies of the 2^R table
  l <- rep(list(factor(1:2,levels=c(1,2),labels=c(1,2))), nrat)
  pi.label<-expand.grid(l)
  
  data.f_<-as.data.frame.table(table(data.frame(2-datas))/n)
  names(pi.label)<-names(data.f_)[1:nrat]
  all.pi<-as.data.frame(pi.label)
  
  data.f<-merge(data.f_,as.data.frame(pi.label),all.y=TRUE)
  data.f[is.na(data.f)] <- 0
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  fre<-data.f[,nrat+1]
  
  
  alpha<-fre*n+number
  #Step 2: determine the weights lambda to calculate po
  
  index<-seq(1,length(fre))
  
  ind_1<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==1])
  ind_2<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==2])
  
  lambda<-as.data.frame.table(table(c(c(ind_1),c(ind_2))))[,2]/npair
  lambdak<-max(lambda)
  lambda1<-min(lambda)
  nu<-sum(alpha)
  
  if (nrat>2){
  El<-sum((alpha*lambda))/nu
  Vl<-(nu*sum(alpha*lambda^2)-sum((alpha*lambda))^2)/(nu^2*(nu+1))
  
  p<-((El-lambda1)/(lambdak-lambda1))*((El-lambda1)*(lambdak-El)/Vl-1)
  q<-((lambdak-El)/(lambdak-lambda1))*((El-lambda1)*(lambdak-El)/Vl-1)
  
  mean.beta4<-lambda1+(lambdak-lambda1)*p/(p+q)
  
  var.beta4<-p*q*(lambdak-lambda1)^2/((p+q)^2*(p+q+1))
  
  qua.beta4<-qBeta.4P(quantiles, lambda1, lambdak, p, q, lower.tail = TRUE)
  
  
  result<-matrix(c(mean.beta4,sqrt(var.beta4),qua.beta4),nrow=1)
  colnames(result) <- c("Po", "SE(Po)",rep(NA,length(qua.beta4)))
  rownames(result) <- ""
  return(result)
  }
  if (nrat==2){
    p<-sum((alpha*c(1,0,0,1)))
     
    q<-nu-p
    
    mean.beta<-p/(p+q)
    
    var.beta<-p*q/((p+q)^2*(p+q+1))
    
    qua.beta<-qbeta(quantiles, p, q, lower.tail = TRUE)
    
    
    result<-matrix(c(mean.beta,sqrt(var.beta),qua.beta),nrow=1)
    colnames(result) <- c("Po", "SE(Po)",rep(NA,length(qua.beta)))
    rownames(result) <- ""
    return(result)
    
  }
}


RA<-c(0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,0)
RB<-c(0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,1,0,0,0)
RC<-c(0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1,0,1,0)
RD<-c(0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,1,0,1,0)
RE<-c(0,0,0,0,1,0,1,1,1,1,1,1,1,1,0,1,1,1,0,0)

grant<-data.frame(cbind(RA,RB,RC,RD,RE))

datas<-grant[,1:2]



pa.bayes<-function(data,size,quantiles){
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  number<-1/2^nrat
  
  l <- rep(list(factor(1:2,levels=c(1,2),labels=c(1,2))), nrat)
  pi.label<-expand.grid(l)
  
  data.f_<-as.data.frame.table(table(data.frame(2-data))/n)
  names(pi.label)<-names(data.f_)[1:nrat]
  all.pi<-as.data.frame(pi.label)
  
  data.f<-merge(data.f_,as.data.frame(pi.label),all.y=TRUE)
  data.f[is.na(data.f)] <- 0
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  fre<-data.f[,nrat+1]
  
  index<-seq(1,length(fre))
  
  ind_1<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==1])
  ind_2<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==2])
  ind_3<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==2])
  ind_4<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==1])
  
  N2 <-size
  
  alpha<-fre*n+number
  
  w<- lapply(alpha, function(x) rgamma(N2,shape=x,rate=1))
  
  t<-Reduce("+", w)
  
  z<-sapply(w, function(x) x/t)
  
  if (npair>1){
    p11<-rowSums(z[,c(ind_1)])/(npair)
    p22<-rowSums(z[,c(ind_2)])/(npair)
    p12<-rowSums(z[,c(ind_3)])/(npair)
    p21<-rowSums(z[,c(ind_4)])/(npair)
    
  }
  
  if (npair==1){
    p11<-z[,c(ind_1)]/(npair)
    p22<-z[,c(ind_2)]/(npair)
    p12<-z[,c(ind_3)]/(npair)
    p21<-z[,c(ind_4)]/(npair)
    
  }
  pa<-(2*p11)/(p11+p12+p21+p11)
  
  mean.pa<-mean(pa)
  
  var.pa<-var(pa)
  
  q.pa<-quantile(pa,quantiles)
  
  result<-matrix(c(mean.pa,sqrt(var.pa),q.pa),nrow=1)
  colnames(result) <- c("Pa", "SE(Pa)",rep(NA,length(q.pa)))
  rownames(result) <- ""
  return(result)
  
}




pca.bayes<-function(data,size,quantiles){
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  number<-1/2^nrat
  
  l <- rep(list(factor(1:2,levels=c(1,2),labels=c(1,2))), nrat)
  pi.label<-expand.grid(l)
  
  data.f_<-as.data.frame.table(table(data.frame(2-data))/n)
  names(pi.label)<-names(data.f_)[1:nrat]
  all.pi<-as.data.frame(pi.label)
  
  data.f<-merge(data.f_,as.data.frame(pi.label),all.y=TRUE)
  data.f[is.na(data.f)] <- 0
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  fre<-data.f[,nrat+1]
  
  index<-seq(1,length(fre))
  
  ind_1<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==1])
  ind_2<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==2])
  ind_3<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==2])
  ind_4<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==1])
  
  N2 <-size
  
  alpha<-fre*n+number
  
  w<- lapply(alpha, function(x) rgamma(N2,shape=x,rate=1))
  
  t<-Reduce("+", w)
  
  z<-sapply(w, function(x) x/t)
  
  if (npair>1){
    p11<-rowSums(z[,c(ind_1)])/(npair)
    p22<-rowSums(z[,c(ind_2)])/(npair)
    p12<-rowSums(z[,c(ind_3)])/(npair)
    p21<-rowSums(z[,c(ind_4)])/(npair)
    
  }
  
  if (npair==1){
    p11<-z[,c(ind_1)]/(npair)
    p22<-z[,c(ind_2)]/(npair)
    p12<-z[,c(ind_3)]/(npair)
    p21<-z[,c(ind_4)]/(npair)
    
  }
  pca<-(p11)/(p11+p12+p21)
  
  mean.pca<-mean(pca)
  
  var.pca<-var(pca)
  
  q.pca<-quantile(pca,quantiles)
  
  result<-matrix(c(mean.pca,sqrt(var.pca),q.pca),nrow=1)
  colnames(result) <- c("Pca", "SE(Pca)",rep(NA,length(q.pca)))
  rownames(result) <- ""
  return(result)
  
}


agree.plot<-function(r1,r2){
  
  table22<-addmargins(prop.table(table(r1,r2)))*100
  
  plot(1,1,xlim=c(0,100),ylim=c(0,100),type="n",xlab="",ylab="rater 1",xaxs = "i", yaxs = "i",xaxt = "n",yaxt="n")
  axis(3,at=seq(0,100,20))
  axis(2,at=seq(0,100,20),lab=seq(100,0,-20))
  mtext(text = "Rater 2",side = 3,line = 2)
  segments(0,table22[2,3],100,table22[2,3])
  segments(100-table22[3,2],0,100-table22[3,2],100)
  rect(0,100,table22[1,1],100-table22[1,1],col=1)
  rect(100-table22[2,2],table22[2,2],100,0,col=1)
  segments(0,100,100,0,col=2)
}



correlo1<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
for (i in 1:nrat){
  for (j in 1:nrat){
 
    data_temp<-as.data.frame(data[,c(i,j)])
    
      result[i,j]<-po(data_temp)[1]
    }
  }

colnames(result)<-colnames(data)
rownames(result)<-colnames(data)

corrplot(result, method="circle",addCoef.col = rgb(1,1,1, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)

}

correlo2<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-pa(data_temp)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="circle",addCoef.col = rgb(1,1,1, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}



pa_mean<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-pa(data_temp)[1]
    }
  }
  mean(result[upper.tri(result, diag = FALSE)])
 
}

correlo3<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-pca(data_temp)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="circle",addCoef.col = rgb(1,1,1, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}

pca_mean<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-pca(data_temp)[1]
    }
  }
  mean(result[upper.tri(result, diag = FALSE)])
  
}

prop_mean<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-table(as.data.frame(data[,c(i,j)]))
      
      result[i,j]<-data_temp[2,2]+data_temp[1,2]+data_temp[2,1]
    }
  }
 return(result)
  
}
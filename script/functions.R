require(betafunctions)

many1<-function(data.raw){
  datas<-cbind(rowSums(data.raw),ncol(data.raw)-rowSums(data.raw))
  
  #data are counts in each category, not raw data!
  N<-nrow(datas)            #number of subjects
  ncat<-ncol(datas)         #number of categories
  nrat<-rowSums(datas)      #number of raters per subject
  idd2<-seq(1,nrow(datas))  #subject ID
  
  
  #OBSERVED AGREEMENT
  oprim_i<-rowSums(datas*(datas-1))/(nrat*(nrat-1))
  P_o<-mean(oprim_i)
  
  #MARGINAL PROBABILITY DISTRIBUTION
  p_j<-colMeans(sweep(datas,1,nrat,"/"))
  eprim_i<-rowSums(sweep(datas,2,p_j,"*"))/nrat
  
  #EXPECTED AGREEMENT
  P_e1<-mean(eprim_i)		
  
  
  #STANDARD ERROR
  d1<-sum(((1-P_e1)*oprim_i-2*(1-P_o)*eprim_i)^2)/N
  d2<-(P_o*P_e1-2*P_e1+P_o)^2
  var1<-(d1-d2)/(N*(1-P_e1)^4)
  kappa_1<-(P_o-P_e1)/(1-P_e1)
  
  return(c(kappa_1,sqrt(var1)))
  
}




many2<-function(datas){
  

  idd2<-seq(1,nrow(datas))
  ncat<-2
  #CREATE THE COMPLETE CASE DATASET
  data_<-na.omit(cbind(datas,idd2))
  data2<-data_[,1:ncol(datas)]
 
  
  #CREATE FACTORS WITH THE CLASSIFICATION OF THE OBSERVERS (not usefull)
  all<-c(data2)
  all.f<-factor(all,levels=names(table(all)))
  # ncat<-nlevels(all.f)    		#number of categories of the scale
  
  K<-max(idd2)						    #number of clusters
  nrat<-ncol(datas)					#number of raters
  npair<-nrat*(nrat-1)/2  				  #number of ditincts pairs 


  #CREATE WEIGHT FOR THE ITEMS
  xi<-prop.table(table(idd2))

  #CREATE THE OBSERVED AGREEMENT
  tmp = expand.grid(1:nrat[1],1:(nrat[1]-1))
  tmp2<-tmp[tmp[,1]>tmp[,2],]
  RB<-tmp2[,1]#ID number of the rater for pair i
  RA<-tmp2[,2]#ID number of the rater for pair i
  
  a<-mapply(agree<-function(x,y){ifelse(data2[,x]==data2[,y],1,0)},tmp2[,1],tmp2[,2])
  
  if (npair>1){Po<-unname(as.list(aggregate(a, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1]))}
  if (npair==1){Po<-replicate(npair,matrix(0,ncol=K,nrow=1), simplify=FALSE);Po[[1]]<-unname(aggregate(a, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1])}
  
  #CREATE THE MARGINALS AT THE ITEM LEVEL
  p <- replicate(ncol(datas), matrix(0, ncol = max(idd2), nrow = ncat),simplify = FALSE)
  for (i in 1:nrat) {
    for (s in 1:K) {
      p[[i]][, s] <- prop.table(table(factor(data2[, i][idd2==s], levels = c(0,1))))
    }
  }
  
  #CREATE THE MATRIX OMEGA
  Omega<-(diag(K)-xi%*%matrix(1,ncol=K,nrow=1))%*%diag(xi^2)%*%(diag(K)-matrix(1,ncol=1,nrow=K)%*%t(xi))
  
  
  #CREATE THE MATRIX V
  V<-matrix(NA,ncol=(npair+nrat*ncat),nrow=(npair+nrat*ncat))
  
  #var-cov for observed agreement
  tmp = expand.grid(1:npair,1:npair)
  V[1:npair,1:npair]<-mapply(function(i, j) (K^2/(K-1))*matrix(Po[[i]],nrow=1)%*%Omega%*%t(matrix(Po[[j]],nrow=1)),tmp[,1],tmp[,2])
  
  
  #var-cov for observed agreement versus marginals
  tmp = expand.grid(1:nrat,1:npair)
  V[1:npair,(npair+1):(npair+nrat*ncat)]<-t(matrix(mapply(function(i, j)(K^2/(K-1))*matrix(Po[[i]],nrow=1)%*%Omega%*%t(p[[j]]),tmp[,2],tmp[,1]),ncol=npair))
  
  #var-cov for marginals
  for (i in 1:nrat){
    for (j in (i:nrat)){
      V[(npair+1+(i-1)*ncat):(npair+i*ncat),(npair+1+(j-1)*ncat):(npair+j*ncat)]<-(K^2/(K-1))*p[[i]]%*%Omega%*%t(p[[j]])
    }
  }
  
  for (i in (2:(npair+nrat*ncat))){
    for (j in 1:(i-1)){
      V[i,j]<-V[j,i]
    }
  }
  
  
  #CREATE THE MARGINALS OVERALL
  M<-replicate(nrat,matrix(0,ncol=ncat,nrow=1), simplify=FALSE)
  M<-lapply(1:nrat,f1<-function(idd){prop.table(table(factor(data2[,idd],levels=c(0,1))))})
  
  
  
  #CREATE THE MATRIX J
  J<-matrix(0,2*npair,(npair+nrat*ncat))
  
  if (npair>1){
    diag(J[1:npair,1:npair])<-1
  }
  if (npair==1){
    J[1:npair,1:npair]<-1
  }
  
  i<-0
  for (i in (0:(npair-1))){
    
    J[(npair+i+1),(npair+1+(RA[i+1]-1)*ncat):(npair+RA[i+1]*ncat)]<-t(M[[RB[i+1]]])
    J[(npair+i+1),(npair+1+(RB[i+1]-1)*ncat):(npair+RB[i+1]*ncat)]<-t(M[[RA[i+1]]])                                                         
    
  }
  
  #COMPUTE VARPSY
  varpsy<-J%*%V%*%t(J)
  
  #COMPUTE THE VECTOR OF EXPECTED AGREEMENT
  
  Pe<-matrix(NA,ncol=1,nrow=npair)
  P0<-matrix(NA,ncol=1,nrow=npair)
  
  for (i in 0:(npair-1)){
    
    P0[i+1]<-sum(Po[[i+1]]*matrix(xi,ncol=K,nrow=1))               
    Pe[i+1]<-t(M[[RA[i+1]]])%*%M[[RB[i+1]]]
    
  }
  
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR THE GENERAL Po and Pe
  
  L<-matrix(0,2,2*npair)
  
  L[1,1:npair]<-1/npair
  L[2,(npair+1):(2*npair)]<-1/npair
  
  
  var_bkappa<-L%*%varpsy%*%t(L)
  
  
  
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR THE KAPPAS
  
  s<-matrix(c(1/(1-mean(Pe)),(mean(P0)-1)/(1-mean(Pe))^2),ncol=2)
  
  var_kappa<-s%*%var_bkappa%*%t(s)*(K-1)/K^2#correction factor for one unit per cluster
  
  
  kappa<-matrix((mean(P0)-mean(Pe))/(1-mean(Pe)),ncol=1)
  
  return(c(kappa,sqrt(var_kappa)))
  
  
}

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

correlo4<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-many1(data_temp)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="circle",addCoef.col = rgb(1,1,1, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}

correlo5<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-many2(data_temp)[1]
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





size.pom<-function(a.level=0.05,w=0.20,pm=0.5,po=0.8,nrat=2,e_po=0,e_pm=0,nb_s=40000,quant=100){
  
  min_po<-1-2*min(pm,1-pm)
  
  if (a.level > 1 | a.level < 0)
    stop("the significance level should be between 0 and 1")
  if (po> 1 |po<min_po)
  {writeLines(paste("Stop: po should be less than ",1));stop(paste("Stop: po should be more than ",round(min_po,2)))}
  
  z<-qnorm(1-a.level/2)
  
  if (nrat==2){
    N<-4*z^2*po*(1-po)/w^2
   
    return(ceiling(N))
    
  }
  
  if (nrat>2){
    
    
      npair<-nrat*(nrat-1)/2
      r<-seq(0,nrat,1)
      agree_2<-(r*(r-1)+(nrat-r)*(nrat-r-1))/2
      
      E<-matrix(ncol=nrat+1,byrow=TRUE,data=c(agree_2,r,rep(1,nrat+1)))
      F<-c(po*npair,pm*nrat,1)
      G<-diag(nrat+1)
      H<-rep(0,nrat+1)
      
      library(limSolve)
      scenario<- round(xsample(E = E, F = F, G = G, H = H, type = "mirror",iter=nb_s, burninlength = 1000)$X,2)
      
      #l <- rep(list(seq(0,1,0.1)), nrat+1)
      #scenario<-expand.grid(l)[rowSums(expand.grid(l))==1,]
      sumpm<-rowSums((1/nrat)*sweep(scenario, MARGIN=2, r, `*`))
      sumpo<-rowSums((1/npair)*sweep(scenario, MARGIN=2,agree_2, `*`))
      valid<-scenario[sumpo<=po+e_po & sumpo>=po-e_po & sumpm<=pm+e_pm & sumpm>=pm-e_pm,]
      po_r<-round(sumpo[sumpo<=po+e_po & sumpo>=po-e_po & sumpm<=pm+e_pm & sumpm>=pm-e_pm],2)
      pm_r<-round(sumpm[sumpo<=po+e_po & sumpo>=po-e_po & sumpm<=pm+e_pm & sumpm>=pm-e_pm],2)
      if (nrow(valid)<1)
        stop("There is no valid scenario, please provide some positive value for tolerance")
      var_po<-(rowSums(sweep(valid, MARGIN=2, (agree_2/npair)^2, `*`)))-po_r^2
      
      N<-ceiling(4*z^2*var_po/w^2)
      #given<-matrix(rep(c(a.level,w,nrat),each=nrow(valid)),ncol=3)
      #result<-data.frame(unique(cbind(given,matrix(pm_r,ncol=1),matrix(po_r,ncol=1),matrix(N,ncol=1),valid)))
      #colnames(result) <- c("alpha level", "width", "nb raters","pm", "po","min N",colnames(valid))
      #rownames(result) <- rep("",nrow(result))
      quantil<-quantile(N, probs = quant/100)
      return(ceiling(quantil))
      #xranges(E, F, G = G, H = H)
  }
}




# 
# w<-0.1
# pm<-0.5
# po<-0.8
# a.level<-0.05
# e_po<-0
# e_pm<-0
# quant<-100
# 
# size.pom(w=0.1,pm=0.5,po=0.8,a.level=0.05,e_po=0,e_pm=0,quant=100,nrat=2,nb_s=40000)

size.pam<-function(a.level=0.05,w=0.20,pm=0.5,pa=0.8,nrat=2,e_pa=0,e_pm=0,nb_s=40000,quant=100){
  
  min_p11<-pm-min(pm,1-pm)
  max_p11<-pm
  
  min_pa<-2*min_p11/(pm+pm)
  max_pa<-2*max_p11/(pm+pm)
  
  if (a.level > 1 | a.level < 0) 
    stop("the significance level should be between 0 and 1")
  if (pa> max_pa |pa<min_pa) 
  {writeLines(paste("Stop: pa should be less than ",round(max_pa,2)));stop(paste("Stop: pa should be more than ",round(min_pa,2)))}
  
  z<-qnorm(1-a.level/2)
  
  if (nrat==2){
    N<-4*z^2*pa*(1-pa)*(2-pa)/(w^2*(pm+pm))
   return(ceiling(N))
    
  }
  
  if (nrat>2){
    
   
      npair<-nrat*(nrat-1)/2
      r<-seq(0,nrat,1)
      
      agree1_2<-(r*(r-1))/2             #number of pairs agreeing on cat1
      agree1_3<-(r*(r-1)*(r-2))/6       #number of triplets agreeing on cat1
      agree1_4<-(r*(r-1)*(r-2)*(r-3))/24#number of quartets agreeing on cat1
      
      
      E<-matrix(ncol=nrat+1,byrow=TRUE,data=c(agree1_2,r,rep(1,nrat+1)))
      F<-c(pm*pa*npair,pm*nrat,1)
      G<-diag(nrat+1)
      H<-rep(0,nrat+1)
      
      library(limSolve)
      scenario<- round(xsample(E = E, F = F, G = G, H = H, type = "mirror",iter=nb_s, burninlength = 1000)$X,2)
      
      
      A<-2*rowSums(sweep(scenario, MARGIN=2, agree1_2, `*`))
      B<-(nrat-1)*rowSums(sweep(scenario, MARGIN=2, r, `*`))
      
      pa_expected<-A/B
      
      valid<-scenario[pa_expected<=pa+e_pa & pa_expected>=pa-e_pa & B/(nrat*(nrat-1))<=pm+e_pm & B/(nrat*(nrat-1))>=pm-e_pm,]
      
      if (nrow(valid)<1) 
        stop("There is no valid scenario, please provide some positive value for tolerance")
      
      B_r<-(nrat-1)*rowSums(sweep(valid, MARGIN=2, r, `*`))
      A_r<-2*rowSums(sweep(valid, MARGIN=2, agree1_2, `*`))
      
      pa_r<-A_r/B_r
      
      p3_r<-rowSums(sweep(valid, MARGIN=2, agree1_3, `*`))
      p4_r<-rowSums(sweep(valid, MARGIN=2, agree1_4, `*`))
      
      var_pa<-((2-(nrat-1)*pa_r)*(pa_r*(1-(nrat-1)*pa_r)+12*p3_r/B_r)+24*p4_r/B_r)/B_r
      
      N<-ceiling(4*z^2*var_pa/w^2)
      
      # given<-matrix(rep(c(a.level,w,nrat),each=nrow(valid)),ncol=3)
      #result<-data.frame(unique(cbind(given,matrix(B_r/((nrat-1)*4),ncol=1),matrix(pa_r,ncol=1),matrix(N,ncol=1),valid)))
      #colnames(result) <- c("alpha level", "width", "nb raters","pm", "pa","min N",colnames(valid))
      #rownames(result) <- rep("",nrow(result))
      quantil<-quantile(N, probs = quant/100)
      
     
      return(ceiling(quantil))
    
    
    
  }
}


size.kid<-function(a.level=0.05,w=0.20,pm=0.5,ki=0.8,nrat=2,e_ki=0,e_pm=0,scenario=NULL,nb_s=40000,quant=100){
  
  min_po<-1-2*min(pm,1-pm)
  min_ki<-(min_po-(pm^2+(1-pm)^2))/(1-(pm^2+(1-pm)^2))
  
  if (a.level > 1 | a.level < 0)
    stop("the significance level should be between 0 and 1")
  if (ki> 1 |ki<min_ki)
  {writeLines(paste("Stop: po should be less than ",1));stop(paste("Stop: po should be more than ",round(min_ki,2)))}
  
  z<-qnorm(1-a.level/2)
  
  
  if (nrat>=2){
    
    if (is.null(scenario)){
      npair<-nrat*(nrat-1)/2
      r<-seq(0,nrat,1)
      agree_2<-(r*(r-1)+(nrat-r)*(nrat-r-1))/2
      pexp<-pm^2+(1-pm)^2
      
      E<-matrix(ncol=nrat+1,byrow=TRUE,data=c(agree_2,r,rep(1,nrat+1)))
      F<-c((ki*(1-pexp)+pexp)*(npair),pm*nrat,1)
      G<-diag(nrat+1)
      H<-rep(0,nrat+1)
      
      library(limSolve)
      scenario<- round(xsample(E = E, F = F, G = G, H = H, type = "mirror",iter=nb_s, burninlength = 1000)$X,2)
      
      #l <- rep(list(seq(0,1,0.1)), nrat+1)
      #scenario<-expand.grid(l)[rowSums(expand.grid(l))==1,]
      sumpm<-rowSums((1/nrat)*sweep(scenario, MARGIN=2, r, `*`))
      sumpo<-rowSums((1/npair)*sweep(scenario, MARGIN=2,agree_2, `*`))
      sumpe<-sumpm^2+(1-sumpm)^2
      sumki<-(sumpo-sumpe)/(1-sumpe)
      
      
      valid<-scenario[sumki<=ki+e_ki & sumki>=ki-e_ki & sumpm<=pm+e_pm & sumpm>=pm-e_pm,]
      
      ki_r<-round(sumki[sumki<=ki+e_ki & sumki>=ki-e_ki & sumpm<=pm+e_pm & sumpm>=pm-e_pm],2)
      pm_r<-round(sumpm[sumki<=ki+e_ki & sumki>=ki-e_ki & sumpm<=pm+e_pm & sumpm>=pm-e_pm],2)
      po_r<-round(sumpo[sumki<=ki+e_ki & sumki>=ki-e_ki & sumpm<=pm+e_pm & sumpm>=pm-e_pm],2)
      pe_r<-round(sumpe[sumki<=ki+e_ki & sumki>=ki-e_ki & sumpm<=pm+e_pm & sumpm>=pm-e_pm],2)
      
      sumnr1<-rowSums(sweep(valid, MARGIN=2,r, `*`))
      sumnr2<-rowSums(sweep(valid, MARGIN=2,r^2, `*`))
      sumnr3<-rowSums(sweep(valid, MARGIN=2,r^3, `*`))
      sumnr4<-rowSums(sweep(valid, MARGIN=2,r^4, `*`))
      
      
      poi_2<-(4*sumnr4-8*nrat*sumnr3+4*nrat*(2*nrat-1)*sumnr2-4*nrat^2*(nrat-1)*sumnr1+nrat^2*(nrat-1)^2)/(nrat^2*(nrat-1)^2)
      
      pei_2<-((2*pm_r-1)^2*sumnr2+2*nrat*(1-pm_r)*(2*pm_r-1)*sumnr1+nrat^2*(1-pm_r)^2)/(nrat^2)
      
      popei<-(2*(2*pm_r-1)*sumnr3+2*nrat*(2-3*pm_r)*sumnr2+nrat*(4*nrat*pm_r-3*nrat-2*pm_r+1)*sumnr1+nrat^2*(nrat-1)*(1-pm_r))/(nrat^2*(nrat-1))
      
      
      if (nrow(valid)<1)
        stop("There is no valid scenario, please provide some positive value for e_ki or e_pm")
      
      term2<-(po_r*pe_r-2*pe_r+po_r)^2
      
      term11<-(1-pe_r)^2*poi_2
      
      term12<--4*(1-pe_r)*(1-po_r)*popei
      
      term13<-4*(1-po_r)^2*pei_2
      
      
      var_ki<-(term11+term12+term13-term2)/((1-pe_r)^4)
      
      N<-ceiling(4*z^2*var_ki/w^2)
      #given<-matrix(rep(c(a.level,w,nrat),each=nrow(valid)),ncol=3)
      #result<-data.frame(unique(cbind(given,matrix(pm_r,ncol=1),matrix(po_r,ncol=1),matrix(N,ncol=1),valid)))
      #colnames(result) <- c("alpha level", "width", "nb raters","pm", "po","min N",colnames(valid))
      #rownames(result) <- rep("",nrow(result))
      quantil<-quantile(N, probs = quant/100)
      
      return(ceiling(quantil))
      #xranges(E, F, G = G, H = H)
    }
    
    if (!is.null(scenario)){
      npair<-nrat*(nrat-1)/2
      r<-seq(0,nrat,1)
      agree_2<-(r*(r-1)+(nrat-r)*(nrat-r-1))/2
      #l <- rep(list(seq(0,1,0.1)), nrat+1)
      #scenario<-expand.grid(l)[rowSums(expand.grid(l))==1,]
      pm_r<-rowSums((1/nrat)*sweep(scenario, MARGIN=2, r, `*`))
      po_r<-rowSums((1/npair)*sweep(scenario, MARGIN=2,agree_2, `*`))
      pe_r<-pm_r^2+(1-pm_r)^2
      ki_r<-(po_r-pe_r)/(1-pe_r)
      
      sumnr1<-rowSums(sweep(scenario, MARGIN=2,r, `*`))
      sumnr2<-rowSums(sweep(scenario, MARGIN=2,r^2, `*`))
      sumnr3<-rowSums(sweep(scenario, MARGIN=2,r^3, `*`))
      sumnr4<-rowSums(sweep(scenario, MARGIN=2,r^4, `*`))
      
      poi_2<-(4*sumnr4-8*nrat*sumnr3+4*nrat*(2*nrat-1)*sumnr2-4*nrat^2*(nrat-1)*sumnr1+nrat^2*(nrat-1)^2)/(nrat^2*(nrat-1)^2)
      pei_2<-((2*pm_r-1)^2*sumnr2+2*nrat*(1-pm_r)*(2*pm_r-1)*sumnr1+nrat^2*(1-pm_r)^2)/(nrat^2)
      popei<-(2*(2*pm_r-1)*sumnr3+2*nrat*(2-3*pm_r)*sumnr2+nrat*(4*nrat*pm_r-3*nrat-2*pm_r+1)*sumnr1+nrat^2*(nrat-1)*(1-pm_r))/(nrat^2*(nrat-1))
      
      term2<-(po_r*pe_r-2*pe_r+po_r)^2
      term11<-(1-pe_r)^2*poi_2
      term12<--4*(1-pe_r)*(1-po_r)*popei
      term13<-4*(1-po_r)^2*pei_2
      
      var_ki<-(term11+term12+term13-term2)/((1-pe_r)^4)
      
      N<-ceiling(4*z^2*var_ki/w^2)
      given<-matrix(rep(c(a.level,w,nrat),each=nrow(scenario)),ncol=3)
      result<-data.frame(unique(cbind(given,matrix(pm_r,ncol=1),matrix(round(ki_r,3),ncol=1),matrix(N,ncol=1),scenario)))
      colnames(result) <- c("alpha level", "width", "nb raters","pm", "ki","min N",colnames(scenario))
      #rownames(result) <- rep("",nrow(result))
      return(result)
      #xranges(E, F, G = G, H = H)
    }
  }
}



size.pcam<-function(a.level=0.05,w=0.20,pm=0.5,pca=0.8,nrat=2,e_pca=0,e_pm=0,nb_s=40000,quant=100){
  
  min_p11<-pm-min(pm,1-pm)
  max_p11<-pm
  
  pa<-2*pca/(1+pca)
  
  min_pa<-2*min_p11/(pm+pm)
  max_pa<-2*max_p11/(pm+pm)
  
  min_pca<-min_pa/(2-min_pa)
  max_pca<-max_pa/(2-max_pa)
  
  if (a.level > 1 | a.level < 0) 
    stop("the significance level should be between 0 and 1")
  if (pca> max_pca |pca<min_pca) 
  {writeLines(paste("Stop: pca should be less than ",round(max_pca,2)));stop(paste("Stop: pca should be more than ",round(min_pca,2)))}
  
  z<-qnorm(1-a.level/2)
  
  if (nrat==2){
    
    
    N<-4*z^2*4*pa*(1-pa)*(2-pa)/(w^2*(pm+pm)*(2-pa)^4)
    return(ceiling(N))
    
  }
  
  if (nrat>2){
    
    
      npair<-nrat*(nrat-1)/2
      r<-seq(0,nrat,1)
      
      agree1_2<-(r*(r-1))/2             #number of pairs agreeing on cat1
      agree1_3<-(r*(r-1)*(r-2))/6       #number of triplets agreeing on cat1
      agree1_4<-(r*(r-1)*(r-2)*(r-3))/24#number of quartets agreeing on cat1
      
      
      E<-matrix(ncol=nrat+1,byrow=TRUE,data=c(agree1_2,r,rep(1,nrat+1)))
      F<-c(pm*pa*npair,pm*nrat,1)
      G<-diag(nrat+1)
      H<-rep(0,nrat+1)
      
      library(limSolve)
      scenario<- round(xsample(E = E, F = F, G = G, H = H, type = "mirror",iter=nb_s, burninlength = 1000)$X,2)
      
      
      A<-2*rowSums(sweep(scenario, MARGIN=2, agree1_2, `*`))
      B<-(nrat-1)*rowSums(sweep(scenario, MARGIN=2, r, `*`))
      
      pa_expected<-A/B
      
      pca_expected<-pa_expected/(2-pa_expected)
      
      valid<-scenario[pca_expected<=pca+e_pca & pca_expected>=pca-e_pca & B/(nrat*(nrat-1))<=pm+e_pm & B/(nrat*(nrat-1))>=pm-e_pm,]
      
      if (nrow(valid)<1) 
        stop("There is no valid scenario, please provide some positive value for tolerance")
      
      B_r<-(nrat-1)*rowSums(sweep(valid, MARGIN=2, r, `*`))
      A_r<-2*rowSums(sweep(valid, MARGIN=2, agree1_2, `*`))
      
      pa_r<-A_r/B_r
      
      pca_r<-pa_r/(2-pa_r)
      
      p3_r<-rowSums(sweep(valid, MARGIN=2, agree1_3, `*`))
      p4_r<-rowSums(sweep(valid, MARGIN=2, agree1_4, `*`))
      
      var_pa<-((2-(nrat-1)*pa_r)*(pa_r*(1-(nrat-1)*pa_r)+12*p3_r/B_r)+24*p4_r/B_r)/B_r
      
      var_pca<-var_pa*4/(2-pa_r)^4
      
      N<-ceiling(4*z^2*var_pca/w^2)
      
      # given<-matrix(rep(c(a.level,w,nrat),each=nrow(valid)),ncol=3)
      #result<-data.frame(unique(cbind(given,matrix(B_r/((nrat-1)*4),ncol=1),matrix(pa_r,ncol=1),matrix(N,ncol=1),valid)))
      #colnames(result) <- c("alpha level", "width", "nb raters","pm", "pa","min N",colnames(valid))
      #rownames(result) <- rep("",nrow(result))
      quantil<-quantile(N, probs = quant/100)
     
      return(ceiling(quantil))

  }
}


Prob_calc<-function(prop,kappa0,raters){
  Nprob<-raters+1
  Proba<-rep(NA,Nprob)
  lambda<-rep(NA,Nprob)
  sum_pi<-rep(NA,Nprob)
  
  mat_pi<-matrix(NA, nrow=Nprob,ncol=Nprob)
  mat_P<-matrix(NA, nrow=Nprob,ncol=Nprob)
  
  for (j in 1:Nprob){
    for (i in 1:Nprob){
      mat_pi[i,j]<-prop^(j-i)
    }
    sum_pi[j]<-sum(mat_pi[1:(j-1),j])
  }
  
  lambda[1]<-1
  lambda[2]<-prop
  for (k in 2:(Nprob-1)){
    lambda[k+1]<-prop^(k)+kappa0*(1-prop)*sum_pi[k]
  }
  
  for (j in 0:(Nprob-1)){
    for (i in 0:(Nprob-1)){
      mat_pi[i+1,j+1]<-(-1)^i*choose((raters-j), i)*lambda[i+j+1]
    }
  }
  
  
  for (i in 0:(Nprob-1)){
    Proba[i+1]<-choose(raters,i)*sum(mat_pi[1:(raters-i+1),(i+1)])
  }
  return(Proba)
}

chi_sq<-function(kappa0,kappa1,N,prop,raters){
  sum((N*Prob_calc(prop,kappa0,raters)-N*Prob_calc(prop,kappa1,raters))^2/(N*Prob_calc(prop,kappa1,raters)))
}


size.kappai<-function(kappa0,w,prop,raters,alpha=0.05){
  
  crit<-qchisq(1-alpha,1)
  
  kappaL<-kappa0-1.96*w/2
  kappaU<-kappa0+1.96*w/2
  
  N <- 10
  resultsl <- 0
  resultsu <- 0
  while ((abs(resultsl - 0.001) < crit) || (abs(resultsu -0.001) < crit)) {
    N <- N + 1
    resultsl <- chi_sq(kappa0,kappaL,N,prop,raters)
    resultsu <- chi_sq(kappa0,kappaU,N,prop,raters)
    if (is.infinite(resultsu)) {
      resultsu <- 0
    }
    if (is.infinite(resultsl)) {
      resultsl <- 0
    }
    
  }
  return(N)
}

# 
# grant<-read.csv("C:/Users/sophie.vanbelle/surfdrive/ONDERWIJS/OSLO/agreeshare/data/grant.csv",header=TRUE,na.string="NA");

many1.boot<-function(data,ind){
  
  many1(data[ind,])[1]
}

many2.boot<-function(data,ind){
  
  many2(data[ind,])[1]
}

po.boot<-function(data,ind){
  
  po(data[ind,])[1]
}


pa.boot<-function(data,ind){
  
  pa(data[ind,])[1]
}

pca.boot<-function(data,ind){
  
  pca(data[ind,])[1]
}


ki.boot<-function(data,ind){
  
  many1(data[ind,])[1]
}

kc.boot<-function(data,ind){
  
  many2(data[ind,])[1]
}
# library(boot)
# set.seed(123)
# boot.sample<-boot(grant, pa.boot,5000)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# plot(density(boot.sample$t))
# 
# a<-rep(NA,5000)
# for (i in 1:5000){
# a[i]<-mean(boot.sample$t[1:i])
# }
# 
# plot(seq(1:5000),a)
# 
# 

kc.bayes<-function(data,quantiles){
  
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  
  number<-1/(2^ncol(data))
  
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
  
  ind<-sapply(1:nrat,function(x) index[data.f[,x]==1])
  
  N2 <-30000
  
  alpha<-fre*n+number
  
  w<- lapply(alpha, function(x) rgamma(N2,shape=x,rate=1))
  
  t<-Reduce("+", w) 
  
  z<-sapply(w, function(x) x/t)
  
  if (npair>1){
    
    p11<-rowSums(z[,c(ind_1)])/(npair)
    p22<-rowSums(z[,c(ind_2)])/(npair)
    p12<-rowSums(z[,c(ind_3)])/(npair)
    p21<-rowSums(z[,c(ind_4)])/(npair)
    
    po<-p11+p22
    
    
    marge_<-sapply(1:nrat,function(x) rowSums(z[,c(ind[,x])]))  
    
    pe_pair<-matrix(NA,nrow=N2,ncol=npair)
    
    for (b in 1:npair){
      pe_pair[,b]<-marge_[,rat.pair[b,1]]*marge_[,rat.pair[b,2]]+(1- marge_[,rat.pair[b,1]])*(1-marge_[,rat.pair[b,2]])
    }
    
    pe<-rowMeans(pe_pair)
    
    kappa<-(po-pe)/(1-pe)
    
  }
  
  if (npair==1){
    p11<-z[,c(ind_1)]/(npair)
    p22<-z[,c(ind_2)]/(npair)
    p12<-z[,c(ind_3)]/(npair)
    p21<-z[,c(ind_4)]/(npair)
    
    po<-p11+p22
    
    marge_<-sapply(1:nrat,function(x) rowSums(z[,c(ind[,x])]))  
    
    pe_pair<-matrix(NA,nrow=N2,ncol=npair)
    
    for (b in 1:npair){
      pe_pair[,b]<-marge_[,rat.pair[b,1]]*marge_[,rat.pair[b,2]]+(1- marge_[,rat.pair[b,1]])*(1-marge_[,rat.pair[b,2]])
    }
    pe<-rowMeans(pe_pair)
    
    kappa<-(po-pe)/(1-pe)
  }
  
  mean.kc<-mean(kappa)
  
  var.kc<-var(kappa)
  
  q.kc<-quantile(kappa,quantiles)
  
  result<-matrix(c(mean.kc,sqrt(var.kc),q.kc),nrow=1)
  colnames(result) <- c("Kc", "SE(Kc)",rep(NA,length(q.kc)))
  rownames(result) <- ""
  return(result)
  
}


ki.bayes<-function(data,quantiles){
  
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  
  number<-1/(2^ncol(data))
  
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
  
  ind<-sapply(1:nrat,function(x) index[data.f[,x]==1])
  
  N2 <-30000
  
  alpha<-fre*n+number
  
  w<- lapply(alpha, function(x) rgamma(N2,shape=x,rate=1))
  
  t<-Reduce("+", w) 
  
  z<-sapply(w, function(x) x/t)
  
  if (npair>1){
    
    p11<-rowSums(z[,c(ind_1)])/(npair)
    p22<-rowSums(z[,c(ind_2)])/(npair)
    p12<-rowSums(z[,c(ind_3)])/(npair)
    p21<-rowSums(z[,c(ind_4)])/(npair)
    
    po<-p11+p22
    
    
    marge_<-sapply(1:nrat,function(x) rowSums(z[,c(ind[,x])]))  
    
    pe_pair<-matrix(NA,nrow=N2,ncol=npair)
    
    for (b in 1:npair){
      pe_pair[,b]<-((marge_[,rat.pair[b,1]]+marge_[,rat.pair[b,2]])/2)^2+((1-marge_[,rat.pair[b,1]]+1-marge_[,rat.pair[b,2]])/2)^2
    }
    
    pe<-rowMeans(pe_pair)
    
    kappa<-(po-pe)/(1-pe)
    
  }
  
  if (npair==1){
    p11<-z[,c(ind_1)]/(npair)
    p22<-z[,c(ind_2)]/(npair)
    p12<-z[,c(ind_3)]/(npair)
    p21<-z[,c(ind_4)]/(npair)
    
    po<-p11+p22
    
    marge_<-sapply(1:nrat,function(x) rowSums(z[,c(ind[,x])]))  
    
    pe_pair<-matrix(NA,nrow=N2,ncol=npair)
    
    for (b in 1:npair){
      pe_pair[,b]<-((marge_[,rat.pair[b,1]]+marge_[,rat.pair[b,2]])/2)^2+((1-marge_[,rat.pair[b,1]]+1-marge_[,rat.pair[b,2]])/2)^2
    }
    pe<-rowMeans(pe_pair)
    
    kappa<-(po-pe)/(1-pe)
  }
  
  mean.ki<-mean(kappa)
  
  var.ki<-var(kappa)
  
  q.ki<-quantile(kappa,quantiles)
  
  result<-matrix(c(mean.ki,sqrt(var.ki),q.ki),nrow=1)
  colnames(result) <- c("Ki", "SE(Ki)",rep(NA,length(q.ki)))
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
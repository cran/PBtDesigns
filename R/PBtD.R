###################################################################################################################################
#' Partially Balanced t-Designs (PBtDesigns)
#'
#' @param v Number of treatments
#' @param Series Series of Partially Balanced t-Designs
#'
#' @return Three series are given for generating of partially balanced t-designs namely Series 1, Series 2 and Series 3.
#'
#' Series 1 are designs having equal block sizes and with treatment structure 4(t + 1).
#'
#' Series 2 are designs having equal block sizes and with treatment as a prime number.
#'
#' Series 3 consists of designs with unequal block sizes and with treatment structure n(n-1)/2.
#'
#' This function generates partially balanced t-designs along with their parameters, information matrices, average variance factors and canonical efficiency factors.
#' @export
#'@importFrom utils combn
#' @examples
#' library(PBtDesigns)
#' PBtD(7,2)
#' @description This package contains functions named PBtD() for generating partially balanced t-designs along with their parameters, information matrices, average variance factors and canonical efficiency factors.
#'@references
#'1) Karmakar, S., Varghese, C., Jaggi, S. & Harun, M. (2021)< DOI:10.1080/03610918.2021.2008436>. " Partially Balanced t-designs ".
#'
#'2) Raghavarao, D. & Zhou, B. (1998)<https://doi.org/10.1080/03610929808832657> " Universal optimality of UE 3-designs for a competing effects model ".
#'
#'3) Karmakar, S., Varghese, C., Jaggi, S. & Harun, M. (2022)." Partially Balanced t-designs with unequal block sizes ".

###################################################################################################################################
PBtD<-function(v,Series){
  if(Series==1){
    v=v-1
    ###from BIBD
    t=(v-3)/4
    v=4*t+3
    if(v==3){
      pr=2
    }
    if(v==5){
      pr=2
    }
    if(v==7){
      pr=3
    }
    if(v==11){
      pr=2
    }
    if(v==13){
      pr=2
    }
    if(v==17){
      pr=3
    }
    if(v==19){
      pr=2
    }
    if(v==23){
      pr=5
    }
    if(v==29){
      pr=3
    }
    if(v==31){
      pr=3
    }

    seq=seq(2,(v-1),by=2)
    ele=c()
    for(i in seq){
      x=(pr)^i
      ele=c(ele,x)
    }
    ele=ele%%v
    ele=sort(ele)
    mat=c()
    for(i in ele){
      x=seq(i,i+v-1,1)
      mat=cbind(mat,x)
    }
    mat=mat%%v
    colnames(mat)<-NULL
    mat[mat==0]=v
    mat11=mat
    mat1=cbind(mat,(max(mat)+1))

    ############second step
    all=c(1:v)
    mat2=matrix(0,nrow=nrow(mat1),ncol=ncol(mat1))
    for(i in 1:nrow(mat2)){
      mat2[i,]<-setdiff(all,mat11[i,])
    }
    mat2
    final=rbind(mat1,mat2)
    message("t-design")
    rownames=c()
    for(i in 1:nrow(final)){
      x=paste0("Block ", i)
      x<-noquote(x)
      rownames<- rbind(rownames,x)
    }
    rownames(final)<-c(rownames)
    print("                   Treatmets ",quote=F)
    print(final)
  }

  #####################################################################################################################################

  ##########From mols
  if(Series==2){
    ini=c(1,2,3)
    add=c(0,1,2)
    times=(v-1)/2
    inirow=c(ini)
    for(i in 2:times){
      x=ini+add
      inirow=c(inirow,x)
      add=add+c(0,1,2)
    }
    allini=inirow
    allini1=allini
    for(i in 1:(v-1)){
      x=allini1+i
      allini=rbind(allini,x)
    }
    final=allini%%v
    final[final==0]<-v
    finalmat=c()
    for(i in 1:times){
      finalmat=rbind(finalmat,final[,((3*i)-2):((3*i))])
    }
    #######################row names
    final=finalmat
    message("t-design")
    rownames=c()
    for(i in 1:nrow(final)){
      x=paste0("Block ", i)
      x<-noquote(x)
      rownames<- rbind(rownames,x)
    }
    rownames(final)<-c(rownames)
    print("       Treatmets ",quote=F)

    print(final)
  }
  ##################################################################################################################
  ##From Triangular association
  if(Series==3){
    m=2
    n1=(1+sqrt(1-4*1*(-2*v)))/2
    n2=(1-sqrt(1-4*1*(-2*v)))/2
    if(n1>0){
      n=n1
    }else{
      n=n2
    }
    #####################

    ##################Start from here
    matt=matrix(0,nrow=n,ncol=n)
    ss<-c(1:choose(n,2))

    for(i in 1:(n-1)){
      matt[i,(i+1):n]<-ss[1:(n-i)]
      ss=setdiff(ss,ss[1:(n-i)])
    }
    for(i in 1:ncol(matt)){
      matt[,i]<-matt[i,]
    }
    matt
    #############
    ss1<-c(1:choose(n,2))
    com<-c()
    for(i in 1:(n-2)){
      x=c(combn(ss1[1:(n-i)],2))
      ss1=setdiff(ss1,ss1[1:(n-i)])
      com<-c(com,x)
    }
    com=matrix(com,nrow=length(com)/2,byrow=T)
    ############
    intsec=c()
    for(i in 1:nrow(com)){
      a=com[i,1]
      b=com[i,2]
      pa=which(t(matt)==a)%%n
      pa[pa==0]<-n
      pa=pa[1]
      pb=which((matt)==b)%%n
      pb[pb==0]<-n
      pb=pb[1]
      intsec=c(intsec,matt[pa,pb])
    }
    first=cbind(com,intsec)
    ####################################
    remaining=c()
    for(i in 1:nrow(com)){
      a=com[i,1]
      b=com[i,2]
      pos.a=which(matt==a,arr.ind = T)
      pos.b=which(matt==b,arr.ind = T)
      allrows=unique(c(pos.a[,1],pos.b[,1]))
      remaining=rbind(remaining,setdiff(c(1:v),c(matt[c(allrows),])))
    }
    ####################
    second=cbind(first,remaining)
    ##############
    matt[matt==0]=NA
    nmatt=matrix(NA,nrow=n,ncol=n-1)
    for(i in 1:n){
      x=sort(matt[i,])
      nmatt[i,]<-x
    }
    if(v==10){
      design=rbind(second,nmatt)
      colnames(design)=NULL
      message("t-design")
      print("                   Treatmets ",quote=F)
      rownames=c()
      for(i in 1:nrow(design)){
        x=paste0("Block ", i)
        x<-noquote(x)
        rownames<- rbind(rownames,x)
      }
      rownames(design)<-c(rownames)
      print(design)
    }else{
      na.mat=matrix(NA,nrow=n,ncol=ncol(second)-n+1)
      nmatt<-cbind(nmatt,na.mat)
      ##############
      design=rbind(second,nmatt)
      ###################row and column name
      rownames=c()
      for(i in 1:nrow(design)){
        x=paste0("Block ", i)
        x<-noquote(x)
        rownames<- rbind(rownames,x)
      }
      rownames(design)<-c(rownames)

      ##################
      colnames(design)=NULL
      message("t-design")
      print("                   Treatmets ",quote=F)
      print.table(design,na.print="")
    }
    final=design
    final[is.na(final)]=0
  }
  ##########

  ###############################canonical efficiency factor#########################################################################################

  ############PARAMETERS
  message("Parameters Of the t-Design")

  if(Series==3){
    if(v!=10){
      cat(c("\n","Number Of treatments (v) =",max(final),"\n","Number Of blocks (b) =",nrow(final),"\n","Number Of replications (r) =",length(which(final==max(final))),"\n",
            "Block Size (k1) =",ncol(final),"\n", "Block Size (k2) =",n-1))
    }else{
      cat(c("\n","Number Of treatments (v) =",max(final),"\n","Number Of blocks (b) =",nrow(final),"\n","Number Of replication (r) =",length(which(final==max(final))),"\n",
            "Block size (k) =",ncol(final)))
    }
  }
  if(Series==1){
    cat(c("\n","Number Of treatments (v) =",max(final),"\n","Number Of blocks (b) =",nrow(final),"\n","Number Of replications (r) =",length(which(final==max(final))),"\n",
          "Block size (k) =",ncol(final)),"\n","Pairwise appearance of first associates(lambda_1)=",2*t,"\n","Pairwise appearance of second associates(lambda_2)=",2*t+1,"\n", "Number of blocks in which triplets of same group are appearing (delta_1)=", t-1, "\n","Number of blocks in which triplets of same group are appearing (delta_2)=", t)
  }
  if(Series==2){
    cat(c("\n","Number Of treatments (v) =",max(final),"\n","Number Of blocks (b) =",nrow(final),"\n","Number Of replications (r) =",length(which(final==max(final))),"\n",
          "Block Size (k) =",ncol(final)),"\n","Pairwise occurance of treatments (lambda_)=",3,"\n", "Number of blocks in which triplets of same group are appearing (delta_1)=", 1, "\n","Number of blocks in which triplets of same group are appearing (delta_2)=", 0)
  }

  ##############R(replication) matrix
  rep_r<-length(final[which(final==final[1,1])])
  R_matrix<-diag(rep_r,max(final))
  #################

  ################(block size) matrix
  block_size<-ncol(final)
  K_matrix<-diag(block_size,nrow(final))
  if(Series==3 && v!=10){
    kk=nrow(final)-n
    for(i in (kk+1):(kk+n)){
      K_matrix[i,i]<-(n-1)
    }
  }
  #############incidence matrix(N-v*b)
  t.incidence<-function(bibd){
    bibd
    b=nrow(bibd)
    kk=ncol(bibd)
    v<-max(bibd)
    cc<-c(1:v)

    incident<-matrix(0,nrow=v, ncol=b)
    vv=1
    while(vv<=v){
      ###########raw position of a element
      x<-which(bibd %in% c(cc[vv]))

      #########################position identify only
      k=1
      while(k<=length(x)){
        if(x[k]%%b!=0){
          x[k]<-x[k]%%b
        }else{
          x[k]<-b
        }

        k=k+1
      }
      ################

      ss=1
      while(ss<=length(x)){
        incident[vv,x[ss]]<-1
        ss=ss+1
      }
      vv=vv+1
    }
    N_prime<-t(incident)
  }
  N_matrix<-t(t.incidence(final))
  #########
  cat("\n")
  ####################C matrix
  C_matrix<-R_matrix-(N_matrix%*%MASS::ginv(K_matrix)%*%t(N_matrix))
  cat("\n")
  message("C Matrix")
  cat("\n")
  print(C_matrix)
  cat("\n")
  ############
  ########P matrix generation
  v=max(final)
  p_matrix<-matrix(0,nrow=choose(v,2),ncol=max(final))
  #########
  elepos<-t(combn(v,2))
  #####
  for(i in 1:nrow(p_matrix)){
    p_matrix[i,(elepos[i,])]<-c(1,-1)
  }
  ########## Variance covariance part
  variances<-(p_matrix)%*%MASS::ginv(C_matrix)
  variances<-variances%*%t(p_matrix)
  #######variance part
  var<-diag(variances)
  ########## avg variance
  Avg_var<-mean(var)
  message("Average Variance Factor")
  cat("\n")
  print(Avg_var)
  cat("\n")
  #####cef
  CEF<-(rep_r)/2*Avg_var
  #######canonical efficiency factor
  eigen<-c(eigen(C_matrix)$values)

  ##############
  p<-c()
  for(i in eigen){
    if(i>(10^(-6))){
      p<-c(p,i)
    }
  }
  #########harmonic mean
  harmonic_mean<-1/(mean(1/p))
  ####canonical efficiency factor (need to change replication)
  cannonical_efficiency<-(1/(rep_r))*harmonic_mean
  ############
  message("Cannonical Efficiency Factor")
  cat("\n")
  print(cannonical_efficiency)
}













#MENYESUAIKAN DENGAN CODING MB NURUL
dataa=read.delim("clipboard",header=TRUE)
#matrik sigular:Matriks singular merupakan matriks yang nilai determinannya nol, sedangkan nilai determinan matriks non singular tidak sama dengan nol.
#https://en-m-wikipedia-org.translate.goog/wiki/Singular_value_decomposition?_x_tr_sl=auto&_x_tr_tl=en&_x_tr_hl=id&_x_tr_pto=wapp 

ginverse = function(x,eps=1e-016)
{
  x<-as.matrix(x) #mendefeni bahwa x adalah matriks
  xsvd<-svd(x) #bahwa x merupakan Singular Value Decomposition (pemfaktoran matriks dengan mengurai suatu matriks ke dalam tiga matriks unit yaitu U, V, dan sebuah matriks diagonal S yang berisi faktor skala yang disebut dengan nilai singular.)
  diago<-xsvd$d[xsvd$d>eps] #membuat matriks diagonal
  if(length(diago)==1) 
  {
    xplus<-as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u[,1])/diago) #
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)]%*%diag(1/diago)%*%t(xsvd$u[,1:length(diago)])#R menyediakan sebuah fungsi yaitu diag() untuk mengakses nilai-nilai pada diagonal utama sebuah matriks.    }
    return(xplus)
  }
}
# Nilai increment: https://www.duniailkom.com/tutorial-belajar-c-jenis-jenis-operator-increment-dan-decrement-bahasa-c/ 
#Titik knot dari nilai quantil jadi kita 
#P:polinomial atau suku banyak (juga ditulis sukubanyak) adalah pernyataan matematika yang melibatkan jumlahan perkalian pangkat dalam satu atau lebih variabel dengan koefisien.
quant<-function(pred,P)#P:Banyaknya Knot + 1 (k+1)
{
  r<-quantile(pred,seq(0,1,by=1/P))#seq(from = 1, to = 10) # 1 sampai 10 dengan increment(kenaikan) 1 (default by = 1)
  return(r)
}

#truncated: calon titik not yang baru
trun<-function(pred,k,m)
{
  pred[pred<k]<-k #jika nilai var kurang dari titik not maka titik Knot = 0, karena k-k
  b<-(pred-k)^m #jika lebih atau sama titik knot maka nilainya tetap (x-e)^m (m:orde)
  return(b)#return untuk menyimpan data yang telah ada
}

matrikx<-function(pred1,pred2,m,k)
{
  predbaru1<-(unique(pred1))#unique untuk mencari titik knotnya
  predbaru2<-(unique(pred2))
  n<-length(pred1)#panjang variabel prediktor, lenght hanya 1 karena panjang data sama 
  w1<-quant(predbaru1,k+1)#Calon titik knot yang akan digunakan karena quant tadi titik knot
  w2<-quant(predbaru2,k+1) #quant<-function(pred,P)
  z1<-matrix(0,n,m+1) #matrix 
  for(i in 1:(m+1))
  {
    z1[,i]<-pred1^(i-1)#berakitana dengan banyaknya orde
  }
  z2<-matrix(0,n,k)
  for(j in 1:k)
  {
    z2[,j]<-trun(pred1,w1[j+1],m)
  }
  z3<-matrix(0,n,m)
  for(i in 1:m)
  {
    z3[,i]<-pred2^((i+1)-1) #(i+4)
  }
  z4<-matrix(0,n,k)
  for(j in 1:k)
  {
    z4[,j]<-trun(pred2,w2[j+1],m)
  }
  x<-cbind(z1,z2,z3,z4)
  return(x)
}
#matrikx(dataa$x1,dataa$x2,1,1)
matrikd<-function(m,k)
{
  d1<-matrix(0,      2*m+1,      2*m+2*k+1) 
  d2<-matrix(0,      2*k,        2*m+1)
  d3<-diag(2*k)
  d4<-cbind(d2,d3)
  d<-rbind(d1,d4)
  return(d)
}
#matrikd(1,1)

beta<-function(respon,pred1,pred2,m,k,lam)
{ 
  n<-length(respon)
  y<-as.vector(respon)
  x<-matrikx(pred1,pred2,m,k)
  d<-matrikd(m,k)
  b1<-(t(x)%*%x)+(n*lam*d)
  b2<-ginverse(b1)
  beta<-b2%*%t(x)%*%y
  print(beta)
  return(beta)
}


#beta(dataa$y,dataa$x1,dataa$x2,1,1,0.1)

ytopi<-function(respon,pred1,pred2,m,k,lam)
{
  n<-length(respon)
  y<-as.vector(respon)
  x<-matrikx(pred1,pred2,m,k)
  d<-matrikd(m,k)
  b1<-(t(x)%*%x)+(n*lam*d)
  b2<-ginverse(b1)
  beta<-b2%*%t(x)%*%y
  Hlam<-x%*%b2%*%t(x)
  ytopi<-x%*%beta
  return(ytopi)
}
#ytopi(dataa$y,dataa$x1,dataa$x2,1,1,0.1)
Hlam<-function(respon,pred1,pred2,m,k,lambda) 
{
  n<-length(respon)
  y<-as.vector(respon)
  h<-lambda
  x<-matrikx(pred1,pred2,m,k)
  d<-matrikd(m,k)
  f1<-(t(x)%*%x)+(n*h*d)
  f2<-ginverse(f1)
  beta<-f2%*%t(x)%*%y
  Hlambda<-x%*%f2%*%t(x)
  return(Hlambda)
}
#Hlam(dataa$y,dataa$x1,dataa$x2,1,1,0.1)

gcv<-function(respon,pred1,pred2,m,k,lam)
{
  n<-length(respon)
  y<-as.vector(respon)
  x<-matrikx(pred1,pred2,m,k)
  d<-matrikd(m,k)
  b1<-(t(x)%*%x)+(n*lam*d)
  b2<-ginverse(b1)
  beta<-b2%*%t(x)%*%y
  Hlam<-x%*%b2%*%t(x)
  ytopi<-x%*%beta
  MSE<-(t(y-ytopi)%*%(y-ytopi))/n
  GCV<-MSE/(1-((1/n)*sum(diag(Hlam))))^2
  nilai<-cbind(GCV,MSE)
  return(nilai)
}
#gcv(dataa$y,dataa$x1,dataa$x2,1,1,0.1)

#Untuk mencari Lambda optimal
carilambda<-function(respon,pred1,pred2,m,k,bb,ba,ic)
{
  y<-as.vector(respon)
  n<-length(respon)
  x<-matrikx(pred1,pred2,m,k)
  d<-matrikd(m,k)
  w1<-quant(pred1,k+1)
  w2<-quant(pred2,k+1)
  z1<-matrix(0,n,m+1)
  for(a in 1:k)
  {
    cat("titik knot prediktor 1 ke-[",a,"]=",w1[a+1],"\n")
    cat("titik knot prediktor 2 ke-[",a,"]=",w2[a+1],"\n")
  }
  lambda<-seq(bb,ba,ic)
  nk<-length(lambda)
  GCV<-rep(0,nk) #membuat vektor dengan pengulangan 0 sampai nk kali
  MSE<-rep(0,nk)
  for(i in 1:nk)
  {
    b1<-(t(x)%*%x)+(n*lambda[i]*d)
    b2<-ginverse(b1)
    betatopi<-b2%*%t(x)%*%y
    aps<-x%*%b2%*%t(x)
    ytopi<-x%*%betatopi
    MSE<-(t(y-ytopi)%*%(y-ytopi))/n
    GCV[i]<-MSE/(1-((1/n)*sum(diag(aps))))^2
  }
  s<-matrix(c(lambda,GCV),length(lambda),2)
  GCVmin<-min(GCV)
  lambdaopt<-s[s[,2]==min(GCV),1]
  plot(lambda,GCV,type="l")
  c<-cbind(lambdaopt,GCVmin)
  return(c)
}
carilambda(dataa$y,dataa$x1,dataa$x2,1,1,0.1,100,0.1)
#ADA ERROR DI W2 SCRIPT DILUAR BATAS
carioptimal<-function(respon,pred1,pred2)
{
  bb<-as.numeric(readline("Masukkan batas bawah lambda:")) #batas bawah
  ba<-as.numeric(readline("Masukkan batas atas lambda:"))  #batas atas
  ic<-as.numeric(readline("Masukkan nilai increment lambda:"))   #
  n<-length(respon)   #banyaknya data
  y<-as.vector(respon) #variabel respon
  x1<-as.vector(pred1)   #variabel prediksi
  x2<-as.vector(pred2)   #variabel prediksi
  i<-1   #jumlah iterasi
  repeat
  {
    k<-2   #jumlah knot
    repeat
    {
      hasil1<-carilambda(y,x1,x2,i,1,bb,ba,ic)   #hasil1 adalah
      lambda1<-hasil1[,1]
      gcv1<-hasil1[,2]
      m1<-matrix(0,1,4)#kemungkinan nanti ada 8 jika variabel x nya ada 4, karena ada matrik trun dan matriks knots
      m1[1,]<-c(i,1,lambda1,gcv1)
      m2<-matrix(0,(k-1),4)
      for(n in 1:(k-1))
      {
        hasil2<-carilambda(y,x1,x2,i,n+1,bb,ba,ic)
        lambda<-hasil2[,1]
        gcv<-hasil2[,2] 
        m2[n,]<-c(i,n+1,lambda,gcv)
      }
      m3<-rbind(m1,m2)
      m4<-matrix(0,1,4)
      if(m3[k,4]<m3[k-1,4])
      {
        m4[1,]<-m3[k,]
      }
      else
      {
        m4[1,]<-m3[k-1,]
      }
      if(m3[k,4]>(m3[k-1,4]))break
      k<-k+1
    }
    if(i==1)
    {
      mgcv<-m4[,4]
      mgcvlama<-mgcv
    }
    else
    {
      mgcv<-m4[,4]
      if(mgcv>mgcvlama)break
      mgcvlama<-mgcv
    }
    i<-i+1
    mopt<-m4
  }
  colnames(mopt)=c("Orde","Jumlah Knot","Lambda","GCV")
  return(mopt)
}
carioptimal(dataa$y,dataa$x1,dataa$x2)


#Estimasi Parameter PSpline
pspline.dua<-function(respon,pred1,pred2)
{
  y<-as.vector(respon)
  x1<-as.vector(pred1)
  x2<-as.vector(pred2)
  n<-length(y)
  optimal<-carioptimal(y,x1,x2)
  m<-optimal[,1]
  k<-optimal[,2]
  h<-optimal[,3]
  predictorbaru1<-unique(x1)
  predictorbaru2<-unique(x2)
  w1<-quant(predictorbaru1,k+1)
  w2<-quant(predictorbaru2,k+1)
  cat("orde:",m,"\n")
  cat("lambda:",h,"\n")
  cat("jumlah knot:",k,"\n")
  #cat("quantile(",1/(k+1),")=",w1,"\n")
  #cat("quantile(",1/(k+1),")=",w2,"\n")
  for(i in 1:k)
  {
    cat("titik knots prediktor 1[",i,"]=",w1[i+1],"\n")
  }
  cat("quantile(",1/(k+1),")=",w1,"\n")
  #cat diatas diulang sampe banyak variabel x jadi sampe 4 kali
  for(i in 1:k)
  {
    cat("titik knots prediktor 2[",i,"]=",w2[i+1],"\n")
  }
  cat("quantile(",1/(k+1),")=",w2,"\n")
  hasil<-gcv(y,x1,x2,m,k,h)
  GCV<-hasil[,1]
  MSE<-hasil[,2]
  beta<-beta(y,x1,x2,m,k,h)
  fstar<-ytopi(y,x1,x2,m,k,h)
  error<-y-fstar
  cat("\nGCV=",(GCV))
  cat("\nMSE=",format(MSE),"\n")
  cat("\n-------------------------------------------------------------")
  cat("\n   x1  x2  Y      (Y*topi)       error")
  cat("\n-------------------------------------------------------------")
  for(i in 1:n)
    cat("\n",x1[i],"     ",x2[i],"    ",y[i],"     ",fstar[i],"    ",error[i])
  cat("\n-------------------------------------------------------------")
  for(i in 1:(2*m+2*k+1))
    cat("\n nilai beta[",i,"]=",beta[i])
  Xurut<-seq(1, length(y), by=1)
  win.graph()
  plot(Xurut,y,xlim=c(min(Xurut),max(Xurut)),ylim=c(min(y),max(y)),xlab="predictor",ylab="respon") 
  title("Fungsi Penalize untuk dua Prediktor")
  par(new=T)
  plot(Xurut,fstar,type="l",xlim=c(min(Xurut),max(Xurut)),ylim=c(min(y),max(y)),xlab="prediktor",ylab="respon")
}
pspline.dua(dataa$y,dataa$x1,dataa$x2)
#carioptimal(dataa$Y,dataa$X1,dataa$X2, dataa$X3)
#carioptimal(dataa$Y,dataa$X2)
#carioptimal(dataa$Y,dataa$X3)
#carioptimal(dataa$Y,dataa$X4)
#pspline.empat(dataa$Y,dataa$X1,dataa$X2,dataa$X3,dataa$X4)



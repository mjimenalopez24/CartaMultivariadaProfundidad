#### Calcular la matriz de los datos ####

library(mvtnorm)


## Fase I
datos1<-function(n,p,pi,independientes=c("SI","NO"),cambio){
  mu<-c()
  mu<-c(rep(0,p))
  Sigma<-matrix(nrow=p,ncol=p)
  if(independientes=="SI"){ 
    for(i in 1:p){
      for(j in 1:p){
        if(i==j){
          Sigma[i,j]<-1
        }
        else{
          Sigma[i,j]<-0
        }
      }
    }
  }
  if(independientes=="NO"){
    for(i in 1:p){
      for(j in 1:p){
        if(i==j){
          Sigma[i,j]<-1
        }
        else{
          Sigma[i,j]<-0.9
        }
      }
    }
  }
  data<-matrix(nrow=n,ncol=p)
  if(pi==0){
    data<-rmvnorm(n=n, mean=mu, sigma=Sigma)
  }
  if(pi!=0){
    mu1<-c(rep(cambio,p))
    delta<-t(mu1)%*%(solve(Sigma))%*%mu1
    data<-rbind(rmvnorm(n=n*(1-pi), mean=mu, sigma=Sigma),rmvnorm(n=n*pi, mean=mu1, sigma=Sigma))
  }
  return(list(datos=data,ncp=delta))
}

## Fase II
datos2<-function(n,p,independientes=c("SI","NO"),cambio){
  mu<-c()
  mu<-c(rep(cambio,p))
  Sigma<-matrix(nrow=p,ncol=p)
  if(independientes=="SI"){
    for(i in 1:p){
      for(j in 1:p){
        if(i==j){
          Sigma[i,j]<-1
        }
        else{
          Sigma[i,j]<-0
        }
      }
    }
  }
  if(independientes=="NO"){
    for(i in 1:p){
      for(j in 1:p){
        if(i==j){
          Sigma[i,j]<-1
        }
        else{
          Sigma[i,j]<-0.9
        }
      }
    }
  }
  data<-matrix(nrow=n,ncol=p)
  delta<-t(mu)%*%(solve(Sigma))%*%mu
  data<-rmvnorm(n=n, mean=mu, sigma=Sigma)
  return(list(datos=data,ncp=delta))
}

#### Calcular la media de los datos ####


# primero hay que encontrar mu_MM
# toca calcular la matriz jacobiana
# ver si es minimo
library(robustX)
library(pracma)
library(MASS)

iter<-1000
set.seed(1234)

p<-30
n<-50
T_SR<-matrix(NA,nrow=iter,ncol=n)
T_SR2<-c(rep(0,iter))
beta<-c(rep(0,iter))
changes<-c(0.2,0.5,1,1.2,1.5,2,2.5)
uno_beta<-rep(0,length(changes))
nocentralidad<-rep(0,length(changes))

t<-0
for(change in changes){
  t<-t+1
for(r in 1:iter){
  daticos<-datos1(50,p,0.2,"SI",change)
  X<-daticos$datos
  #nocentralidad<-daticos$ncp
mu_MM <- L1median(X)$estimate

e <- c(rep(1,ncol(X)))
nu <- (t(mu_MM)%*%e)/ncol(X)
nue<-nu%*%e
nue<-matrix(nue)

# tambien toca estimar \eta y para eso se necesita A
 # por si acaso A no es invertible
vector_etas<-c(rep(0,nrow(X)))
for(i in 1:nrow(X)){
  A <- matrix(NA,nrow=ncol(X),ncol=ncol(X))
  B <- matrix(NA,nrow=ncol(X),ncol=ncol(X))
  ip <- diag (1, nrow = ncol(X)) # matriz identidad p
  v <- t(X[i,])
  v <- matrix(v)
  A <-(1/norm(v,type="2"))*(ip-(v%*%t(v))/(norm(v,type="2"))^2)
  B <- (v%*% t(v))/(norm(v,type="2"))^2
  inversa<-ginv(A)
  matriz<-(1/nrow(X))*(inversa%*%B%*%inversa)
  vector_etas[i] <- (sum(diag(matriz)))/(norm(mu_MM-nu%*%e))^2
}
eta<-mean(vector_etas)
nu<-mean(nu)
mu_SH<-(1-eta)*mu_MM+(eta*nu)*e

# luego S_Sh
# pero esta matriz no es de tamaño pxp sino nxn pero pues es la estimacion de la matriz de varianzas y covarianza, no?


listamatrices<-list()
cant_matrices<-(ncol(X))*(ncol(X))
for(k in 1:cant_matrices){
  producto<-matrix(NA,nrow=ncol(X),ncol=ncol(X))
  for(i in 1:ncol(X)){
    for(j in 1:ncol(X)){
      u<-t(X[i,])
      v<-t(X[j,])
      u<-matrix(u)
      v<-matrix(v)
      producto<-(u-mu_SH)%*%t(v-mu_SH)
    }
  }
  listamatrices[[k]]<-producto
}

tamano_matriz <- dim(listamatrices[[1]])
mediana_matrices <- matrix(0, nrow = tamano_matriz[1], ncol = tamano_matriz[2])
for (i in 1:tamano_matriz[1]) {
  for (j in 1:tamano_matriz[2]) {
    valores_elemento <- sapply(listamatrices, function(matriz) matriz[i, j])
    mediana_matrices[i, j] <- median(valores_elemento)
  }
}
S_SH <- 2.198*mediana_matrices

Id<-diag(1,nrow=ncol(X),ncol=ncol(X))
etanu<-eta*nu
etanu<-mean(etanu)
Sigma_SH<-(1-eta)*S_SH+etanu*Id

# Media SR

# Distancia Mahalanobis
dist_mah<-c(rep(0,nrow(X)))
for(i in 1:nrow(X)){
  u<-X[i,]
  u<-matrix(u)
  dist_mah[i]<-(t(u-mu_SH))%*%Sigma_SH%*%(u-mu_SH)
}
  
w<-c(rep(0,nrow(X)))
cuantil<-qchisq(0.975,nrow(X)+1)
for(i in 1:nrow(X)){
  if(dist_mah[i]<cuantil){
    w[i]<-1
  }else{
    w[i]<-0
  }
}
  
suma1<-c(rep(0,ncol(X)))
suma2<-0
for(i in 1:nrow(X)){
  u<-t(X[i,])
  u<-matrix(u)
  suma1<-suma1+w[i]*u
  suma2<-suma2+w[i]
}
mu_SR<-(1/suma2)*suma1


suma1<-matrix(rep(0,(ncol(X))*(ncol(X))),nrow=ncol(X),ncol=ncol(X))
suma2<-0
for(i in 1:nrow(X)){
  u<-t(X[i,])
  u<-matrix(u)
  suma1<-suma1+w[i]*((u-mu_SR)%*%t(u-mu_SR))
  suma2<-suma2+w[i]
}
Sigma_SR<-(1/suma2)*suma1


for(i in 1:nrow(X)){
  u<-X[i,]
  u<-matrix(u)
  T_SR[r,i]<-(t(u-mu_SR))%*%Sigma_SR%*%(u-mu_SR)
}
daticos2<-datos2(1,p,"SI",change)
X2<-daticos2$datos
#nocentralidad2<-daticos2$ncp
X2<-t(X2)
T_SR2[r]<-(t(X2-mu_SR))%*%Sigma_SR%*%(X2-mu_SR)
p<-ncol(X)
n<-nrow(X)
if(p==2){
  a3<-9.210
  a1<-204
  a2<-0.850
}
if(p==6){
    a3<-16.812
    a1<-327500
    a2<-2.253
}
if(p==10){
  a3<-23.209
  a1<-4319000
  a2<-2.644
}
if(p==30){
  a3<-50.892
  a1<-14890000
  a2<-2.438
}
LC<-a3+(a1/(n^a2))

beta[r] <- round(mean(T_SR2[r] < LC),3)
}
uno_beta[t]<-mean(1-beta)
nocentralidad[t]<-daticos2$ncp
}


##### Funciones de periferia a partir de distintas funciones de profundidad
library(DepthProc)

daticos<-datos1(50,2,0.2,"SI")
daticos_nuevos<-datos2(50,2,"SI")

# Proyeccion
proyeccion<- 1/depth(daticos, method = "Projection") # es robusta

lc_p<-quantile(proyeccion, probs = 0.5) + 1.5*IQR(proyeccion)
plot(1:50,proyeccion, ylab = "Periferia basada en projection depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase I")
abline(h=lc_p, lwd = 1, lty = 2)
which(proyeccion>lc_p)

proyeccion_nuevos<-1/depth(u=daticos_nuevos,X=daticos, method = "Projection")
plot(1:50,proyeccion_nuevos, ylab = "Periferia basada en projection depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase II")
abline(h=lc_p, lwd = 1, lty = 2)

# Euclidiana
euclidiana<-1/depth(daticos, method = "Euclidean")

lc_e<-quantile(euclidiana, probs = 0.5) + 1.5*IQR(euclidiana)
plot(1:50,euclidiana, ylab = "Periferia basada en euclidean depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase I")
abline(h=lc_e, lwd = 1, lty = 2)
which(euclidiana>lc_e)

euclidiana_nuevos<-1/depth(u=daticos_nuevos,X=daticos, method = "Euclidean")
plot(1:50,euclidiana_nuevos, ylab = "Periferia basada en euclidean depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase II")
abline(h=lc_e, lwd = 1, lty = 2)

# Simplicial
library(ddalpha)
simplicial<-1/depth.simplicial(x=daticos, data=daticos)

lc_s<-quantile(simplicial, probs = 0.5) + 1.5*IQR(simplicial)
plot(1:50,simplicial, ylab = "Periferia basada en simplicial depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase I")
abline(h=lc_s, lwd = 1, lty = 2)
which(simplicial>lc_s)

simplicial_nuevos<-1/depth.simplicial(x=daticos_nuevos, data=daticos)
plot(1:50,tukey_nuevos, ylab = "Periferia basada en simplicial depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase II")
abline(h=lc_s, lwd = 1, lty = 2)

# Spatial

spatial<-1/depth.spatial(x=daticos_nuevos, data=daticos, mah.estimate = "MCD")

lc_sp<-quantile(spatial, probs = 0.5) + 1.5*IQR(spatial)
plot(1:50,spatial, ylab = "Periferia basada en spatial depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase I")
abline(h=lc_sp, lwd = 1, lty = 2)
which(spatial>lc_sp)

spatial_nuevos<-1/depth.simplicial(x=daticos_nuevos, data=daticos)
plot(1:50,spatial_nuevos, ylab = "Periferia basada en spatial depth", type = "b", pch = 20, 
     xlab = "Observaciones", main = "Carta de control para p=2 y contaminación=2 \n Fase II")
abline(h=lc_sp, lwd = 1, lty = 2)


##### beta #######

# Cuando las variables son independientes:

iter <- 1000
p<-30
n<-50

beta <- rep(0, iter)
changes<-c(0.2,0.5,1,1.2,1.5,2,2.5)
uno_beta_spatial<-rep(0,length(changes))
nocentralidad_spatial<-rep(0,length(changes))

set.seed(1234)
t<-0
for(change in changes){
  t<-t+1
for (i in 1:iter){
  daticos<-datos1(n,p,0.2,"SI",change)
  X<-daticos$datos
  spatial<-1/depth.spatial(x=X, data=X)
  lc_sp<-quantile(spatial, probs = 0.5) + 0.5*IQR(spatial)
  
  daticos_nuevos<-datos2(1,p,"SI",change)
  X2 <- daticos_nuevos$datos
  spatial_nuevos<-1/depth.spatial(x=X2, data=X)
  beta[i] <- round(mean(spatial_nuevos < lc_sp),3)
}
uno_beta_spatial[t] <- mean(1-beta)
nocentralidad_spatial[t]<-daticos_nuevos$ncp
}

par(mar = c(3,3,1,1), mgp = c(1.75,0.75,0))
plot(nocentralidad,uno_beta, main="p = 10",
     xlab = "Parámetro de no centralidad Fase II",
     ylab = "Potencia", type = "o", pch = 19, col = "red",
     ylim = c(0,1))
lines(nocentralidad_spatial,uno_beta_spatial, 
      type = "o", pch = 19, col = "blue")
legend("topleft", legend = c(expression(T[SR]^2), expression(O[Spatial])), 
       col = c('red','blue'), lty = 1, lwd = 2, bty = "n")



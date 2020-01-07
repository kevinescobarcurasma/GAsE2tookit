#' Run the program
#'
#' more detailed description
#'
#' @param none
#'
#' @return none
#'
#' @examples
#' Run_GAs()
#'
#' @export
Run_GAs<-function(){
creac_poblac<-function(number_ind1,genes1){
  poblac1<-matrix(0,number_ind1,genes1)
  for(i in 1:number_ind1){
    for(j in 1:genes1) {
      x<-runif(1)
      if(x[1]<0.5){
        poblac1[i,j]<-0 }
      else{
        poblac1[i,j]<-1 }
    }
  }
  return(poblac1)}
decodifica<-function(individuo1,NT11){
  genes1<-poblac[individuo1,]
  diametros1<-NULL
  for(i in 1:NT11){
    diametros1[i]<-switch(as.character(genes1[1+3*(i-1)]*2^2+genes1[2+3*(i-1)]*2+genes1[3+3*(i-1)]),"0"=4,"1"=6,"2"=8,"3"=10,"4"=12,"5"=14,"6"=16,"7"=18)
  }

  return(diametros1)}

simula_hidrau<-function(diametros1,NT11,NN11){
  for (i in 1:NT11){
    diam_mm<-diametros1[i]*25.4   # Convercion de pulg a mm para el diametro de al tuberia
    ENsetlinkvalue(i, 0, diam_mm)
  }
  ttt<-1
  tiempo=3600*0
  numnudos<-ENgetcount(0)
  presion1<-rep(0,NN11)
  ENopenH()
  ENinitH(0)
  while(ttt!=0){
    tt<-ENrunH()
    if(tiempo==tt){
      for(i in 1:numnudos){
        tiponudo  <- ENgetnodetype(i)
        if (tiponudo==0){ #
          presion1[i]<- ENgetnodevalue(i,"EN_PRESSURE")
        }
      }
    }
    ttt<-ENnextH()
  }
  ENcloseH()
  return(presion1)
}


fitness<-function(NT1,NN1,number_ind1){
  P_min=10
  Cp=20000000
  #presiones<-c(10.46,11.95,14.44,19.73,21.97,28.14,19.66)
  presiones_min<-rep(P_min,NN1)
  list_fa<-NULL
  for(i in 1:number_ind1){
    costo<-0
    jj<-1
    diametros<-decodifica(i,NT1)
    for(j in diametros){ #los diametros estasn en pulg
      costo<-costo+longitudes[jj]*switch(as.character(j),"4"=3.333,"6"=4.444,"8"=5.555,"10"=6.666,"12"=7.777,"14"=9.999,"16"=11.111,"18"=13.333)  #Precios en S/ segun al diametro por ml.
      jj<-jj+1}
    presiones<-simula_hidrau(diametros,NT1,NN1)
    list_fa[i]<-1/(costo+round(max(abs(presiones-presiones_min))*Cp,4))
  }

  return(list_fa)}

selecc<-function(number_ind1,list_fa){
  fa_1<-NULL
  jj<-1
  for(i in list_fa){
    fa_1[jj]<-i/sum(list_fa)
    jj<-jj+1
  }

  fa_acumu1<-cumsum(fa_1)
  list_selecc1<-NULL
  for(j in 1:number_ind1){
    jj<-1
    x<-runif(1)
    for (i in fa_acumu1){
      if(x[1]<i){
        list_selecc1[j]<-jj
        break
      }
      jj<-jj+1}}
  return(list_selecc1)}
cruza<-function(Pc1,genes1,number_ind1,list_selecc1){
  new_number<-number_ind1/2
  new_poblac<-matrix(0,number_ind1,genes1)
  for (i in 1:new_number){
    x<-runif(1)
    #print(x)
    if(x[1]<=Pc1){
      parent1<-poblac[list_selecc1[2*i-1],]
      parent2<-poblac[list_selecc1[2*i],]
      #print(parent1)
      #print(parent2)
      y<-round(runif(1,min=0,max=genes1-1),0)
      #print(y)
      descend1<-parent1
      descend2<-parent2
      descend1[(y[1]+1):genes1]<-parent2[(y[1]+1):genes1]
      descend2[(y[1]+1):genes1]<-parent1[(y[1]+1):genes1]
      new_poblac[2*i-1,]<-descend1
      new_poblac[2*i,]<-descend2
      #print(descend1)
      #print(descend2)
    }
    else{
      new_poblac[2*i-1,]<-poblac[list_selecc1[2*i-1],]
      new_poblac[2*i,]<-poblac[list_selecc1[2*i],]
      #print(new_poblac[2*i-1,])
      #print(new_poblac[2*i,])
    }
  }

  poblac<-new_poblac
  return(poblac)}

muta<-function(Pm1,genes1,number_ind1){
  for (i in 1:number_ind1){
    x<-runif(1)
    if (x[1]<=Pm1){
      y<-round(runif(1,min=1,max=genes1),0)
      if(poblac[i,y[1]]==0){
        poblac[i,y[1]]<-1 }
      else{
        poblac[i,y[1]]<-0 }
    }
  }
  return(poblac)}

library(epanet2toolkit)
inp <- file.path( find.package("epanet2toolkit"), "extdata","Net11.inp")
rpt <- file.path( find.package("epanet2toolkit"), "extdata","Net11.rpt")
ENopen( inp, rpt)
misdatos<-read.table("data.txt", header=TRUE, skip=1)
NN<-0
NT<-0
Pc<-misdatos["Pc"][1,1]
Pm<-misdatos["Pm"][1,1]
numnudos<-ENgetcount(0)
NT<-ENgetcount(2)
for (i in 1:numnudos){
  tiponodo<-ENgetnodetype(i)
  if (tiponodo==0){
    NN<-NN+1
  }
}
cromosomas<-3
genes<-NT*cromosomas
number_ind<-misdatos["I_N"][1,1]
generation_number<-misdatos["G_N"][1,1]
#longitudes<<-c(100,800,200,800,200) # las long estan en ml.
longitudes<-rep(0,NT)
for (i in 1:NT){
  longitudes[i]<-ENgetlinkvalue(i,1)
}
  poblac<-creac_poblac(number_ind,genes)
  print(poblac)
  YY<-rep(0,generation_number)
  for (generation in (1:generation_number)){
    print(generation)
    list_f_a<-fitness(NT,NN,number_ind)
    list_selecc<-selecc(number_ind,list_f_a)
    poblac<-cruza(Pc,genes,number_ind,list_selecc)
    poblac<-muta(Pm,genes,number_ind)
    YY[generation]<-mean(list_f_a)*1000000000 #Donde YY es una lista de los valores de la funcion aptitud para cada generacion
  }
  i<-which.max(list_f_a)
  diametros<-decodifica(i,NT)
  presiones<-simula_hidrau(diametros,NT,NN)
  print(diametros)
  print(presiones)
  N_N<-1:NN
  resultados<-data.frame(N_N, presiones)
  Para_grafica<-data.frame(1:generation_number, YY)
  write.table(resultados, file=paste(getwd(),"/Reporte_resultados.txt",sep=""), row.names=F )
  write.table(Para_grafica, file=paste(getwd(),"/valores_Grafica.txt",sep=""), row.names=F )
  #row.names=F
  ENclose()
}

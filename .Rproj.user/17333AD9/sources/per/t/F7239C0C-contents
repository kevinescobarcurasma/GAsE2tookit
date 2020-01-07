#' Initialize the program
#'
#' more detailed description
#'
#' @param none
#'
#' @return none
#'
#' @examples
#' Open_epanet2toolkit()
#'
#' @export
Open_epanet2toolkit<-function(){
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


}

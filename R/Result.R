#' Result report
#'
#' more detailed description
#'
#' @param none
#'
#' @return none
#'
#' @examples
#' Results_report()
#'
#' @export
Results_report<-function(){
  write.table(resultados, file=paste(getwd(),"/Reporte_resultados.txt",sep=""), row.names=F )
  write.table(Para_grafica, file=paste(getwd(),"/valores_Grafica.txt",sep=""), row.names=F )
  #row.names=F
  ENclose()
}

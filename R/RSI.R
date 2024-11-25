

#' Radiosensitivity Index
#'
#' @param x An gene expression matrix with rownames<-gene symbol and colnames<-sample ID
#'
#' @return RSI
#' @export
#'
#' @examples
#' @export
calRSI <- function(x) {
  RSIlist <-c("AR",
              "JUN",
              "STAT1",
              "PRKCB",
              "RELA",
              "ABL1",
              "SUMO1",
              "CDK1",
              "HDAC1",
              "IRF1")
  dataRSI<-x[RSIlist,]
  RSI<- NULL
  for (i in colnames(dataRSI) ) {
    RANKi<-rank(dataRSI[,i])
    RSIi<- -0.0098009*RANKi[1]+0.0128283*RANKi[2]+0.0254552*RANKi[3]-0.0017589*RANKi[4]-0.0038171*RANKi[5]+0.1070213*RANKi[6]-0.0002509*RANKi[7]-0.0092431*RANKi[8]-0.0204469*RANKi[9]-0.0441683*RANKi[10]
    outab<-rbind(i,RSIi)
    RSI <-cbind(RSI,outab)
  }
}



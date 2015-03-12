#'Get a list with two named items "x" and "survModel" (the return value of encapsulating a call to survreg in tryCatch'.)
#''
#''
#' @param Sam is a copy of a dataframe created for EGRET analysis with only records for items with non-zero weights.
#' @param weight are the weights associated with each record in Sam
#''
#' @import survival
#''
#' @return either NA, NULL, or a list with two named items "x" and "survModel" (the return value of encapsulating a call to survreg in tryCatch'). "x" and "survModel" are both normally of the class "survreg".
#''
#' @export
#''
#' @examples
#''
#GetX_and_survModel<-function(ConcLow, ConcHigh, DecYear, LogQ, SinDY, CosDY, Sam,weight){
GetX_and_survModel<-function(Sam, weight){

    x <- tryCatch({
      survModel<-survreg(Surv(log(ConcLow),log(ConcHigh),type="interval2") ~ 
                           DecYear+LogQ+SinDY+CosDY,data=Sam,weights=weight,dist="gaus")
#      return(list("x"=x,"survModel"=survModel))
      

    }, warning=function(w) {
      return(NA)
    }, error=function(e) {
      message(e, "Error")
      return(NULL)
    })
    if(exists("survModel")) {
      return(list("x"=x,"survModel"=survModel))
    }
}
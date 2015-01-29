#' Estimation process for the WRTDS (Weighted Regressions on Time, Discharge, and Season)
#'
#' This one function does a jack-knife cross-validation of a WRTDS model, fits the surface
#' (concentration as a function of discharge and time), 
#' estimates daily values of concentration and flux, and flow normalized values. 
#' It returns a named list with the following dataframes: Daily, INFO, Sample, and the matrix: surfaces.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param windowY numeric specifying the half-window width in the time dimension, in units of years, default is 7
#' @param windowQ numeric specifying the half-window width in the discharge dimension, units are natural log units, default is 2
#' @param windowS numeric specifying the half-window with in the seasonal dimension, in units of years, default is 0.5
#' @param minNumObs numeric specifying the miniumum number of observations required to run the weighted regression, default is 100
#' @param minNumUncen numeric specifying the minimum number of uncensored observations to run the weighted regression, default is 50
#' @param edgeAdjust logical specifying whether to use the modified method for calculating the windows at the edge of the record.  
#' The modified method tends to reduce curvature near the start and end of record.  Default is TRUE.
#' @keywords water-quality statistics
#' @export
#' @return eList named list with Daily, Sample, and INFO dataframes, along with the surfaces matrix.
#' Any of these values can be NA, not all EGRET functions will work with missing parts of the named list eList.
#' @examples
#' eList <- Choptank_eList
#' \dontrun{EGRETreturn <- modelEstimation(eList)
#' Daily <- EGRETreturn$Daily
#' Sample <- EGRETreturn$Sample
#' INFO <- EGRETreturn$INFO
#' surfaces <- EGRETreturn$surfaces
#' }
modelEstimation<-function(eList, 
                          windowY=7, windowQ=2, windowS=0.5,
                          minNumObs=100,minNumUncen=50, 
                          edgeAdjust=TRUE){

GetWeightsSource <- '
  Rcpp::DataFrame internal_localSample = Rcpp::as<Rcpp::DataFrame>(localSample);
  Rcpp::CharacterVector ColumnNames = internal_localSample.attr("names");
  int DecYearCol = -1;
  int LogQCol = -1;
  int UncenCol = -1;
  int ColumnCount = ColumnNames.length();
  for(int ColIndex = 0; ColIndex < ColumnCount; ColIndex++) {
    if (ColumnNames[ColIndex] == "LogQ") {
      LogQCol = ColIndex;
    }
    if (ColumnNames[ColIndex] == "DecYear") {
      DecYearCol = ColIndex;
    }
    if (ColumnNames[ColIndex] == "Uncen") {
      UncenCol = ColIndex;
    }
  }
  if ((DecYearCol < 0) || (LogQCol < 0) || (UncenCol < 0)) {
    throw std::invalid_argument( "DecYear, LogQ, or Uncen were not columns in the data frame." );
  }
  Rcpp::NumericVector LogQ = Rcpp::as<Rcpp::NumericVector>(internal_localSample[LogQCol]);
  Rcpp::NumericVector Uncen = Rcpp::as<Rcpp::NumericVector>(internal_localSample[UncenCol]);
  Rcpp::NumericVector DecYear = Rcpp::as<Rcpp::NumericVector>(internal_localSample[DecYearCol]);
  Rcpp::CharacterVector RowNames = internal_localSample.attr("row.names");
//  double internal_estY = Rcpp::as<double>(estY);
  double internal_WindowS = Rcpp::as<double>(windowS); 
  double internal_WindowQ = Rcpp::as<double>(windowQ); 
  double internal_WindowY = Rcpp::as<double>(windowY); 
  double internal_minNumObs = Rcpp::as<double>(minNumObs); 
  double internal_minNumUncen = Rcpp::as<double>(minNumUncen); 
  Rcpp::NumericVector Internal_estPtLQ = Rcpp::as<Rcpp::NumericVector>(estPtLQ);
  bool internal_edgeAdjust = Rcpp::as<bool>(edgeAdjust);
  Rcpp::NumericVector internal_estPtYear = Rcpp::as<Rcpp::NumericVector>(estPtYear);
  double internal_DecLow = Rcpp::as<double>(DecLow); 
  double internal_DecHigh = Rcpp::as<double>(DecHigh); 
  int internal_numEstPt = internal_estPtYear.length();

  // In C++, arrays start at 0 whereas in R they start at 1.
//  int internal_i = Rcpp::as<int>(i)-1;

  Rcpp::List WeightList(internal_numEstPt);

  for(int internal_i = 0; internal_i < internal_numEstPt; internal_i++) {
    Rcpp::NumericVector weightsVector = Rcpp::clone(LogQ);
    double internal_tempWindowY = internal_WindowY;
    double internal_tempWindowQ = internal_WindowQ; 
    double internal_tempWindowS = internal_WindowS; 

    double internal_estY = internal_estPtYear[internal_i];

    double distLow = internal_estY-internal_DecLow;
    double distHigh = internal_DecHigh-internal_estY;
    double distTime = std::min(distLow,distHigh);
    if (internal_edgeAdjust) {
      // I think this should be if (distTime < internal_tempWindowY) 
      // (less than) instead of (less than or equal). However,
      // <= captures the behavior of the original R code.
      if (distTime <= internal_tempWindowY) {
        internal_tempWindowY = ((2 * internal_tempWindowY) - distTime);
      }
    }

    double internal_estLQ = Internal_estPtLQ[internal_i];

    int Count = LogQ.length();
    do {

      int numPosWt = 0;
      int numUncen = 0;
      for(int it = 0; it < Count; it++) {
        if (DecYear[it]-internal_estY <= internal_tempWindowY){
          double diffY = fabs(DecYear[it]-internal_estY);
          double weightY = tricube(diffY,internal_tempWindowY);
          double weightQ = tricube(LogQ[it]-internal_estLQ,internal_tempWindowQ);
          double diffUpper = ceil(diffY);
          double diffLower = floor(diffY);
          double diffSeason = std::min(fabs(diffUpper-diffY),fabs(diffY-diffLower));
          double weightS = tricube(diffSeason,internal_tempWindowS);
          double weight = weightY*weightQ*weightS;
          weightsVector[it] = weight; 
          if (weight > 0) {
            numPosWt++;
            numUncen = numUncen + Uncen[it];
          }
        }
        else {
          weightsVector[it] = 0; 
        }
      
      }
      internal_tempWindowY = internal_tempWindowY*1.1;
      internal_tempWindowQ = internal_tempWindowQ*1.1;
  //  the next line is designed so that if the user sets windowS so it includes
  //  data from all seasons, the widening process leaves it alone  
      if (internal_WindowS<=0.5){
        internal_tempWindowS = std::min(internal_tempWindowS*1.1,0.5);
      }
      else
      {
        internal_tempWindowS = internal_WindowS;
      }
      if ((numPosWt>=internal_minNumObs) && (numUncen>=internal_minNumUncen)) {
        break;	
      }
    } while(1);
//  return weightsVector;
    WeightList[internal_i] = weightsVector;
  }
  return WeightList;
'

srcTriCube = '
double tricube(double d, double h)
{
  double ad = fabs(d);
  double w = pow(1- pow(ad/h,3),3);
  if (w < 0)
  {
    w = 0;
  }
  return w;
};
'

GetWeights <- cxxfunction(signature(estPtYear="numeric", estPtLQ="numeric",   
  DecLow="numeric", DecHigh="numeric", localSample="data.frame", 
  windowY="numeric", windowQ="numeric",
  windowS="numeric", minNumObs="numeric", minNumUncen="numeric",
  edgeAdjust="logical"), 
  GetWeightsSource, include = srcTriCube, plugin="Rcpp")

  # this code is a wrapper for several different functions that test the model, fit a surface,
  #  estimate daily values and flow normalized daily values
  #  and organize these into monthly results
  #  it returns several data frames
  #  all of the data frames are given their "standard" names
  #
  localINFO <- getInfo(eList)
  localSample <- getSample(eList)
  localDaily <- getDaily(eList)
  
  numDays <- length(localDaily$DecYear)
  DecLow <- localDaily$DecYear[1]
  DecHigh <- localDaily$DecYear[numDays]
  
  cat("\n first step running estCrossVal may take about 1 minute")
  Sample1<-estCrossVal(numDays,DecLow,DecHigh, localSample, GetWeights, 
                       windowY, windowQ, windowS, minNumObs, minNumUncen,
                       edgeAdjust)

  surfaceIndexParameters<-surfaceIndex(localDaily)
  localINFO$bottomLogQ<-surfaceIndexParameters[1]
  localINFO$stepLogQ<-surfaceIndexParameters[2]
  localINFO$nVectorLogQ<-surfaceIndexParameters[3]
  localINFO$bottomYear<-surfaceIndexParameters[4]
  localINFO$stepYear<-surfaceIndexParameters[5]
  localINFO$nVectorYear<-surfaceIndexParameters[6]
  localINFO$windowY<-windowY
  localINFO$windowQ<-windowQ
  localINFO$windowS<-windowS
  localINFO$minNumObs<-minNumObs
  localINFO$minNumUncen<-minNumUncen
  localINFO$numDays <- numDays
  localINFO$DecLow <- DecLow
  localINFO$DecHigh <- DecHigh
  localINFO$edgeAdjust <- edgeAdjust
  
  cat("\nNext step running  estSurfaces with survival regression:\n")
  surfaces1<-estSurfaces(eList, GetWeights, 
                         windowY, windowQ, windowS, minNumObs, minNumUncen, edgeAdjust)

  eList <- as.egret(Daily=localDaily, 
                    Sample=Sample1,
                    INFO=localINFO,
                    surfaces=surfaces1)
  
  Daily1<-estDailyFromSurfaces(eList)
  
  eList <- as.egret(Daily=Daily1, 
               Sample=Sample1,
               INFO=localINFO,
               surfaces=surfaces1)
  
  return(eList)
  
}

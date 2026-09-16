#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector GARCH_ES(NumericVector THETA, NumericVector y ,double empiricalES, NumericVector ES, int RowsOfy, int esPredict)
{
  int i;
  
  /* Initialize output variables */
  ES[0] = empiricalES;
  
  /* Start the loop */
  for(i = 1; i < RowsOfy; i++)
  {
    // Symmetric Absolute Value
    ES[i] = -sqrt(THETA[0] + THETA[1]*pow(ES[i-1],2) + THETA[2]*pow(y[i-1],2));
    
  }
  
  if (esPredict == 1){
    ES[RowsOfy] = -sqrt(THETA[0] + THETA[1]*pow(ES[RowsOfy-1],2) + THETA[2]*pow(y[RowsOfy-1],2));
  } 
  return ES;
}





#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector AS_ES(NumericVector THETA, NumericVector y ,double empiricalES, NumericVector ES, int RowsOfy, int esPredict)
{
  int i;
  
  /* Initialize output variables */
  ES[0] = empiricalES;
  
  /* Start the loop */
  for(i = 1; i < RowsOfy; i++)
  {
    // Symmetric Absolute Value
    ES[i] = THETA[0] + THETA[1] * ES[i-1] + THETA[2] * y[i-1]*(y[i-1]>0) - THETA[3]*y[i-1]*(y[i-1]<0);
    
  }
  
  if (esPredict == 1){
    ES[RowsOfy] = THETA[0] + THETA[1] * ES[RowsOfy-1] + THETA[2] * y[RowsOfy-1]*(y[RowsOfy-1]>0) - THETA[3]*y[RowsOfy-1]*(y[RowsOfy-1]<0);
  } 
  return ES;
}





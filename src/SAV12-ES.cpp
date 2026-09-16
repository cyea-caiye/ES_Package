#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector SAV12_ES(NumericVector THETA, NumericVector y ,double empiricalES, NumericVector ES, int RowsOfy, int esPredict)
{
  int i;
  
  /* Initialize output variables */
  ES[0] = empiricalES;
  ES[1] = THETA[0] + THETA[1]*ES[0] + THETA[2]*(y[0]*(y[0]>0) - y[0]*(y[0]<0));
  /* Start the loop */
  for(i = 2; i < RowsOfy; i++)
  {
    // Symmetric Absolute Value
    ES[i] = THETA[0] + THETA[1] * ES[i-1] + THETA[2] * (y[i-1]*(y[i-1]>0) - y[i-1]*(y[i-1]<0)) + THETA[3] * (y[i-2]*(y[i-2]>0) - y[i-2]*(y[i-2]<0));
    
  }
  
  if (esPredict == 1){
      ES[RowsOfy] = THETA[0] + THETA[1] * ES[RowsOfy-1] + THETA[2] * (y[RowsOfy-1]*(y[RowsOfy-1]>0) - y[RowsOfy-1]*(y[RowsOfy-1]<0)) + THETA[3] * (y[RowsOfy-2]*(y[RowsOfy-2]>0) - y[RowsOfy-2]*(y[RowsOfy-2]<0));
  } 
  return ES;
}





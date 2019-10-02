#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix multinomialSupport(int n, int k) {

  int size = Rf_choose(n+k-1,k-1);
  int p = 0;
  NumericMatrix mat(size, k);
  int j = 1;

  for(int i=0; i < n; i++)
    {
    p = Rf_choose(i+k-2, k-2);
    NumericMatrix::Sub sub = m(Range(j,j+p-1), Range(0,k-2));
    sub = multinomialSupport(k - 1, i);
    j = j + p;
    }

  mat( _ , k-1) = n - rowSums(mat);

  return mat;

}



#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Support of a multinomial distribution
//'
//' @param k The number of categories.
//' @param n The sample size.
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
arma::Mat<int> multinomial_support(int k, int n)
{

  if(k == 1)
  {
    arma::Mat<int> mat(1, 1);
    mat(0,0) = n;
    return(mat);
  }

  int size = Rf_choose(n+k-1,k-1);
  int p;
  arma::Mat<int> mat(size, k, fill::zeros);
  int j = 0;

  for(int i=0; i < n; i++)
  {
    p = Rf_choose(i+k-1, k-2);
    mat.submat(j+1, 0, j+p, k-2) = multinomial_support(k - 1, i+1);
    j = j + p;
  }

  mat.col(k-1) = n - arma::sum(mat, 1);

  return(mat);

}




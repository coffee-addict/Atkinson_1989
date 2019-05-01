#include <bits/stdc++.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
namespace blas = boost::numeric::ublas;

#define FOR(i,s,e) for(int i=(s);(i)<(int)(e);(i)++)
#define REP(i,e) FOR(i,0,e)
#define each(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();it++)
#define all(o) (o).begin(), (o).end()
#define pb(x) push_back(x)
#define mp make_pair
#define mt make_tuple
#define t0(t) get<0>((t))
#define t1(t) get<1>((t))
#define t2(t) get<2>((t))

typedef long long ll;

int main() {
  int n = 3;
  blas::matrix<double> A(n,n);  
  blas::vector<double> x(n);  //results (initially b in Ax=b)
  blas::permutation_matrix<double> p(n);
puts("A(before):");
  A(0,0) = -5;
  A(0,1) = 3;
  A(0,2) = 4;
  A(1,0) = 10;
  A(1,1) = -8;
  A(1,2) = -9;
  A(2,0) = 15;
  A(2,1) = 1;
  A(2,2) = 2;
  x(0) = 1;
  x(1) = 2;
  x(2) = 3;
  REP(i,n) {
    REP(j,n) printf("%.5f ", A(i,j));
    puts("");
  }
  blas::lu_factorize(A, p);
  lu_substitute(A,p,x);
puts("A(after):");
  REP(i,n) {
    REP(j,n) printf("%.5f ", A(i,j));
    puts("");
  }
  puts("p:");
  REP(i,n) printf("%.5f ", p(i));
  puts("");
  puts("x:");
  REP(i,n) printf("%.10f ", x(i));
  puts("");

  return 0;
}

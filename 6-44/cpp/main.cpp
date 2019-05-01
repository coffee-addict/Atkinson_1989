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

const double PI = acos(-1.0);
const int M = 16;
const int T = 5;

double calc_u(double t, double x) {
  return exp(-t*0.1) * sin(x*PI);
}

int main () {
  vector<vector<double>> res;
  res.pb(vector<double>({1, 2, 3, 4, 5}));
  int m = 4;
  while (m <= M) {
    double dx = 1.0/m;
    double dt = dx;
    double c = PI*PI*5*dt/(dx*dx);
    blas::vector<double> u_prev(m+1), x(m+1);
    u_prev[0] = u_prev[m] = 0;
    FOR(i,1,m) {
      x[i] = dx*i;
      u_prev[i] = calc_u(0,x[i]);
    }
    //mat_a u_prev = mat_b x
    blas::matrix<double> mat_a(m+1, m+1), mat_b(m+1, m+1);
    mat_a(0,0) = mat_a(m,m) = 1;
    mat_b(0,0) = mat_b(m,m) = 1;
    FOR(i,1,m) {
      if (i > 0) mat_a(i,i-1) = -c;
      if (i < m) mat_a(i,i+1) = -c;
      mat_a(i,i) = c*2 - 1;
      if (i > 0) mat_b(i,i-1) = c;
      if (i < m) mat_b(i,i+1) = c;
      mat_b(i,i) = -(c*2 + 1);
    }
//printf("m = %d, mat_b:\n", m);
//puts("mat_b:");
//REP(i,m+1) {
//  REP(j,m+1) printf("%.3f ", mat_b(i,j));
//  puts("");
//}
    //LU factorization for mat_b
    blas::permutation_matrix<double> p(m+1);
    blas::lu_factorize(mat_b, p);
    //iterate & solve x
    vector<double> tmp;
    double t = 0;
    int i = 0;
    while (i <= T*m) {
      x = blas::prod(mat_a, u_prev);
      blas::lu_substitute(mat_b,p,x);
      u_prev = x;
      i++;
      t = dt*i;
      //store results
      if (i%m == 0) {
        double error = 0;
        REP(j,m-1) error = max(error, abs(calc_u(t,(double)(j+1)/m) - x[j+1]));
        tmp.pb(error);
      }
    }
    res.pb(tmp);
    m *= 2;
  }
  //output
  puts("t,err(m=4),err(m=8),err(m=16)");
  REP(j,res[0].size()) REP(i,res.size())
    printf("%.5f%c", res[i][j], res.size()==i+1 ? '\n' : ' ');
}

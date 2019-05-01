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

const double EPS = 1e-5;
const double EPS_ZERO = 1e-12;
const int K = 1000;

//the format of input data
  //n m s
  //a11 a12 ... a1m
  //...
  //an1 an2 ... anm
  //b1 b2 ... bm

//A is assumed symmetory and positive definite in this program

double norm_infinity(const vector<double> &v) {
  double res = 0;
  REP(i,v.size()) res = max(res, abs(v[i]));
  return res;
}

double prd_inn(const vector<double> &u, const vector<double> &v) {
  double res = 0;
  REP(i,u.size()) res += u[i]*v[i];
  return res; 
}

vector<double> prd_mat_v(const vector<vector<double>> &mat, 
                         const vector<double> &v) {
  int n = (int)mat.size();
  int m = (int)v.size();
  vector<double> res(n);
  REP(i,n) REP(j,m) res[i] += mat[i][j]*v[j];
  return res;
}

vector<double> set_result(const vector<vector<double>> &mat,
  const vector<double> &r, const blas::vector<double> &x, 
  const vector<double> &sim) {
  vector<double> res;
  res.pb(norm_infinity(r));
  int n = (int)x.size();
  vector<double> vd(x.size());
  REP(i,n) vd[i] = x(i) - sim[i];
  res.pb(norm_infinity(vd));
  res.pb(sqrt(prd_inn(vd, prd_mat_v(mat, vd))));
  return res;
}

int main () {
  int n, m, s;
  scanf("%d%d%d ", &n, &m, &s);
  vector<vector<double>> mat, sim(2), p(2), r(2);
  vector<double> b(m);
  mat.resize(n); REP(i,n) mat[i].resize(m);
  REP(i,2) {
    sim[i].resize(m);
    p[i].resize(m);
    r[i].resize(m);
  }
  REP(i,n) REP(j,m) scanf("%lf ", &(mat[i][j]));
  REP(i,m) scanf("%lf ", &(b[i]));
  //use ublas to get results of direct methods(LU factorization)
  blas::vector<double> x(m);
  blas::matrix<double> mat_lu(n,m);
  REP(i,n) REP(j,m) mat_lu(i,j) = mat[i][j];
  REP(i,m) x(i) = b[i];
  blas::permutation_matrix<double> pm(m);
  blas::lu_factorize(mat_lu, pm);
  blas::lu_substitute(mat_lu,pm,x);
  //Conjugate gradient method
  sim[0][0] = 1;
  vector<vector<double>> res;
  vector<double> mat_p0 = prd_mat_v(mat, sim[0]);
  REP(i,m) p[0][i] = r[0][i] = b[i] - mat_p0[i];
  int converged = 0;
  double eps_exit = norm_infinity(b);
  REP(k,K) {
    int pk = k&1, nk = (k+1)&1;
    double a = prd_inn(r[pk], p[pk]);
    vector<double> mat_p = prd_mat_v(mat, p[pk]);
    double denom = prd_inn(p[pk], mat_p);
    if (abs(denom) < EPS_ZERO) {
      printf("%dth iteration, 0 division detected computing alpha.\n", k+1);
      break;
    }
    a /= denom;
    REP(i,m) sim[nk][i] = sim[pk][i] + a*p[pk][i];
    REP(i,m) r[nk][i] = r[pk][i] - a*mat_p[i];
    double b = - prd_inn(r[nk], mat_p);
    denom = prd_inn(p[pk], mat_p);
    if (abs(denom) < EPS_ZERO) {
      printf("%dth iteration, 0 division detected computing beta.\n", k+1);
      break;
    }
    b /= denom;
    REP(i,m) p[nk][i] = r[nk][i] + b*p[pk][i];
    res.pb(set_result(mat,r[nk],x,sim[nk]));
    if (res[k][0] < eps_exit*EPS) {
printf("eps_expt*EPS = %f\n", eps_exit*EPS);
      converged = 1;
      break;
    }
  }
  puts("k,|r_k|_inf,|x-x_k|_inf,|x-x_k|_A");
  REP(k,res.size()) {
    printf("%d", (int)(k+1));
    REP(j,res[k].size()) printf(",%.3e", res[k][j]);
    puts("");
  }
  puts("Solutions:");
  printf("LU: ");
  REP(i,x.size()) printf("%.5f ", x[i]);
  puts("");
  printf("CG: ");
  REP(i,sim[(sim.size()+1)&1].size())
    printf("%.5f ", sim[(sim.size()+1)&1][i]);
  puts("");
  if (converged)
    printf("converged by CG method at %dth iteration\n", (int)res.size());
  else
    printf("does not converse by CG method in %d iterations.\n", K);
  return 0;
}

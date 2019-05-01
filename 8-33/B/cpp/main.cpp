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

const double EPS = 1e-10;
const double PI = acos(-1.0);
const int M = 16;
const int K = 1000;
double tu[M+1][M+1], sim[2][M+1][M+1], inorm[2][K];

double calc_u(double x, double y) {
  return exp(x*PI) * cos(y*PI);
}

int main () {
  //true values
  REP(i,M+1) REP(j,M+1) tu[i][j] = calc_u(1.0*j/M, 1.0*i/M);
  //initial data on the boundaries
  REP(k,2) REP(j,M+1) {
    sim[k][j][0] = tu[j][0];
    sim[k][0][j] = tu[0][j];
    sim[k][j][M] = tu[j][M];
    sim[k][M][j] = tu[M][j];
  }
  //Gauss-Seidel method
  REP(i,2) REP(k,K) inorm[i][k] = -1;
  int converged = 0;
  int pi = 0, ni = 1;
  REP(k,K) {
    pi = k&1, ni = (k+1)&1;
    FOR(i,1,M) FOR(j,1,M) {
      sim[ni][i][j] = sim[ni][i-1][j];
      sim[ni][i][j] += sim[pi][i+1][j];
      sim[ni][i][j] += sim[ni][i][j-1];
      sim[ni][i][j] += sim[pi][i][j+1];
      sim[ni][i][j] *= 0.25;
    }
    inorm[0][k] = inorm[1][k] = 0;
    FOR(i,1,M) FOR(j,1,M) {
      inorm[0][k] = max(inorm[0][k], abs(tu[i][j]-sim[ni][i][j]));
      inorm[1][k] = max(inorm[1][k], abs(sim[ni][i][j]-sim[pi][i][j]));
    }
    if (inorm[0][k] < EPS || (k && abs(inorm[0][k]-inorm[0][k-1]) < EPS)) {
      converged = 1;
      break;
    }
  }
  //output
  printf("k,max|u-u^(k)|,max|u^(k+1)-u^(k)|,");
  puts("max|u^(k+1)-u^(k)|/max|u^(k)-u^(k-1)|");
  printf("%d,%.3e,%.3e,%.3e\n", 1, inorm[0][0], 0.0, 0.0);
  FOR(k,1,K) {
    if (inorm[0][k] < 0) break;
    double a, b;
    a = abs(inorm[0][k] - inorm[0][k-1]);
    if (k>1) b = abs(inorm[0][k-1] - inorm[0][k-2]);
    printf("%d,%.3e,%.3e,%.3e\n", k+1, inorm[0][k], a, k>1 ? a/b : 0);
  }
  if (!converged) printf("does not converse in %d iterations.\n", K);
  return 0;
}

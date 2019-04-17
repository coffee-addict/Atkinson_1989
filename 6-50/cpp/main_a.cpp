#include <bits/stdc++.h>

using namespace std;

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
const int NF = 3; //# of kinds of functions
const int M = 3;  //# of kinds of step sizes
const int N = 500; //(10-0)/.125
int n[M] = {10, 50, 500};
double h;
double a[NF] = {0, 0, 0};
double b[NF] = {5, 5, 5}; 
double lipschitz_c[NF] = {-1, -10, 50};
double lambda[NF] = {-1, -10, -50};
double xs[N+1], ys[N+1];
double res[3][N+1]; //y_h, |y_h - y|, bound(|y(t_n)-y_h(t_n)|)

//y(x, y) 
double fn_y_d0(int i_f, double x) {
  return sin(x) + cos(x);
}

//y' = f(x,y)
double fn_y_d1(int i_f, double x, double y) {
  double res = lambda[i_f]*y;
  res += cos(x)*(1.0-lambda[i_f]);
  res += sin(x)*(1.0+lambda[i_f]);
  return res;
}

void init(int i_f, int i_div) {
  h = (b[i_f]-a[i_f])/n[i_div];
  memset(xs, 0, sizeof(xs));
  memset(ys, 0, sizeof(ys));
  memset(res, 0, sizeof(res));
  REP(k,n[i_div]) {
    xs[k] = a[i_f] + h*k;
    ys[k] = fn_y_d0(i_f, xs[k]);
  }
  xs[n[i_div]] = b[i_f];
  ys[n[i_div]] = fn_y_d0(i_f, b[i_f]);
}

void calc_ode_runge_kutta(int i_f, int i_div) {
  res[0][0] = ys[0];
  REP(i,n[i_div]) {
    double v1 = fn_y_d1(i_f, xs[i], res[0][i]);
    double v2 = fn_y_d1(i_f, xs[i]+h*0.5, res[0][i]+v1*h*0.5);
    double v3 = fn_y_d1(i_f, xs[i]+h*0.5, res[0][i]+v2*h*0.5);
    double v4 = fn_y_d1(i_f, xs[i]+h, res[0][i]+v3*h);
    double v = (v1+v2*2+v3*2+v4);
    res[0][i+1] = res[0][i] + v*h/6;
    res[1][i+1] = abs(ys[i+1] - res[0][i+1]);
    res[2][i] = (ys[i+1]-ys[i])/h;
    res[2][i] -= v/6;
    res[2][i] = i ? max(abs(res[2][i]), res[2][i-1]) : abs(res[2][i]);
  }
  REP(i,n[i_div]) {
    res[2][i] *= exp(lipschitz_c[i_f]*(xs[i]-xs[0])-1);
    res[2][i] /= lipschitz_c[i_f];
  }
  res[2][n[i_div]] = res[2][n[i_div]-1];
}

void output(int i_f, int i_div) {
  puts("func,h,x,y_h,y,y-y_h,bound");
  REP(i,n[i_div]+1) {
    int skip = 1;
    REP(j,5) {
      double diff = abs(xs[i]-(double)(j+1));
      if (diff < 0.001) skip = 0;
    }
    if (skip) continue;
    printf("%d,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n", 
           i_f, h, xs[i], res[0][i], ys[i], res[1][i], res[2][i]);
  }
}

int main() {
  REP(i,NF) REP(j,M) {
    init(i,j);
    calc_ode_runge_kutta(i,j);
    output(i,j);
  }
  return 0;
}

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
const int NF = 1; //# of kinds of functions
const int M = 3;  //# of kinds of step sizes
const int N = 64; //(4-0)/0.0625
int n[M] = {16, 32, 64};
double h;
double a[NF] = {0};
double b[NF] = {4}; 
double lipschitz_c[NF] = {1};
double xs[N+1], ys[N+1], ys_d2[N+1];
double res[3][N+1]; //y_h, |y_h - y|, bound(|y(t_n)-y_h(t_n)|)

double f(double x, double y) {
  return x*x-y;
}

double y(double x) {
  return x*x - x*2 + 2 - exp(-x);
}

double y_d2(double x) {
  return 2.0 - exp(-x);
}

double (*f_p[])(double, double) = {f};
double (*y_p[])(double) = {y};
double (*y_d2_p[])(double) = {y_d2};

void init(int i, int j) {
  h = (b[i]-a[i])/n[j];
  memset(xs, 0, sizeof(xs));
  memset(ys, 0, sizeof(ys));
  memset(ys_d2, 0, sizeof(ys_d2));
  memset(res, 0, sizeof(res));
  REP(k,n[j]) {
    xs[k] = a[i] + h*k;
    ys[k] = y_p[i](xs[k]);
    ys_d2[k] = y_d2_p[i](xs[k]);
  }
  xs[n[j]] = b[i];
  ys[n[j]] = y_p[i](b[i]);
  ys_d2[n[j]] = y_d2_p[i](b[i]);
}

void calc_euler(int i, int j) {
  res[0][0] = ys[0];
  REP(k,n[j]) {
    res[0][k+1] = res[0][k] + h*f_p[i](xs[k],ys[k]);
    res[1][k+1] = ys[k+1] - res[0][k+1];
    res[2][k+1] = max(abs(ys_d2[k+1]), abs(ys_d2[k]))*h*0.5;
    res[2][k+1] *= exp(lipschitz_c[i]*(xs[k+1]-xs[0])) - 1;
    res[2][k+1] /= lipschitz_c[i];
  }
}

void output(int i, int j) {
  puts("func,h,x,y_h,y,|y - y_h|,bound");
  REP(k,n[j]+1) 
    printf("%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n", 
           i, h, xs[k], res[0][k], ys[k], res[1][k], res[2][k]);
}

int main() {
  REP(i,NF) REP(j,M) {
    init(i,j);
    calc_euler(i,j);
    output(i,j);
  }
  return 0;
}

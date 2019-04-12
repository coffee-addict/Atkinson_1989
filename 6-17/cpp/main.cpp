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
const int NF = 4; //# of kinds of functions
const int M = 2;  //# of kinds of step sizes
const int N = 40; //(10-0)/.25
const int I = 3;  //max # of iterations
int n[M] = {20, 40};
double h;
double a[NF] = {0, 0, 0, 0};
double b[NF] = {10, 10, 10, 10}; 
double lipschitz_c[NF] = {20, 0.25, 1, 1};
double xs[N+1], ys[N+1], ys_d3[N+1];
double res[3][N+1]; //y_h, |y_h - y|, bound(|y(t_n)-y_h(t_n)|)

double f(int i_f, double x, double y) {
  double res = 0;
  switch(i_f) {
    case 0:
      res = -y*y;
      break;
    case 1:
      res = y*0.25*(1.0 - y*0.05);
      break;
    case 2:
      res = -y + cos(x)*2;
      break;
    case 3:
      res = y - sin(x)*2;
      break;
    default:
      res = -9999999999;
      break;
  }
  return res;
}

double fn_y_d0(int i_f, double x) {
  double res = 0;
  switch(i_f) {
    case 0:
      res = 1.0/(1.0 + x);
      break;
    case 1:
      res = 20.0/(1.0 + exp(-x*0.25)*19);
      break;
    case 2:
      res = cos(x) + sin(x);
      break;
    case 3:
      res = cos(x) + sin(x);
      break;
    default:
      res = -9999999999;
      break;
  }
  return res;
}

double fn_y_d3(int i_f, double x) {
  double res = 0, tmp = 0;
  switch(i_f) {
    case 0:
      res = 1;
      REP(i,4) res *= (1.0 + x);
      res = -6.0/res;
      break;
    case 1:
      tmp = 1;
      res = 361.0 + exp(x*0.5) - exp(x*0.25)*76;
      res *= exp(x*0.25)*95;
      REP(i,4) tmp *= (exp(x*0.25) + 19);
      res /= 16*tmp;
      REP(i,4) res *= (1.0 + x);
      res = -6.0/res;
      break;
    case 3:
      res = -cos(x) + sin(x);
      break;
    case 4:
      res = -cos(x) + sin(x);
      break;
    default:
      res = -9999999999;
      break;
  }
  return res;
}

void init(int i_f, int i_div) {
  h = (b[i_f]-a[i_f])/n[i_div];
  memset(xs, 0, sizeof(xs));
  memset(ys, 0, sizeof(ys));
  memset(ys_d3, 0, sizeof(ys_d3));
  memset(res, 0, sizeof(res));
  REP(k,n[i_div]) {
    xs[k] = a[i_f] + h*k;
    ys[k] = fn_y_d0(i_f, xs[k]);
    ys_d3[k] = fn_y_d3(i_f, xs[k]);
  }
  xs[n[i_div]] = b[i_f];
  ys[n[i_div]] = fn_y_d0(i_f, b[i_f]);
  ys_d3[n[i_div]] = fn_y_d3(i_f, b[i_f]);
}

void calc_ode_trapezoidal(int i_f, int i_div, int i_itr) {
  res[0][0] = ys[0];
  REP(i,n[i_div]) {
    //initial guess(midpoint method)
    res[0][i+1] = i ? res[0][i-1] + f(i_f, xs[i], res[0][i])*h*2 : res[0][0];
    //iterations
    REP(j,i_itr+1) {
      res[0][i+1] = f(i_f, xs[i], res[0][i]) + f(i_f, xs[i+1], res[0][i+1]);
      res[0][i+1] += res[0][i+1]*h*0.5 + res[0][i];
      res[1][i+1] = ys[i+1] - res[0][i+1];
    }
    res[2][i+1] = max(abs(ys_d3[i+1]), abs(ys_d3[i]));
    res[2][i+1] *= (exp((xs[i+1]-xs[0])) - 1)/lipschitz_c[i_f];
    res[2][i+1] *= h*h/12;
  }
}

void output(int i_f, int i_div) {
  puts("func,h,x,y_h,y,|y - y_h|,bound");
  REP(i,n[i_div]+1) 
    printf("%d,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", 
           i_f, h, xs[i], res[0][i], ys[i], res[1][i], res[2][i]);
}

int main() {
  REP(i,NF) REP(j,M) REP(k,I) {
if (i!=0 || j!=0 || k!=0) continue;
    init(i,j);
    calc_ode_trapezoidal(i,j,k);
    output(i,j);
  }
  return 0;
}

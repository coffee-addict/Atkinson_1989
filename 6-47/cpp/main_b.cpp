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
const int M = 4;  //# of kinds of step sizes
const int N = 160; //(10-0)/.0625
int n[M] = {20, 40, 80, 160};
double h;
double a[NF] = {0};
double b[NF] = {10};
double y_in[NF] = {0};
double xs[N+1];
double res[2][N+1]; //y_h, error estimate

//y' = f(x,y)
double fn_y_d1(int i_f, double x, double y) {
  double res = 0;
  switch(i_f) {
    case 0:
      res = x-y*y;
      break;
    default:
      res = -9999999999;
      break;
  }
  return res;
}

double fn_y_d2(int i_f, double x, double y) {
  double res = 0;
  switch(i_f) {
    case 0:
      res = 1.0 - y*2*(x - y*y);
      break;
    default:
      res = -9999999999;
      break;
  }
  return res;
}

//y^(3)
double fn_y_d3(int i_f, double x, double y) {
  double res = 0;
  switch(i_f) {
    case 0:
      res = (y + (x - y*y)*(x - y*y*3))*(-2.0);
      break;
    default:
      res = -9999999999;
      break;
  }
  return res;
}

//y^(4)
double fn_y_d4(int i_f, double x, double y) {
  double res = 0;
  switch(i_f) {
    case 0:
      res = x*3 - y*y*6;
      res += y*4*(x-y*y)*(y*y*3 - 2);
      res *= -2.0;
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
  memset(res, 0, sizeof(res));
  REP(k,n[i_div]) xs[k] = a[i_f] + h*k;
  xs[n[i_div]] = b[i_f];
}

void calc_ode_runge_kutta(int i_f, int i_div) {
  res[0][0] = y_in[i_f];
  REP(i,n[i_div]) {
    double v1 = fn_y_d1(i_f, xs[i], res[0][i]);
    double v2 = fn_y_d1(i_f, xs[i]+h*0.5, res[0][i]+v1*h*0.5);
    double v3 = fn_y_d1(i_f, xs[i]+h*0.5, res[0][i]+v2*h*0.5);
    double v4 = fn_y_d1(i_f, xs[i]+h, res[0][i]+v3*h);
    double v = (v1+v2*2+v3*2+v4);
    res[0][i+1] = res[0][i] + v*h/6;
  }
}

void output(int i_f, int i_div) {
  puts("func,h,x,y_h,estimate");
  REP(i,n[i_div]+1) 
    printf("%d,%.5e,%.5e,%.5e\n", i_f, h, xs[i], res[0][i]);
}

int main() {
  REP(i,NF) REP(j,M) {
    init(i,j);
    calc_ode_runge_kutta(i,j);
    output(i,j);
  }
  return 0;
}

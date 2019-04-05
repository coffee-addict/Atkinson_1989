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
const int NF = 5;
const int M = 9;
const int N = 1<<M;
double a[NF] = {0, 0, -4, 0, 0};
double b[NF] = {1, 1, 4, acos(-1.0)*2, acos(-1.0)}; 
double xs[N+1], ys[N+1];
double res[2][M+1];

double fa(double x) {
  return exp(-x*x);
}

double fb(double x) {
  return pow(x,2.5);
}

double fc(double x) {
  return 1.0 / (1.0 + x*x);
}

double fd(double x) {
  return 1.0 / (2.0 + cos(x));
}

double fe(double x) {
  return exp(x) * cos(x*4);
}

double (*fn_p[])(double) = {fa, fb, fc, fd, fe};

void init(int i) {
  double h = (b[i] - a[i])/N;
  memset(res, 0, sizeof(res));
  REP(j,N+1) {
    xs[j] = a[i] + h*j;
    ys[j] = fn_p[i](xs[j]);
  }
}

int main() {
  REP(i,NF) {
    init(i);
    res[0][0] = (ys[0] + ys[N])*0.5;
    int m = 1, n = 1;
    while (m <= M) {
      res[0][m] = res[0][m-1];
      REP(j,n) res[0][m] += ys[N*(j*2+1)/(n*2)];
      m++;
      n <<= 1;
    }
    REP(j,M+1) res[0][j] *= (b[i]-a[i])/(1<<j);
    printf("case (%c):\n", 'a'+i);
    puts("n e r");
    FOR(j,2,M) res[1][j] = (res[0][j] - res[0][j-1])/(res[0][j+1] - res[0][j]);
    FOR(j,1,M+1) printf("%d %.10f %.10f\n", 1<<j, res[0][j], res[1][j]);
  }
  return 0;
}

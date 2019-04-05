#include <bits/stdc++.h>

using namespace std;

#define FOR(i,s,e) for(int (i)=(s);(i)<(int)(e);(i)++)
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

double f(double x) {
  if (x >= 0.0) 
    return sqrt(x);
  else
    return -sqrt(-x);
}

double ff(double x) {
  double res = 0.5;
  if (x >= 0.0) 
    res /= sqrt(x);
  else {
    res *= -1.0;
    res /= sqrt(-x);
  }
  return res;
}

int main() {
  double x = 0.5;
  REP(i,10) {
    printf("(x, f(x)) = %f, %f\n", x, f(x));
    x = x - f(x)/ff(x);
  }
  return 0;
}

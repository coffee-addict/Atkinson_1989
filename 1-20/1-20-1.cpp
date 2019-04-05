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

int main() {
  int d = 100;
  double x0, y0;
  scanf("%lf%lf ", &x0, &y0);
  printf("x = %e, y = %e\n", x0, y0);
  puts("p, res");
  printf("%d,%e\n", 1, x0+y0);
  double x = x0;
  double y = y0;
  FOR(p,2,d+1) {
    x *= x0;
    y *= y0;
    double z = abs(x0);
    z *= pow(pow(y/x, (double)p) + 1.0, 1.0/(double)p);
    printf("%d,%e\n", p, z);
  }
  return 0;
}

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
  double x = 1.6, y = 1.2;
  REP(i,5) {
    printf("i = %d: (x, y) = (%.6f %.6f)\n", i, x, y);
    x = (x*2 + 5.0/x)/4.0;
    y = (y*2 + 3.0/y)/4.0;
  }
  return 0;
}

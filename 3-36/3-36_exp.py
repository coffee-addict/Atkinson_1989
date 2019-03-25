#!/usr/local/bin/python3
import math
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# the number of points to calculate sup-norm
k = 10**3

### define the function and derivatives
s = 0
e = 1
x = np.arange(s, e, (e-s)/k)
spline_fn_type = 'exp'
y0 = np.exp(x)
y1 = np.exp(x)
y2 = np.exp(x)
y3 = np.exp(x)

with open('x.dat', mode='w', encoding='utf-8') as f_x:
    for j, e in enumerate(x):
        f_x.write('{0}\n'.format(x[j]))

### set bc_type
spline_bc_type = 'a-3.7.20'
#spline_bc_type = 'b-not-a-knot'
#spline_bc_type = 'c-natural'
m = 8 #n <= 2^(m-1)
for i in range(3,m):
    n = 2**i
    xs = np.arange(s, e, e/n)
    y_sample = np.sin(xs)
### BCs
#    #for a
#    df_s = np.power(s,0.5)*1.5
#    df_e = np.power(e,0.5)*1.5
#    cs = CubicSpline(xs, y_sample, bc_type=((1,df_s), (1,df_e)))
    #for knot-a-knot bc
    cs = CubicSpline(xs, y_sample)
#    #for other c
#    cs = CubicSpline(xs, y_sample, bc_type='natural')
    e0 = cs(x) - y0
    e1 = cs(x, 1) - y1
    e2 = cs(x, 2) - y2
    e3 = cs(x, 3) - y3
    with open(spline_fn_type + '_' + spline_bc_type + '_' + f'{n:03}' + '.dat', mode='w', encoding='utf-8') as f:
        f.write(f'{n:03}-e0,{n:03}-e1,{n:03}-e2,{n:03}-e3\n')
        for j, e in enumerate(x):
            f.write('{0},{1},{2},{3}\n'.format(e0[j], e1[j], e2[j], e3[j]))

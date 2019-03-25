import math
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

PI = math.pi

### define the function and derivatives
s = 0
e = PI*0.5
k = 10**3
x = np.arange(s, e, e/k)
spline_fn_type = 'sin'
y0 = np.sin(x)
y1 = np.cos(x)
y2 = -np.sin(x)
y3 = -np.cos(x)

### set bc_type
spline_bc_type = 'a-3.7.20'
#spline_bc_type = 'not-a-knot'
#spline_bc_type = 'natural'
m = 8 #n <= 2^(m-1)
for i in range(3,m):
    n = 2**i
    xs = np.arange(s, e, e/n)
    y_sample = np.sin(xs)
    ### BCs
    #for a
    df_s = np.cos(s)
    df_e = np.cos(e)
    cs = CubicSpline(xs, y_sample, bc_type=((1,df_s), (1,df_e)))
    #for knot-a-knot bc
#    cs = CubicSpline(xs, y_sample)
    #for other c
#    cs = CubicSpline(xs, y_sample, bc_type=spline_bc_type)
    e0 = cs(x) - y0
    e1 = cs(x, 1) - y1
    e2 = cs(x, 2) - y2
    e3 = cs(x, 3) - y3
    with open(spline_fn_type + '_' + spline_bc_type + '_' + f'{n:03}' + '.dat', mode='w', encoding='utf-8') as f:
        f.write(f'x,{n:03}-e0,{n:03}-e1,{n:03}-e2,{n:03}-e3\n')
        for j, e in enumerate(x):
            f.write('{0},{1},{2},{3},{4}\n'.format(x[j], e0[j], e1[j], e2[j], e3[j]))

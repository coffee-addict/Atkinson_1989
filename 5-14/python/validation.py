#!/usr/local/bin/python3

import math
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

PI = math.pi

def fa(xs):
    return np.exp(-xs*xs)

def fb(xs):
    return np.float_power(xs,2.5)

def fc(xs):
    return 1.0 / (1.0 + xs*xs)

def fd(xs):
    return 1.0 / (2.0 + np.cos(xs))

def fe(xs):
    return np.exp(xs) * np.cos(xs*4)

NF = 5
M = 9
N = 1<<M
a = np.array([0, 0, -4, 0, 0])
b = np.array([1, 1, 4, PI*2, PI])
fn = np.array([fa, fb, fc, fd, fe])

def calc_gauss_legendre(i, n):
    xs, ws = np.polynomial.legendre.leggauss(n+1)
    xs = (xs + 1)*(b[i] - a[i])*0.5 + a[i]
    ys = np.array(fn[i](xs))
    return sum(ws * ys)*(b[i] - a[i])*0.5

for i in range(NF):
    res = np.zeros((2,M+1))
    m = 1
    n = 1<<m
    while m <= M:
        res[0][m] = calc_gauss_legendre(i,n)
        m += 1
        n *= 2
    print('case ({0}):'.format(chr(ord('a')+i)))
    print('n e r')
    for j in np.arange(2,M):
        res[1][j] = (res[0][j] - res[0][j-1]) / (res[0][j+1] - res[0][j])
    for j in np.arange(1,M+1):
        print('{0} {1} {2}'.format(int(1<<j), res[0][j], res[1][j]))

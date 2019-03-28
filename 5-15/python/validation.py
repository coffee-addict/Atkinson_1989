#!/usr/local/bin/python3
import math
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

PI = math.pi

def f(xs):
    return 1.0 / (1.0 + xs*xs)

N = int(10)
a = -4.0
b = 4.0

def calc_gauss_legendre(n):
    xs, ws = np.polynomial.legendre.leggauss(n+1)
    xs = (xs + 1)*(b - a)*0.5 + a
    ys = np.array(f(xs))
    return sum(ws*ys)*(b-a)*0.5

n = 1
print('n e')
while n*2 <= N:
    print('{0} {1}'.format(n*2, calc_gauss_legendre(n*2)))
    n += 1

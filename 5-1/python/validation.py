#!/usr/local/bin/python3

import numpy as np
import scipy.integrate as integrate

def f(x):
    return np.exp(-x*x)

a, b = 0, 1
n = 2**9
xs = np.arange(a, b, (b-a)/n)
ys = np.k

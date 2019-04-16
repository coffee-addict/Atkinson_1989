#!/usr/local/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt



#def f(t):
#    return np.exp(-t) * np.cos(2.*np.pi*t)

def main():
    in_data = np.genfromtxt('result.dat', delimiter=',')
    print(in_data)
#    print(in_data.size)
#    print(in_data.shape)
#    print(in_data[303][3])

#    x = in_data[:,2]
#    y = in_data[:,3]

#    t0 = np.arange(in_data, 5., 0.1)
#    t1 = np.arange(0., 5., 0.2)
#
#    plt.figure(1)
#    plt.subplot(211)
#    plt.plot(x, y, "bo")
#    plt.plot(t1, f(t1), "k")
#    plt.subplot(212)
#    plt.plot(t1, np.cos(2.*np.pi*t1), "r--")
    plt.show()

if __name__ == '__main__':
    main()



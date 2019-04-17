#!/usr/local/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt



#def f(t):
#    return np.exp(-t) * np.cos(2.*np.pi*t)

def main():
    fig, axs = plt.subplots(3,1)
    for i in range(3):
        x = []
        yt = []
        yh = []
        for j in range(3):
            in_file = 'results_' + str(i) + str(j) + '.dat'
            in_data = np.genfromtxt(in_file, delimiter=',')
            x.append(in_data[:,2])
            yh.append(in_data[:,3])
            if j==2:
                yt= in_data[:,4]
        axs[i].plot(x[2],yt, ls='-')
#        axs[i].plot(x[0],yh[0], ls='--', marker=5)
        axs[i].plot(x[1],yh[1], ls='-.', marker='.')
        axs[i].plot(x[2],yh[2], ls=':')
        axs[i].set_xlabel('x')
        axs[i].set_ylabel('y')
#        axs[i].legend(('true values', 'h = .5', 'h = .25', 'h = .125'))
        axs[i].legend(('true values', 'h = .25', '.125'))
    axs[0].set_title('lambda = -1')
    axs[1].set_title('lambda = -10')
    axs[2].set_title('lambda = -50')
#        line_yt = plt.plot(x[2], yt)
#        line_yh0 = plt.plot(x[0], yh[0])
#        line_yh1 = plt.plot(x[1], yh[1])
#        line_yh2 = plt.plot(x[2], yh[2])
#        plt.setp(line_yt, c='black', lw='1.0')
#        plt.setp(line_yh0, c='blue', lw='1.0', ls='-.', marker='.')
#        plt.setp(line_yh1, c='r', lw='1.0', ls='--', marker='+')
#        plt.setp(line_yh2, c='g', lw='1.0', ls=':', marker='')
#        plt.xlabel('x')
#        plt.ylabel('y')
#        plt.legend(('true values', 'lambda = -1', 'lambda = -10', 'lambda = -50'))
    plt.title('The graphs of the numerical results in problem 6-50', 
              y=-0.5)
    plt.show()

if __name__ == '__main__':
    main()



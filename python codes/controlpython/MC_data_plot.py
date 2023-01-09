#import slycot
import matplotlib.pyplot as plt
import numpy as np
N = 4000
t = np.arange(0,N/2000,1/2000) # t = 0:1/2000:500/2000
output = np.loadtxt('pi_stepout.dat')
input = np.loadtxt('pi_stepin.dat')
plt.figure()
plt.plot(t,output)
plt.title("Output")
plt.xlabel("n")
plt.ylabel("y[n]")
plt.figure()
plt.plot(t,input)
plt.title("Input")
plt.xlabel("n")
plt.ylabel("x[n]")
plt.show()
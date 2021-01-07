import matplotlib.pyplot as plt
import numpy as np
import csv
import math
from scipy import optimize

tstart = 10
tfinal = 99
Nt  = 90
Nx = 1000
DIR='_output'
file_num = np.linspace(tstart,tfinal,Nt,dtype='int')
sum_energy=np.zeros(Nt)
sum_pressure2=np.zeros(Nt)
for i in range(Nt):
    file_name = DIR+"/fort.q00%d"%file_num[i]
    pressure,velocity = np.loadtxt(file_name, skiprows=6, usecols=(0,1), unpack=True)
    energy = velocity*velocity
    pressure2 = pressure*pressure
    sum_energy[i] = np.log(np.sum(energy)/Nx)
    sum_pressure2[i] = np.sum(pressure2)/Nx


def test_func(x, a, b):
    return a + b * x

params, params_covariance = optimize.curve_fit(test_func, file_num, sum_energy,
                                               p0=[2, 2])

print(params)
slope = str(params[1])
fo = open("Slope.dat","at")
fo.write(slope+ "\n")
#np.savetxt(fo,slope)
fo.close()
plt.figure(figsize=(6, 4))
plt.xlabel('Time(t)', fontsize=10)
plt.ylabel('log($b_y^2$)', fontsize=10)
plt.plot(file_num,sum_energy, label='Wave energy')
#plt.plot(file_num,sum_pressure2, label='Pressure')
plt.plot(file_num, test_func(file_num, params[0], params[1]),
         label='Exponential fit')

plt.legend(loc='best')

plt.show()
#plt.plot(file_num,sum_pressure2)
#plt.show()

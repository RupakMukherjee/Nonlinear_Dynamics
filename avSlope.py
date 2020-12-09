import math
import numpy as np

sigma=0.5
epsilon=0.05
k=5

file_name = "Slope.dat"
slope = np.loadtxt(file_name, unpack=True)

avSlope = np.average(slope)

fo = open("SlopePlot.dat","at")
#fo.write(sigma, epsilon, k, avSlope+ "\n")
fo.write('%f %f %d %g\n' % (sigma, epsilon, k, avSlope))
#np.savetxt(fo,slope)
fo.close()

from evtk.hl import gridToVTK
import numpy as np

file_name = "SlopePlot"
sigma, epsilon, k, avSlope = np.loadtxt(file_name+".dat", unpack=True)

sigma0=0.1
sigmaN=0.5
epsilon0=0.01
epsilonN=0.05
k0=1
kN=4
Nx=5
Ny=5
Nz=5

z_slice = 4

SIG = np.linspace(sigma0, sigmaN, Nx, dtype='float32')
EPSA = np.linspace(epsilon0, epsilonN, Ny, dtype='float32')
KAPPA = np.linspace(k0, kN, Nz, dtype='float32')

surface = avSlope.reshape((Nx,Ny,Nz))

surfaceFlag = np.require(surface, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])

gridToVTK("./"+file_name, SIG, EPSA, KAPPA, pointData = {file_name : surfaceFlag})

x_vals = SIG
y_vals = EPSA

Data = surface[:,:,z_slice]

x_axis = []
y_axis = []

for i in range (len(x_vals)):
    for j in range (len(y_vals)):
        x_axis=np.append(x_axis,x_vals[i])
        y_axis=np.append(y_axis,y_vals[j])

Data1D = Data.reshape(-1)

add = np.stack((x_axis, y_axis), axis=1)

FData = np.column_stack((add,Data1D))

np.savetxt("ParameterPlot.dat", FData, fmt='%.8e', delimiter='\t', newline='\n', header='', footer='', comments='# ', encoding=None)
print('Finished writing data file ParameterPlot.dat for GNUPlot too!')

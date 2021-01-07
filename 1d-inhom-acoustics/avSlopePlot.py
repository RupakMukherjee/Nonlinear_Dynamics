from evtk.hl import gridToVTK
import numpy as np

file_name = "SlopePlot"
sigma, epsilon, k, avSlope = np.loadtxt(file_name+".dat", unpack=True)

sigma0=0.1
sigmaN=1.0
epsilon0=0.01
epsilonN=0.05
k0=1
kN=10
Nx=10
Ny=5
Nz=10

x_slice = 4
y_slice = 4
z_slice = 6

SIG = np.linspace(sigma0, sigmaN, Nx, dtype='float32')
EPSA = np.linspace(epsilon0, epsilonN, Ny, dtype='float32')
KAPPA = np.linspace(k0, kN, Nz, dtype='float32')

surface = avSlope.reshape((Nx,Ny,Nz))

surfaceFlag = np.require(surface, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])

gridToVTK("./"+file_name, SIG, EPSA, KAPPA, pointData = {file_name : surfaceFlag})

x_vals = SIG
y_vals = EPSA
z_vals = KAPPA

x_axis = []
y_axis = []
z_axis = []

# for i in range (len(x_vals)):
#     for j in range (len(y_vals)):
#         x_axis=np.append(x_axis,x_vals[i])
#         y_axis=np.append(y_axis,y_vals[j])

for i in range (len(x_vals)):
    x_axis=np.append(x_axis,x_vals[i])

for j in range (len(y_vals)):
    y_axis=np.append(y_axis,y_vals[j])

for k in range (len(z_vals)):
    z_axis=np.append(z_axis,z_vals[k])

#Data = surface[:,:,z_slice]
Data = surface[x_slice,:,z_slice]

Data1D = Data.reshape(-1)

#add = np.stack((x_axis, y_axis), axis=1)
add = np.stack((y_axis), axis=0)

FData = np.column_stack((add,Data1D))

np.savetxt("ParameterPlot.dat", FData, fmt='%.8e', delimiter='\t', newline='\n', header='', footer='', comments='# ', encoding=None)
print('Finished writing data file ParameterPlot.dat for GNUPlot too!')

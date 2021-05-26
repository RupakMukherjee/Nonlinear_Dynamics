#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import

def step_Euler_radial(solver,state,dt):
    import numpy as np
    Nx      = 1024
    xmin    = 0
    xmax    = 10*np.pi
    h = (xmax-xmin)/Nx

    nu = 1e-6

    q = state.q

    qstar = np.empty(q.shape)

    for i in range(1,Nx-2):
        qstar[0,i] = (q[0,i+1] - 2*q[0,i] + q[0,i-1])/(h*h)

    qstar[0,0]  = (q[0,1] - 2*q[0,0] + q[0,Nx-1])/(h*h)
    qstar[0,Nx-1] = (q[0,0] - 2*q[0,Nx-1] + q[0,Nx-2])/(h*h)

    q[0,:] = q[0,:] + nu * qstar[0,:]
    q[1,:] = q[1,:]

def acoustics(problem='Brahmananda'):
    """
    This example solves the 1-dimensional variable-coefficient acoustics
    equations in a medium with a single interface.
    """
    from numpy import sqrt, abs
    from clawpack import pyclaw
    from clawpack import riemann
    import numpy as np
    import math

    solver = pyclaw.ClawSolver1D(riemann.acoustics_variable_1D)
    solver.step_source = step_Euler_radial
    solver.source_split = 1

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic

    tFinal  = 100
    tFrames = 100
    Nx      = 1024
    xmin    = 0
    xmax    = 10*math.pi

    sigma=0.3
    mu=0.0
    epsilon=0.03
    k=6

    x = pyclaw.Dimension(xmin,xmax,Nx,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 2
    num_aux = 2
    state = pyclaw.State(domain, num_eqn, num_aux)

    if problem == 'Brahmananda':
        rho  = 1.0
        bulk = epsilon

    xc = domain.grid.x.centers

    for i in range(len(xc)):
        U1             = np.random.rand()
        U2             = np.random.rand()
        GRn            = sigma * sqrt(-2*math.log(U1))*math.cos(2*math.pi*U2) + mu
        z              = 1 + bulk*GRn
        state.aux[0,i] = rho*z   # Impedance
        state.aux[1,i] = z       # Sound speed

    state.q[0,:] = np.sin(k*xc)
    state.q[1,:] = state.q[0,:] + 0.

    claw                  = pyclaw.Controller()
    claw.solution         = pyclaw.Solution(state, domain)
    claw.solver           = solver
    claw.tfinal           = tFinal #100#int(131*lx)
    claw.num_output_times = tFrames #100

    # Solve
    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)

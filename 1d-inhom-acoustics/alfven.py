#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
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
    import press_ran as pran

    solver = pyclaw.ClawSolver1D(riemann.acoustics_variable_1D)

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic

    tFinal  = 100
    tFrames = 100
    Nx      = 1024
    xmin    = 0
    xmax    = 10*math.pi

    sigma=sigma
    mu=mu
    epsilon=epsilon
    k=k

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
    claw.tfinal           = tFinal
    claw.num_output_times = tFrames

    # Solve
    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)

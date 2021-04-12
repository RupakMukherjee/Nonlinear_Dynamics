program magnetosonic

implicit none

include "fftw3.f"

! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 64
integer ( kind = 4 ), parameter :: Nz = 64
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real ( kind = 8 ) x(Nx),z(Nz)
real ( kind = 8 ) ux(Nx,Nz),uz(Nx,Nz),bx(Nx,Nz),bz(Nx,Nz)
real ( kind = 8 ) dux_dt(Nx,Nz),duz_dt(Nx,Nz),dbx_dt(Nx,Nz),dbz_dt(Nx,Nz)

integer ( kind = 8 ) :: planf, planb

integer ( kind = 4 ) i,k,t
REAL ( kind = 8 ) :: time, time_min, time_max, dt, f, rho0, B0, Re, Rm
REAL ( kind = 8 ) :: dx, dz, Lx, Lz, kx, kz, uEnergy, bEnergy, Energy

COMMON/comm/Lx, Lz, dx, dz, dt, rho0, B0, Re, Rm

open(unit=10,file='Energy.dat',status='unknown')
open(unit=15,file='Initial_Grid_Data.dat',status='unknown')

!===================== USER INPUTS ============================================

Lx = 2.0d0*pi
Lz = 2.0d0*pi

dx = Lx/dfloat(Nx)
dz = Lz/dfloat(Nz)

time_min = 0.00d0
time_max = 10.00d0
dt = 0.010d0

f = 4.0d0

Re = 1e3
Rm = 1e3

rho0 = 1.0d0
B0 = 1.0d0

! Grid Generation.
do i = 1, Nx
  x(i) = dfloat(i-1) * dx
  do k = 1, Nz
    z(k) = dfloat(k-1) * dz
    ux(i,k) = dsin(f*x(i))
    uz(i,k) = dcos(f*z(k))
    bx(i,k) = +dcos(f*x(i))*dsin(f*z(k))
    bz(i,k) = -dsin(f*x(i))*dcos(f*z(k))
    write(15,*) 0.0d0, x(i),z(k),ux(i,k),uz(i,k),bx(i,k),bz(i,k)
  enddo
enddo

!======================= MAIN PROGRAM =====================================================

DO time = time_min, time_max, dt

Energy = 0.0d0
uEnergy = 0.0d0
bEnergy = 0.0d0

t = nint(time/dt) - dint(time_min/dt)

!Calculating the time evolution in real-space...

CALL derive(Nx, Nz, Nh, pi, time, ux, uz, bx, bz, dux_dt, duz_dt, dbx_dt, dbz_dt)
CALL rk4(Nx, Nz, Nh, pi, time, ux, uz, bx, bz, dux_dt, duz_dt, dbx_dt, dbz_dt)

!PRINT*, "Printing Data"

! IF ( t /= 0 .and. MOD(t,10) == 0 ) THEN
!   DO i = 1, Nx
!     do k = 1, Nz
!       WRITE (t+100,*) t, x(i), z(k), ux(i,k), uz(i,k), bx(i,k), bz(i,k)
!     enddo
!  ENDDO
!   close(t+100)
! ENDIF

!IF ( MOD(t,100) == 0 ) THEN
  DO i = 1, Nx
    do k = 1, Nz
      uEnergy = uEnergy + (ux(i,k) * ux(i,k) + uz(i,k) * uz(i,k)) / dfloat(Nx*Nz)
      bEnergy = bEnergy + (bx(i,k) * bx(i,k) + bz(i,k) * bz(i,k)) / dfloat(Nx*Nz)
      Energy = uEnergy + bEnergy
    enddo
  ENDDO
  WRITE(10,*) time, uEnergy, bEnergy, Energy!log(bEnergy)!, freq(t)
!ENDIF

ENDDO ! time


contains
!===================================================================================

SUBROUTINE derive(Nx, Nz, Nh, pi, time, ux, uz, bx, bz, dux_dt, duz_dt, dbx_dt, dbz_dt)

implicit none

INTEGER ( kind = 4 ) Nx, Nz, Nh, i, k
REAL ( kind = 8 ) Lx, Lz, dx, dz, dt, rho0, B0, Re, Rm
REAL ( kind = 8 ) pi, time, kx, kz
REAL ( kind = 8 ) ux(Nx,Nz), ux_dum(Nx,Nz), dux_dt(Nx,Nz)
REAL ( kind = 8 ) uz(Nx,Nz), uz_dum(Nx,Nz), duz_dt(Nx,Nz)
REAL ( kind = 8 ) bx(Nx,Nz), bx_dum(Nx,Nz), dbx_dt(Nx,Nz)
REAL ( kind = 8 ) bz(Nx,Nz), bz_dum(Nx,Nz), dbz_dt(Nx,Nz)
REAL ( kind = 8 ) dux_dz(Nx,Nz), duz_dx(Nx,Nz), dbx_dz(Nx,Nz), dbz_dx(Nx,Nz)
REAL ( kind = 8 ) d2ux_dx2(Nx,Nz), d2ux_dz2(Nx,Nz), d2uz_dx2(Nx,Nz), d2uz_dz2(Nx,Nz)
REAL ( kind = 8 ) d2bx_dx2(Nx,Nz), d2bx_dz2(Nx,Nz), d2bz_dx2(Nx,Nz), d2bz_dz2(Nx,Nz)
COMPLEX ( kind = 8 ) uxk(Nh,Nz), uzk(Nh,Nz), bxk(Nh,Nz), bzk(Nh,Nz)
COMPLEX ( kind = 8 ) ikz_uxk(Nh,Nz), ikx_uzk(Nh,Nz), ikz_bxk(Nh,Nz), ikx_bzk(Nh,Nz)
COMPLEX ( kind = 8 ) kx2_uxk(Nh,Nz), kz2_uxk(Nh,Nz), kx2_uzk(Nh,Nz), kz2_uzk(Nh,Nz)
COMPLEX ( kind = 8 ) kx2_bxk(Nh,Nz), kz2_bxk(Nh,Nz), kx2_bzk(Nh,Nz), kz2_bzk(Nh,Nz)

common/comm/Lx, Lz, dx, dz, dt, rho0, B0, Re, Rm

do i = 1, Nx
  do k = 1, Nz
    ux_dum(i,k) = ux(i,k)
    uz_dum(i,k) = uz(i,k)
    bx_dum(i,k) = bx(i,k)
    bz_dum(i,k) = bz(i,k)
  enddo
enddo

call dfftw_plan_dft_r2c_2d_(planf, Nx, Nz, ux_dum, uxk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_2d_(planf, Nx, Nz, uz_dum, uzk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_2d_(planf, Nx, Nz, bx_dum, bxk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_2d_(planf, Nx, Nz, bz_dum, bzk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

do i = 1, Nh
  kx = 2.0d0 * pi * dfloat(i-1) / Lx
  do k = 1, Nz/2
    kz = 2.0d0 * pi * dfloat(k-1) / Lz
    ikz_uxk(i,k) = (0.0d0,1.0d0) * kz * uxk(i,k)
    ikx_uzk(i,k) = (0.0d0,1.0d0) * kx * uzk(i,k)
    ikz_bxk(i,k) = (0.0d0,1.0d0) * kz * bxk(i,k)
    ikx_bzk(i,k) = (0.0d0,1.0d0) * kx * bzk(i,k)
    kx2_uxk(i,k) = - kx * kx * uxk(i,k)
    kz2_uxk(i,k) = - kz * kz * uxk(i,k)
    kx2_uzk(i,k) = - kx * kx * uzk(i,k)
    kz2_uzk(i,k) = - kz * kz * uzk(i,k)
    kx2_bxk(i,k) = - kx * kx * bxk(i,k)
    kz2_bxk(i,k) = - kz * kz * bxk(i,k)
    kx2_bzk(i,k) = - kx * kx * bzk(i,k)
    kz2_bzk(i,k) = - kz * kz * bzk(i,k)
  enddo
  do k = Nz/2+1, Nz
    kz = 2.0d0 * pi * dfloat((k-1) - Nz) / Lz
    ikz_uxk(i,k) = (0.0d0,1.0d0) * kz * uxk(i,k)
    ikx_uzk(i,k) = (0.0d0,1.0d0) * kx * uzk(i,k)
    ikz_bxk(i,k) = (0.0d0,1.0d0) * kz * bxk(i,k)
    ikx_bzk(i,k) = (0.0d0,1.0d0) * kx * bzk(i,k)
    kx2_uxk(i,k) = - kx * kx * uxk(i,k)
    kz2_uxk(i,k) = - kz * kz * uxk(i,k)
    kx2_uzk(i,k) = - kx * kx * uzk(i,k)
    kz2_uzk(i,k) = - kz * kz * uzk(i,k)
    kx2_bxk(i,k) = - kx * kx * bxk(i,k)
    kz2_bxk(i,k) = - kz * kz * bxk(i,k)
    kx2_bzk(i,k) = - kx * kx * bzk(i,k)
    kz2_bzk(i,k) = - kz * kz * bzk(i,k)
  enddo
enddo

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, ikz_uxk, dux_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, ikx_uzk, duz_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, ikz_bxk, dbx_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, ikx_bzk, dbz_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kx2_uxk, d2ux_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kz2_uxk, d2ux_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kx2_uzk, d2uz_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kz2_uzk, d2uz_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kx2_bxk, d2bx_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kz2_bxk, d2bx_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kx2_bzk, d2bz_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_2d_(planb, Nx, Nz, kz2_bzk, d2bz_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

do i = 1, Nx
  do k = 1, Nz
    dux_dz(i,k)   = dux_dz(i,k)   / dfloat(Nx*Nz)
    duz_dx(i,k)   = duz_dx(i,k)   / dfloat(Nx*Nz)
    dbx_dz(i,k)   = dbx_dz(i,k)   / dfloat(Nx*Nz)
    dbz_dx(i,k)   = dbz_dx(i,k)   / dfloat(Nx*Nz)
    d2ux_dx2(i,k) = d2ux_dx2(i,k) / dfloat(Nx*Nz)
    d2ux_dz2(i,k) = d2ux_dz2(i,k) / dfloat(Nx*Nz)
    d2uz_dx2(i,k) = d2uz_dx2(i,k) / dfloat(Nx*Nz)
    d2uz_dz2(i,k) = d2uz_dz2(i,k) / dfloat(Nx*Nz)
    d2bx_dx2(i,k) = d2bx_dx2(i,k) / dfloat(Nx*Nz)
    d2bx_dz2(i,k) = d2bx_dz2(i,k) / dfloat(Nx*Nz)
    d2bz_dx2(i,k) = d2bz_dx2(i,k) / dfloat(Nx*Nz)
    d2bz_dz2(i,k) = d2bz_dz2(i,k) / dfloat(Nx*Nz)
    enddo
enddo

DO i = 1, Nx
  do k = 1, Nz
    ! dux_dt(i,k) = - B0 * dbx_dz(i,k) / rho0 + (d2ux_dx2(i,k) + d2ux_dz2(i,k)) / (rho0 * Re)
    ! duz_dt(i,k) = + B0 * dbz_dx(i,k) / rho0 + (d2uz_dx2(i,k) + d2uz_dz2(i,k)) / (rho0 * Re)
    dux_dt(i,k) = - B0 * (dbx_dz(i,k) - dbz_dx(i,k)) / rho0 + (d2ux_dx2(i,k) + d2ux_dz2(i,k)) / (rho0 * Re)
    duz_dt(i,k) = + (d2uz_dx2(i,k) + d2uz_dz2(i,k)) / (rho0 * Re)
    dbx_dt(i,k) = - B0 * dux_dz(i,k) + (d2bx_dx2(i,k) + d2bx_dz2(i,k)) / Rm
    dbz_dt(i,k) = + B0 * duz_dx(i,k) + (d2bz_dx2(i,k) + d2bz_dz2(i,k)) / Rm
  enddo
enddo

return
end subroutine derive

!====================================================================================

SUBROUTINE rk4(Nx, Nz, Nh, pi, time, ux, uz, bx, bz, dux_dt, duz_dt, dbx_dt, dbz_dt)
implicit none

INTEGER ( kind = 4 ) Nx, Nz, Nh, i, k
REAL ( kind = 8 ) pi, time, Lx, Lz, dx, dz, dt, rho0, B0, Re, Rm
REAL ( kind = 8 ) ux(Nx,Nz),dux_dt(Nx,Nz),k1_ux(Nx,Nz),k2_ux(Nx,Nz),k3_ux(Nx,Nz),k4_ux(Nx,Nz),dum_ux(Nx,Nz)
REAL ( kind = 8 ) uz(Nx,Nz),duz_dt(Nx,Nz),k1_uz(Nx,Nz),k2_uz(Nx,Nz),k3_uz(Nx,Nz),k4_uz(Nx,Nz),dum_uz(Nx,Nz)
REAL ( kind = 8 ) bx(Nx,Nz),dbx_dt(Nx,Nz),k1_bx(Nx,Nz),k2_bx(Nx,Nz),k3_bx(Nx,Nz),k4_bx(Nx,Nz),dum_bx(Nx,Nz)
REAL ( kind = 8 ) bz(Nx,Nz),dbz_dt(Nx,Nz),k1_bz(Nx,Nz),k2_bz(Nx,Nz),k3_bz(Nx,Nz),k4_bz(Nx,Nz),dum_bz(Nx,Nz)

common/comm/Lx, Lz, dx, dz, dt, rho0, B0, Re, Rm

!PRINT*, "RKIdx = 1"

do i = 1, Nx
  do k = 1, Nz
    k1_ux(i,k) = dux_dt(i,k)
    k1_uz(i,k) = duz_dt(i,k)
    dum_ux(i,k) = ux(i,k) + k1_ux(i,k) * dt / 2.0d0
    dum_uz(i,k) = uz(i,k) + k1_uz(i,k) * dt / 2.0d0
  enddo
enddo

do i = 1, Nx
  do k = 1, Nz
    k1_bx(i,k) = dbx_dt(i,k)
    k1_bz(i,k) = dbz_dt(i,k)
    dum_bx(i,k) = bx(i,k) + k1_bx(i,k) * dt / 2.0d0
    dum_bz(i,k) = bz(i,k) + k1_bz(i,k) * dt / 2.0d0
  enddo
enddo

CALL derive(Nx, Nz, Nh, pi, time+dt/2.0d0, dum_ux, dum_uz, dum_bx, dum_bz, dux_dt, duz_dt, dbx_dt, dbz_dt)

!PRINT*, "RKIdx = 2"

do i = 1, Nx
  do k = 1, Nz
    k2_ux(i,k) = dux_dt(i,k)
    k2_uz(i,k) = duz_dt(i,k)
    dum_ux(i,k) = ux(i,k) + k2_ux(i,k) * dt/2.0d0
    dum_uz(i,k) = uz(i,k) + k2_uz(i,k) * dt/2.0d0
  enddo
enddo

do i = 1, Nx
  do k = 1, Nz
    k2_bx(i,k) = dbx_dt(i,k)
    k2_bz(i,k) = dbz_dt(i,k)
    dum_bx(i,k) = bx(i,k) + k2_bx(i,k) * dt/2.0d0
    dum_bz(i,k) = bz(i,k) + k2_bz(i,k) * dt/2.0d0
  enddo
enddo

CALL derive(Nx, Nz, Nh, pi, time+dt/2.0d0, dum_ux, dum_uz, dum_bx, dum_bz, dux_dt, duz_dt, dbx_dt, dbz_dt)

!PRINT*, "RKIdx = 3"

do i = 1, Nx
  do k = 1, Nz
    k3_ux(i,k) = dux_dt(i,k)
    k3_uz(i,k) = duz_dt(i,k)
    dum_ux(i,k) = ux(i,k) + k3_ux(i,k) * dt
    dum_uz(i,k) = uz(i,k) + k3_uz(i,k) * dt
  enddo
enddo

do i = 1, Nx
  do k = 1, Nz
    k3_bx(i,k) = dbx_dt(i,k)
    k3_bz(i,k) = dbz_dt(i,k)
    dum_bx(i,k) = bx(i,k) + k3_bx(i,k) * dt
    dum_bz(i,k) = bz(i,k) + k3_bz(i,k) * dt
  enddo
enddo

CALL derive(Nx, Nz, Nh, pi, time+dt, dum_ux, dum_uz, dum_bx, dum_bz, dux_dt, duz_dt, dbx_dt, dbz_dt)

!PRINT*, "RKIdx = 4"

do i = 1, Nx
  do k = 1, Nz
    k4_ux(i,k) = dux_dt(i,k)
    k4_uz(i,k) = duz_dt(i,k)
    ux(i,k) = ux(i,k) + dt/6.0d0*(k1_ux(i,k) + 2.0d0*k2_ux(i,k) + 2.0d0*k3_ux(i,k) + k4_ux(i,k))
    uz(i,k) = uz(i,k) + dt/6.0d0*(k1_uz(i,k) + 2.0d0*k2_uz(i,k) + 2.0d0*k3_uz(i,k) + k4_uz(i,k))
  enddo
enddo

do i = 1, Nx
  do k = 1, Nz
    k4_bx(i,k) = dbx_dt(i,k)
    k4_bz(i,k) = dbz_dt(i,k)
    bx(i,k) = bx(i,k) + dt/6.0d0*(k1_bx(i,k) + 2.0d0*k2_bx(i,k) + 2.0d0*k3_bx(i,k) + k4_bx(i,k))
    bz(i,k) = bz(i,k) + dt/6.0d0*(k1_bz(i,k) + 2.0d0*k2_bz(i,k) + 2.0d0*k3_bz(i,k) + k4_bz(i,k))
  enddo
enddo

RETURN
end subroutine rk4

!==================================================================================

end program magnetosonic

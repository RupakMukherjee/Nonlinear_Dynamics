program linmhd3d

implicit none

include "fftw3.f"

! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 16
integer ( kind = 4 ), parameter :: Ny = 16
integer ( kind = 4 ), parameter :: Nz = 16
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real ( kind = 8 ) x(Nx),y(Ny),z(Nz)
real ( kind = 8 ) ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) bx(Nx,Ny,Nz),by(Nx,Ny,Nz),bz(Nx,Ny,Nz)
real ( kind = 8 ) dux_dt(Nx,Ny,Nz),duy_dt(Nx,Ny,Nz),duz_dt(Nx,Ny,Nz)
real ( kind = 8 ) dbx_dt(Nx,Ny,Nz),dby_dt(Nx,Ny,Nz),dbz_dt(Nx,Ny,Nz)

integer ( kind = 8 ) :: planf, planb

integer ( kind = 4 ) i,j,k,t
REAL ( kind = 8 ) :: time, time_min, time_max, dt, f, rho0, B0, Re, Rm
REAL ( kind = 8 ) :: dx, dy, dz, Lx, Ly, Lz, kx, ky, kz, uEnergy, bEnergy, Energy

COMMON/comm/Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm

open(unit=10,file='Energy.dat',status='unknown')
open(unit=15,file='Initial_Grid_Data.dat',status='unknown')

!===================== USER INPUTS ============================================

Lx = 2.0d0*pi
Ly = 2.0d0*pi
Lz = 2.0d0*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)
dz = Lz/dfloat(Nz)

time_min = 0.00d0
time_max = 10.00d0
dt = 0.010d0

f = 2.0d0

Re = 1e3
Rm = 1e3

rho0 = 1.0d0
B0 = 1.0d0

! Grid Generation.
do i = 1, Nx
  x(i) = dfloat(i-1) * dx
  do j = 1, Ny
    y(j) = dfloat(j-1) * dy
    do k = 1, Nz
      z(k) = dfloat(k-1) * dz
      ux(i,j,k) = +dcos(f*x(i))*dsin(f*z(k))
      uy(i,j,k) = 0.0d0
      uz(i,j,k) = -dsin(f*x(i))*dcos(f*z(k))
      bx(i,j,k) = +dcos(f*x(i))*dsin(f*z(k))
      by(i,j,k) = 0.0d0
      bz(i,j,k) = -dsin(f*x(i))*dcos(f*z(k))
      write(15,*) 0.0d0, x(i),y(j),z(k),ux(i,j,k),uy(i,j,k),uz(i,j,k),bx(i,j,k),by(i,j,k),bz(i,j,k)
    enddo
  enddo
enddo

!======================= MAIN PROGRAM =====================================================

DO time = time_min, time_max, dt

Energy = 0.0d0
uEnergy = 0.0d0
bEnergy = 0.0d0

t = nint(time/dt) - dint(time_min/dt)

!Calculating the time evolution in real-space...

CALL derive(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, &
            dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)
CALL rk4(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, &
         dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)

!PRINT*, "Printing Data"

! IF ( t /= 0 .and. MOD(t,10) == 0 ) THEN
!   DO i = 1, Nx
!    do j = 1, Ny
!     do k = 1, Nz
!       WRITE (t+100,*) t, x(i), y(j), z(k), ux(i,j,k), uy(i,j,k), uz(i,j,k), bx(i,j,k), by(i,j,k), bz(i,j,k)
!     enddo
!    enddo
!  ENDDO
!   close(t+100)
! ENDIF

!IF ( MOD(t,100) == 0 ) THEN
  DO i = 1, Nx
    do j = 1, Ny
      do k = 1, Nz
        uEnergy = uEnergy + (ux(i,j,k) * ux(i,j,k) + uy(i,j,k) * uy(i,j,k) + uz(i,j,k) * uz(i,j,k)) / dfloat(Nx*Ny*Nz)
        bEnergy = bEnergy + (bx(i,j,k) * bx(i,j,k) + by(i,j,k) * by(i,j,k) + bz(i,j,k) * bz(i,j,k)) / dfloat(Nx*Ny*Nz)
        Energy = uEnergy + bEnergy
      enddo
    enddo
  ENDDO
  WRITE(10,*) time, uEnergy, bEnergy, Energy!log(bEnergy)!, freq(t)
!ENDIF

ENDDO ! time


contains
!===================================================================================

SUBROUTINE derive(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, &
                  dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)

implicit none

INTEGER ( kind = 4 ) Nx, Ny, Nz, Nh, i, j, k
REAL ( kind = 8 ) Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm
REAL ( kind = 8 ) pi, time, kx, ky, kz
REAL ( kind = 8 ) uxRe, uyRe, uzRe, bxRm, byRm, bzRm
REAL ( kind = 8 ) ux(Nx,Ny,Nz), ux_dum(Nx,Ny,Nz), dux_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) uy(Nx,Ny,Nz), uy_dum(Nx,Ny,Nz), duy_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) uz(Nx,Ny,Nz), uz_dum(Nx,Ny,Nz), duz_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) bx(Nx,Ny,Nz), bx_dum(Nx,Ny,Nz), dbx_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) by(Nx,Ny,Nz), by_dum(Nx,Ny,Nz), dby_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) bz(Nx,Ny,Nz), bz_dum(Nx,Ny,Nz), dbz_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) dux_dz(Nx,Ny,Nz), duy_dz(Nx,Ny,Nz), duz_dx(Nx,Ny,Nz)
REAL ( kind = 8 ) dux_dx(Nx,Ny,Nz), duy_dy(Nx,Ny,Nz)
REAL ( kind = 8 ) dbx_dz(Nx,Ny,Nz), dby_dz(Nx,Ny,Nz), dbz_dx(Nx,Ny,Nz)
REAL ( kind = 8 ) dbx_dx(Nx,Ny,Nz), dby_dy(Nx,Ny,Nz), dbz_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) d2ux_dx2(Nx,Ny,Nz), d2ux_dy2(Nx,Ny,Nz), d2ux_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2uy_dx2(Nx,Ny,Nz), d2uy_dy2(Nx,Ny,Nz), d2uy_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2uz_dx2(Nx,Ny,Nz), d2uz_dy2(Nx,Ny,Nz), d2uz_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2bx_dx2(Nx,Ny,Nz), d2bx_dy2(Nx,Ny,Nz), d2bx_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2by_dx2(Nx,Ny,Nz), d2by_dy2(Nx,Ny,Nz), d2by_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2bz_dx2(Nx,Ny,Nz), d2bz_dy2(Nx,Ny,Nz), d2bz_dz2(Nx,Ny,Nz)
COMPLEX ( kind = 8 ) uxk(Nh,Ny,Nz), uyk(Nh,Ny,Nz), uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) bxk(Nh,Ny,Nz), byk(Nh,Ny,Nz), bzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikz_uxk(Nh,Ny,Nz), ikz_uyk(Nh,Ny,Nz), ikx_uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikx_uxk(Nh,Ny,Nz), iky_uyk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikz_bxk(Nh,Ny,Nz), ikz_byk(Nh,Ny,Nz), ikx_bzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikx_bxk(Nh,Ny,Nz), iky_byk(Nh,Ny,Nz), ikz_bzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_uxk(Nh,Ny,Nz), ky2_uxk(Nh,Ny,Nz), kz2_uxk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_uyk(Nh,Ny,Nz), ky2_uyk(Nh,Ny,Nz), kz2_uyk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_uzk(Nh,Ny,Nz), ky2_uzk(Nh,Ny,Nz), kz2_uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_bxk(Nh,Ny,Nz), ky2_bxk(Nh,Ny,Nz), kz2_bxk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_byk(Nh,Ny,Nz), ky2_byk(Nh,Ny,Nz), kz2_byk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_bzk(Nh,Ny,Nz), ky2_bzk(Nh,Ny,Nz), kz2_bzk(Nh,Ny,Nz)

common/comm/Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
      bx_dum(i,j,k) = bx(i,j,k)
      by_dum(i,j,k) = by(i,j,k)
      bz_dum(i,j,k) = bz(i,j,k)
    enddo
  enddo
enddo

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, ux_dum, uxk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, uy_dum, uyk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, uz_dum, uzk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, bx_dum, bxk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, by_dum, byk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, bz_dum, bzk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

do i = 1, Nh
  kx = 2.0d0 * pi * dfloat(i-1) / Lx
  do j = 1, Ny
    ky = 2.0d0 * pi * dfloat(j-1) / Ly
    do k = 1, Nz/2
      kz = 2.0d0 * pi * dfloat(k-1) / Lz
      ikz_uxk(i,j,k) = (0.0d0,1.0d0) * kz * uxk(i,j,k)
      ikz_uyk(i,j,k) = (0.0d0,1.0d0) * kz * uyk(i,j,k)
      ikx_uzk(i,j,k) = (0.0d0,1.0d0) * kx * uzk(i,j,k)
      ikx_uxk(i,j,k) = (0.0d0,1.0d0) * kx * uxk(i,j,k)
      iky_uyk(i,j,k) = (0.0d0,1.0d0) * ky * uyk(i,j,k)
      ikz_bxk(i,j,k) = (0.0d0,1.0d0) * kz * bxk(i,j,k)
      ikz_byk(i,j,k) = (0.0d0,1.0d0) * kz * byk(i,j,k)
      ikx_bzk(i,j,k) = (0.0d0,1.0d0) * kx * bzk(i,j,k)
      ikx_bxk(i,j,k) = (0.0d0,1.0d0) * kx * bxk(i,j,k)
      iky_byk(i,j,k) = (0.0d0,1.0d0) * ky * byk(i,j,k)
      ikz_bzk(i,j,k) = (0.0d0,1.0d0) * kz * bzk(i,j,k)
      kx2_uxk(i,j,k) = - kx * kx * uxk(i,j,k)
      ky2_uxk(i,j,k) = - ky * ky * uxk(i,j,k)
      kz2_uxk(i,j,k) = - kz * kz * uxk(i,j,k)
      kx2_uyk(i,j,k) = - kx * kx * uyk(i,j,k)
      ky2_uyk(i,j,k) = - ky * ky * uyk(i,j,k)
      kz2_uyk(i,j,k) = - kz * kz * uyk(i,j,k)
      kx2_uzk(i,j,k) = - kx * kx * uzk(i,j,k)
      ky2_uzk(i,j,k) = - ky * ky * uzk(i,j,k)
      kz2_uzk(i,j,k) = - kz * kz * uzk(i,j,k)
      kx2_bxk(i,j,k) = - kx * kx * bxk(i,j,k)
      ky2_bxk(i,j,k) = - ky * ky * bxk(i,j,k)
      kz2_bxk(i,j,k) = - kz * kz * bxk(i,j,k)
      kx2_byk(i,j,k) = - kx * kx * byk(i,j,k)
      ky2_byk(i,j,k) = - ky * ky * byk(i,j,k)
      kz2_byk(i,j,k) = - kz * kz * byk(i,j,k)
      kx2_bzk(i,j,k) = - kx * kx * bzk(i,j,k)
      ky2_bzk(i,j,k) = - ky * ky * bzk(i,j,k)
      kz2_bzk(i,j,k) = - kz * kz * bzk(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kz = 2.0d0 * pi * dfloat((k-1) - Nz) / Lz
      ikz_uxk(i,j,k) = (0.0d0,1.0d0) * kz * uxk(i,j,k)
      ikz_uyk(i,j,k) = (0.0d0,1.0d0) * kz * uyk(i,j,k)
      ikx_uzk(i,j,k) = (0.0d0,1.0d0) * kx * uzk(i,j,k)
      ikx_uxk(i,j,k) = (0.0d0,1.0d0) * kx * uxk(i,j,k)
      iky_uyk(i,j,k) = (0.0d0,1.0d0) * ky * uyk(i,j,k)
      ikz_bxk(i,j,k) = (0.0d0,1.0d0) * kz * bxk(i,j,k)
      ikz_byk(i,j,k) = (0.0d0,1.0d0) * kz * byk(i,j,k)
      ikx_bzk(i,j,k) = (0.0d0,1.0d0) * kx * bzk(i,j,k)
      ikx_bxk(i,j,k) = (0.0d0,1.0d0) * kx * bxk(i,j,k)
      iky_byk(i,j,k) = (0.0d0,1.0d0) * ky * byk(i,j,k)
      ikz_bzk(i,j,k) = (0.0d0,1.0d0) * kz * bzk(i,j,k)
      kx2_uxk(i,j,k) = - kx * kx * uxk(i,j,k)
      ky2_uxk(i,j,k) = - ky * ky * uxk(i,j,k)
      kz2_uxk(i,j,k) = - kz * kz * uxk(i,j,k)
      kx2_uyk(i,j,k) = - kx * kx * uyk(i,j,k)
      ky2_uyk(i,j,k) = - ky * ky * uyk(i,j,k)
      kz2_uyk(i,j,k) = - kz * kz * uyk(i,j,k)
      kx2_uzk(i,j,k) = - kx * kx * uzk(i,j,k)
      ky2_uzk(i,j,k) = - ky * ky * uzk(i,j,k)
      kz2_uzk(i,j,k) = - kz * kz * uzk(i,j,k)
      kx2_bxk(i,j,k) = - kx * kx * bxk(i,j,k)
      ky2_bxk(i,j,k) = - ky * ky * bxk(i,j,k)
      kz2_bxk(i,j,k) = - kz * kz * bxk(i,j,k)
      kx2_byk(i,j,k) = - kx * kx * byk(i,j,k)
      ky2_byk(i,j,k) = - ky * ky * byk(i,j,k)
      kz2_byk(i,j,k) = - kz * kz * byk(i,j,k)
      kx2_bzk(i,j,k) = - kx * kx * bzk(i,j,k)
      ky2_bzk(i,j,k) = - ky * ky * bzk(i,j,k)
      kz2_bzk(i,j,k) = - kz * kz * bzk(i,j,k)
    enddo
  enddo
enddo

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_uxk, dux_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_uyk, duy_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_uzk, duz_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_uxk, dux_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_uyk, duy_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_bxk, dbx_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_byk, dby_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_bzk, dbz_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_bxk, dbx_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_byk, dby_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_bzk, dbz_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kx2_uxk, d2ux_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ky2_uxk, d2ux_dy2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kz2_uxk, d2ux_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kx2_uyk, d2uy_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ky2_uyk, d2uy_dy2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kz2_uyk, d2uy_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kx2_uzk, d2uz_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ky2_uzk, d2uz_dy2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kz2_uzk, d2uz_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kx2_bxk, d2bx_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ky2_bxk, d2bx_dy2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kz2_bxk, d2bx_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kx2_byk, d2by_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ky2_byk, d2by_dy2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kz2_byk, d2by_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kx2_bzk, d2bz_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ky2_bzk, d2bz_dy2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, kz2_bzk, d2bz_dz2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      dux_dz(i,j,k)   = dux_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      duy_dz(i,j,k)   = duy_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      duz_dx(i,j,k)   = duz_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      dux_dx(i,j,k)   = dux_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      duy_dy(i,j,k)   = duy_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbx_dz(i,j,k)   = dbx_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      dby_dz(i,j,k)   = dby_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbz_dx(i,j,k)   = dbz_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbx_dx(i,j,k)   = dbx_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      dby_dy(i,j,k)   = dby_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbz_dz(i,j,k)   = dbz_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      d2ux_dx2(i,j,k) = d2ux_dx2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2ux_dy2(i,j,k) = d2ux_dy2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2ux_dz2(i,j,k) = d2ux_dz2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2uy_dx2(i,j,k) = d2uy_dx2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2uy_dy2(i,j,k) = d2uy_dy2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2uy_dz2(i,j,k) = d2uy_dz2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2uz_dx2(i,j,k) = d2uz_dx2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2uz_dy2(i,j,k) = d2uz_dy2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2uz_dz2(i,j,k) = d2uz_dz2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2bx_dx2(i,j,k) = d2bx_dx2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2bx_dy2(i,j,k) = d2bx_dy2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2bx_dz2(i,j,k) = d2bx_dz2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2by_dx2(i,j,k) = d2by_dx2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2by_dy2(i,j,k) = d2by_dy2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2by_dz2(i,j,k) = d2by_dz2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2bz_dx2(i,j,k) = d2bz_dx2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2bz_dy2(i,j,k) = d2bz_dy2(i,j,k) / dfloat(Nx*Ny*Nz)
      d2bz_dz2(i,j,k) = d2bz_dz2(i,j,k) / dfloat(Nx*Ny*Nz)
    enddo
  enddo
enddo

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      uxRe = (d2ux_dx2(i,j,k) + d2ux_dy2(i,j,k) + d2ux_dz2(i,j,k)) / (rho0 * Re)
      uyRe = (d2uy_dx2(i,j,k) + d2uy_dy2(i,j,k) + d2uy_dz2(i,j,k)) / (rho0 * Re)
      uzRe = (d2uz_dx2(i,j,k) + d2uz_dy2(i,j,k) + d2uz_dz2(i,j,k)) / (rho0 * Re)
      bxRm = (d2bx_dx2(i,j,k) + d2bx_dy2(i,j,k) + d2bx_dz2(i,j,k)) / Rm
      byRm = (d2by_dx2(i,j,k) + d2by_dy2(i,j,k) + d2by_dz2(i,j,k)) / Rm
      bzRm = (d2bz_dx2(i,j,k) + d2bz_dy2(i,j,k) + d2bz_dz2(i,j,k)) / Rm
      dux_dt(i,j,k) = uxRe - B0 * dbx_dz(i,j,k) / rho0
      duy_dt(i,j,k) = uyRe - B0 * dby_dz(i,j,k) / rho0
      duz_dt(i,j,k) = uzRe - B0 * (dbx_dx(i,j,k) + dby_dy(i,j,k) + 2.0d0 * dbz_dz(i,j,k)) / rho0
      dbx_dt(i,j,k) = bxRm - B0 * dux_dz(i,j,k)
      dby_dt(i,j,k) = byRm - B0 * duy_dz(i,j,k)
      dbz_dt(i,j,k) = bzRm + B0 * (dux_dx(i,j,k) + duy_dy(i,j,k))
    enddo
  enddo
enddo

return
end subroutine derive

!====================================================================================

SUBROUTINE rk4(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, &
               dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)
implicit none

INTEGER ( kind = 4 ) Nx, Ny, Nz, Nh, i, j, k
REAL ( kind = 8 ) pi, time, Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm
REAL ( kind = 8 ) ux(Nx,Ny,Nz),dux_dt(Nx,Ny,Nz),k1_ux(Nx,Ny,Nz),k2_ux(Nx,Ny,Nz),k3_ux(Nx,Ny,Nz),k4_ux(Nx,Ny,Nz),dum_ux(Nx,Ny,Nz)
REAL ( kind = 8 ) uy(Nx,Ny,Nz),duy_dt(Nx,Ny,Nz),k1_uy(Nx,Ny,Nz),k2_uy(Nx,Ny,Nz),k3_uy(Nx,Ny,Nz),k4_uy(Nx,Ny,Nz),dum_uy(Nx,Ny,Nz)
REAL ( kind = 8 ) uz(Nx,Ny,Nz),duz_dt(Nx,Ny,Nz),k1_uz(Nx,Ny,Nz),k2_uz(Nx,Ny,Nz),k3_uz(Nx,Ny,Nz),k4_uz(Nx,Ny,Nz),dum_uz(Nx,Ny,Nz)
REAL ( kind = 8 ) bx(Nx,Ny,Nz),dbx_dt(Nx,Ny,Nz),k1_bx(Nx,Ny,Nz),k2_bx(Nx,Ny,Nz),k3_bx(Nx,Ny,Nz),k4_bx(Nx,Ny,Nz),dum_bx(Nx,Ny,Nz)
REAL ( kind = 8 ) by(Nx,Ny,Nz),dby_dt(Nx,Ny,Nz),k1_by(Nx,Ny,Nz),k2_by(Nx,Ny,Nz),k3_by(Nx,Ny,Nz),k4_by(Nx,Ny,Nz),dum_by(Nx,Ny,Nz)
REAL ( kind = 8 ) bz(Nx,Ny,Nz),dbz_dt(Nx,Ny,Nz),k1_bz(Nx,Ny,Nz),k2_bz(Nx,Ny,Nz),k3_bz(Nx,Ny,Nz),k4_bz(Nx,Ny,Nz),dum_bz(Nx,Ny,Nz)

common/comm/Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm

!PRINT*, "RKIdx = 1"

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k1_ux(i,j,k) = dux_dt(i,j,k)
      k1_uy(i,j,k) = duy_dt(i,j,k)
      k1_uz(i,j,k) = duz_dt(i,j,k)
      dum_ux(i,j,k) = ux(i,j,k) + k1_ux(i,j,k) * dt / 2.0d0
      dum_uy(i,j,k) = uy(i,j,k) + k1_uy(i,j,k) * dt / 2.0d0
      dum_uz(i,j,k) = uz(i,j,k) + k1_uz(i,j,k) * dt / 2.0d0
    enddo
  enddo
enddo

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k1_bx(i,j,k) = dbx_dt(i,j,k)
      k1_by(i,j,k) = dby_dt(i,j,k)
      k1_bz(i,j,k) = dbz_dt(i,j,k)
      dum_bx(i,j,k) = bx(i,j,k) + k1_bx(i,j,k) * dt / 2.0d0
      dum_by(i,j,k) = by(i,j,k) + k1_by(i,j,k) * dt / 2.0d0
      dum_bz(i,j,k) = bz(i,j,k) + k1_bz(i,j,k) * dt / 2.0d0
    enddo
  enddo
enddo

CALL derive(Nx, Ny, Nz, Nh, pi, time+dt/2.0d0, dum_ux, dum_uy, dum_uz, dum_bx, dum_by, dum_bz, &
            dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)

!PRINT*, "RKIdx = 2"

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k2_ux(i,j,k) = dux_dt(i,j,k)
      k2_uy(i,j,k) = duy_dt(i,j,k)
      k2_uz(i,j,k) = duz_dt(i,j,k)
      dum_ux(i,j,k) = ux(i,j,k) + k2_ux(i,j,k) * dt/2.0d0
      dum_uy(i,j,k) = uy(i,j,k) + k2_uy(i,j,k) * dt/2.0d0
      dum_uz(i,j,k) = uz(i,j,k) + k2_uz(i,j,k) * dt/2.0d0
    enddo
  enddo
enddo

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k2_bx(i,j,k) = dbx_dt(i,j,k)
      k2_by(i,j,k) = dby_dt(i,j,k)
      k2_bz(i,j,k) = dbz_dt(i,j,k)
      dum_bx(i,j,k) = bx(i,j,k) + k2_bx(i,j,k) * dt/2.0d0
      dum_by(i,j,k) = by(i,j,k) + k2_by(i,j,k) * dt/2.0d0
      dum_bz(i,j,k) = bz(i,j,k) + k2_bz(i,j,k) * dt/2.0d0
    enddo
  enddo
enddo

CALL derive(Nx, Ny, Nz, Nh, pi, time+dt/2.0d0, dum_ux, dum_uy, dum_uz, dum_bx, dum_by, dum_bz, &
            dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)

!PRINT*, "RKIdx = 3"

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k3_ux(i,j,k) = dux_dt(i,j,k)
      k3_uy(i,j,k) = duy_dt(i,j,k)
      k3_uz(i,j,k) = duz_dt(i,j,k)
      dum_ux(i,j,k) = ux(i,j,k) + k3_ux(i,j,k) * dt
      dum_uy(i,j,k) = uy(i,j,k) + k3_uy(i,j,k) * dt
      dum_uz(i,j,k) = uz(i,j,k) + k3_uz(i,j,k) * dt
    enddo
  enddo
enddo

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k3_bx(i,j,k) = dbx_dt(i,j,k)
      k3_by(i,j,k) = dby_dt(i,j,k)
      k3_bz(i,j,k) = dbz_dt(i,j,k)
      dum_bx(i,j,k) = bx(i,j,k) + k3_bx(i,j,k) * dt
      dum_by(i,j,k) = by(i,j,k) + k3_by(i,j,k) * dt
      dum_bz(i,j,k) = bz(i,j,k) + k3_bz(i,j,k) * dt
    enddo
  enddo
enddo

CALL derive(Nx, Ny, Nz, Nh, pi, time+dt, dum_ux, dum_uy, dum_uz, dum_bx, dum_by, dum_bz, &
            dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)

!PRINT*, "RKIdx = 4"

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k4_ux(i,j,k) = dux_dt(i,j,k)
      k4_uy(i,j,k) = duy_dt(i,j,k)
      k4_uz(i,j,k) = duz_dt(i,j,k)
      ux(i,j,k) = ux(i,j,k) + dt/6.0d0*(k1_ux(i,j,k) + 2.0d0*k2_ux(i,j,k) + 2.0d0*k3_ux(i,j,k) + k4_ux(i,j,k))
      uy(i,j,k) = uy(i,j,k) + dt/6.0d0*(k1_uy(i,j,k) + 2.0d0*k2_uy(i,j,k) + 2.0d0*k3_uy(i,j,k) + k4_uy(i,j,k))
      uz(i,j,k) = uz(i,j,k) + dt/6.0d0*(k1_uz(i,j,k) + 2.0d0*k2_uz(i,j,k) + 2.0d0*k3_uz(i,j,k) + k4_uz(i,j,k))
    enddo
  enddo
enddo

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      k4_bx(i,j,k) = dbx_dt(i,j,k)
      k4_by(i,j,k) = dby_dt(i,j,k)
      k4_bz(i,j,k) = dbz_dt(i,j,k)
      bx(i,j,k) = bx(i,j,k) + dt/6.0d0*(k1_bx(i,j,k) + 2.0d0*k2_bx(i,j,k) + 2.0d0*k3_bx(i,j,k) + k4_bx(i,j,k))
      by(i,j,k) = by(i,j,k) + dt/6.0d0*(k1_by(i,j,k) + 2.0d0*k2_by(i,j,k) + 2.0d0*k3_by(i,j,k) + k4_by(i,j,k))
      bz(i,j,k) = bz(i,j,k) + dt/6.0d0*(k1_bz(i,j,k) + 2.0d0*k2_bz(i,j,k) + 2.0d0*k3_bz(i,j,k) + k4_bz(i,j,k))
    enddo
  enddo
enddo

RETURN
end subroutine rk4

!==================================================================================

end program linmhd3d

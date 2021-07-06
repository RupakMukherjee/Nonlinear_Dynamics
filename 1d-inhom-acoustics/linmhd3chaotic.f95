program linmhd3d

implicit none

include "fftw3.f"

! Define Grid Size.
INTEGER ( kind = 4 ), parameter :: Nx = 32
INTEGER ( kind = 4 ), parameter :: Ny = 32
INTEGER ( kind = 4 ), parameter :: Nz = 32
INTEGER ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

REAL ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

REAL ( kind = 8 ) x(Nx),y(Ny),z(Nz)
REAL ( kind = 8 ) ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
REAL ( kind = 8 ) bx(Nx,Ny,Nz),by(Nx,Ny,Nz),bz(Nx,Ny,Nz)
REAL ( kind = 8 ) B0x(Nx,Ny,Nz),B0y(Nx,Ny,Nz),B0z(Nx,Ny,Nz)
REAL ( kind = 8 ) B0x_dum(Nx,Ny,Nz),B0y_dum(Nx,Ny,Nz),B0z_dum(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0x_dx(Nx,Ny,Nz),dB0x_dy(Nx,Ny,Nz),dB0x_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0y_dx(Nx,Ny,Nz),dB0y_dy(Nx,Ny,Nz),dB0y_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0z_dx(Nx,Ny,Nz),dB0z_dy(Nx,Ny,Nz),dB0z_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dux_dt(Nx,Ny,Nz),duy_dt(Nx,Ny,Nz),duz_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) dbx_dt(Nx,Ny,Nz),dby_dt(Nx,Ny,Nz),dbz_dt(Nx,Ny,Nz)
COMPLEX ( kind = 8 ) B0xk(Nh,Ny,Nz), B0yk(Nh,Ny,Nz), B0zk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikx_B0xk(Nh,Ny,Nz), ikx_B0yk(Nh,Ny,Nz), ikx_B0zk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) iky_B0xk(Nh,Ny,Nz), iky_B0yk(Nh,Ny,Nz), iky_B0zk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikz_B0xk(Nh,Ny,Nz), ikz_B0yk(Nh,Ny,Nz), ikz_B0zk(Nh,Ny,Nz)

REAL ( kind = 8 ), dimension (:), allocatable :: freq, abs_freq_k
COMPLEX ( kind = 8 ), dimension (:), allocatable :: freq_k

INTEGER ( kind = 8 ) :: planf, planb

INTEGER ( kind = 4 ) i,j,k,t,M
REAL ( kind = 8 ) :: time, time_min, time_max, dt, f, rho0, B0, Re, Rm, epsa, A, B, C
REAL ( kind = 8 ) :: dx, dy, dz, Lx, Ly, Lz, kx, ky, kz, uEnergy, bEnergy, Energy
REAL ( kind = 8 ) :: uxEnergy, uyEnergy, uzEnergy, bxEnergy, byEnergy, bzEnergy

INTEGER ( kind = 8 ), dimension (:), allocatable :: sort

COMMON/coeff/Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm

open(unit=10,file='Energy.dat',status='unknown')
open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
open(unit=30,file='FFT.dat',status='unknown')

!===================== USER INPUTS ============================================

Lx = 2.0d0*pi
Ly = 2.0d0*pi
Lz = 2.0d0*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)
dz = Lz/dfloat(Nz)

time_min = 0.00d0
time_max = 100.00d0
dt = 0.010d0

M = nint(time_max/dt) - int(time_min/dt) + 1
allocate(freq(M))
allocate(abs_freq_k(M/2+1))
allocate(freq_k(M/2+1))
allocate(sort(M/2+1))

f = 1.0d0

Re = 1e3
Rm = 1e3

rho0 = 1.0d0
B0 = 1.0d0

epsa = 0.010d0
 A = 1.0d0
 B = 1.0d0
 C = 1.0d0

! Grid Generation.
do i = 1, Nx
  x(i) = dfloat(i-1) * dx
  do j = 1, Ny
    y(j) = dfloat(j-1) * dy
    do k = 1, Nz
      z(k) = dfloat(k-1) * dz
      ux(i,j,k) = 0.0d0 !dsin(f*z(k))
      uy(i,j,k) = dsin(f*z(k))
      uz(i,j,k) = 0.0d0 !dcos(f*x(i))
      bx(i,j,k) = 0.0d0 !+dcos(f*x(i))*dsin(f*z(k))
      by(i,j,k) = dcos(f*z(k))
      bz(i,j,k) = 0.0d0 !-dsin(f*x(i))*dcos(f*z(k))
      B0x(i,j,k) = epsa*(A*dsin(f*z(k)) + C*dcos(f*y(j)))
      B0y(i,j,k) = epsa*(B*dsin(f*x(i)) + A*dcos(f*z(k)))
      B0z(i,j,k) = B0 + epsa*(C*dsin(f*y(j)) + B*dcos(f*x(i)))
      write(15,*) 0.0d0, x(i),y(j),z(k),ux(i,j,k),uy(i,j,k),uz(i,j,k),bx(i,j,k),by(i,j,k),bz(i,j,k)
    enddo
  enddo
enddo

!======================= MAIN PROGRAM =====================================================

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      B0x_dum(i,j,k) = B0x(i,j,k)
      B0y_dum(i,j,k) = B0y(i,j,k)
      B0z_dum(i,j,k) = B0z(i,j,k)
    enddo
  enddo
enddo

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, B0x_dum, B0xk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, B0y_dum, B0yk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_3d_(planf, Nx, Ny, Nz, B0z_dum, B0zk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

do i = 1, Nh
  kx = 2.0d0 * pi * dfloat(i-1) / Lx
  do j = 1, Ny/2
    ky = 2.0d0 * pi * dfloat(j-1) / Ly
    do k = 1, Nz/2
      kz = 2.0d0 * pi * dfloat(k-1) / Lz
      ikx_B0xk(i,j,k) = (0.0d0,1.0d0) * kx * B0xk(i,j,k)
      ikx_B0yk(i,j,k) = (0.0d0,1.0d0) * kx * B0yk(i,j,k)
      ikx_B0zk(i,j,k) = (0.0d0,1.0d0) * kx * B0zk(i,j,k)
      iky_B0xk(i,j,k) = (0.0d0,1.0d0) * ky * B0xk(i,j,k)
      iky_B0yk(i,j,k) = (0.0d0,1.0d0) * ky * B0yk(i,j,k)
      iky_B0zk(i,j,k) = (0.0d0,1.0d0) * ky * B0zk(i,j,k)
      ikz_B0xk(i,j,k) = (0.0d0,1.0d0) * kz * B0xk(i,j,k)
      ikz_B0yk(i,j,k) = (0.0d0,1.0d0) * kz * B0yk(i,j,k)
      ikz_B0zk(i,j,k) = (0.0d0,1.0d0) * kz * B0zk(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kz = 2.0d0 * pi * dfloat((k-1) - Nz) / Lz
      ikx_B0xk(i,j,k) = (0.0d0,1.0d0) * kx * B0xk(i,j,k)
      ikx_B0yk(i,j,k) = (0.0d0,1.0d0) * kx * B0yk(i,j,k)
      ikx_B0zk(i,j,k) = (0.0d0,1.0d0) * kx * B0zk(i,j,k)
      iky_B0xk(i,j,k) = (0.0d0,1.0d0) * ky * B0xk(i,j,k)
      iky_B0yk(i,j,k) = (0.0d0,1.0d0) * ky * B0yk(i,j,k)
      iky_B0zk(i,j,k) = (0.0d0,1.0d0) * ky * B0zk(i,j,k)
      ikz_B0xk(i,j,k) = (0.0d0,1.0d0) * kz * B0xk(i,j,k)
      ikz_B0yk(i,j,k) = (0.0d0,1.0d0) * kz * B0yk(i,j,k)
      ikz_B0zk(i,j,k) = (0.0d0,1.0d0) * kz * B0zk(i,j,k)
    enddo
  enddo
  do j = Ny/2+1, Ny
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    do k = 1, Nz/2
      kz = 2.0d0*pi*float(k-1)/Lz
      ikx_B0xk(i,j,k) = (0.0d0,1.0d0) * kx * B0xk(i,j,k)
      ikx_B0yk(i,j,k) = (0.0d0,1.0d0) * kx * B0yk(i,j,k)
      ikx_B0zk(i,j,k) = (0.0d0,1.0d0) * kx * B0zk(i,j,k)
      iky_B0xk(i,j,k) = (0.0d0,1.0d0) * ky * B0xk(i,j,k)
      iky_B0yk(i,j,k) = (0.0d0,1.0d0) * ky * B0yk(i,j,k)
      iky_B0zk(i,j,k) = (0.0d0,1.0d0) * ky * B0zk(i,j,k)
      ikz_B0xk(i,j,k) = (0.0d0,1.0d0) * kz * B0xk(i,j,k)
      ikz_B0yk(i,j,k) = (0.0d0,1.0d0) * kz * B0yk(i,j,k)
      ikz_B0zk(i,j,k) = (0.0d0,1.0d0) * kz * B0zk(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      ikx_B0xk(i,j,k) = (0.0d0,1.0d0) * kx * B0xk(i,j,k)
      ikx_B0yk(i,j,k) = (0.0d0,1.0d0) * kx * B0yk(i,j,k)
      ikx_B0zk(i,j,k) = (0.0d0,1.0d0) * kx * B0zk(i,j,k)
      iky_B0xk(i,j,k) = (0.0d0,1.0d0) * ky * B0xk(i,j,k)
      iky_B0yk(i,j,k) = (0.0d0,1.0d0) * ky * B0yk(i,j,k)
      iky_B0zk(i,j,k) = (0.0d0,1.0d0) * ky * B0zk(i,j,k)
      ikz_B0xk(i,j,k) = (0.0d0,1.0d0) * kz * B0xk(i,j,k)
      ikz_B0yk(i,j,k) = (0.0d0,1.0d0) * kz * B0yk(i,j,k)
      ikz_B0zk(i,j,k) = (0.0d0,1.0d0) * kz * B0zk(i,j,k)
    enddo
  enddo
enddo

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_B0xk, dB0x_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_B0yk, dB0y_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_B0zk, dB0z_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_B0xk, dB0x_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_B0yk, dB0y_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_B0zk, dB0z_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_B0xk, dB0x_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_B0yk, dB0y_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_B0zk, dB0z_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      dB0x_dx(i,j,k) = dB0x_dx(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0y_dx(i,j,k) = dB0y_dx(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0z_dx(i,j,k) = dB0z_dx(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0x_dy(i,j,k) = dB0x_dy(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0y_dy(i,j,k) = dB0y_dy(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0z_dy(i,j,k) = dB0z_dy(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0x_dz(i,j,k) = dB0x_dz(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0y_dz(i,j,k) = dB0y_dz(i,j,k) / dfloat(Nx*Ny*Nz)
      dB0z_dz(i,j,k) = dB0z_dz(i,j,k) / dfloat(Nx*Ny*Nz)
    enddo
  enddo
enddo

DO time = time_min, time_max, dt

Energy = 0.0d0
uEnergy = 0.0d0
uxEnergy = 0.0d0
uyEnergy = 0.0d0
uzEnergy = 0.0d0
bEnergy = 0.0d0
bxEnergy = 0.0d0
byEnergy = 0.0d0
bzEnergy = 0.0d0

t = nint(time/dt) - dint(time_min/dt)

!Calculating the time evolution in real-space...

CALL derive(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, B0x, B0y, B0z, &
            dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
            dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)
CALL rk4(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, B0x, B0y, B0z, &
         dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
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

freq(t+1) = 0.0d0

DO i = 1, Nx
    do j = 1, Ny
      freq(t+1) = freq(t+1) + by(i,j,Nz/2) / dfloat(Nx*Ny)
  enddo
enddo

!IF ( MOD(t,100) == 0 ) THEN
  DO i = 1, Nx
    do j = 1, Ny
      do k = 1, Nz
        uxEnergy = uxEnergy + (ux(i,j,k) * ux(i,j,k)) / dfloat(Nx*Ny*Nz)
        uyEnergy = uyEnergy + (uy(i,j,k) * uy(i,j,k)) / dfloat(Nx*Ny*Nz)
        uzEnergy = uzEnergy + (uz(i,j,k) * uz(i,j,k)) / dfloat(Nx*Ny*Nz)
        bxEnergy = bxEnergy + (bx(i,j,k) * bx(i,j,k)) / dfloat(Nx*Ny*Nz)
        byEnergy = byEnergy + (by(i,j,k) * by(i,j,k)) / dfloat(Nx*Ny*Nz)
        bzEnergy = bzEnergy + (bz(i,j,k) * bz(i,j,k)) / dfloat(Nx*Ny*Nz)
        uEnergy = uEnergy + (ux(i,j,k) * ux(i,j,k) + uy(i,j,k) * uy(i,j,k) + uz(i,j,k) * uz(i,j,k)) / dfloat(Nx*Ny*Nz)
        bEnergy = bEnergy + (bx(i,j,k) * bx(i,j,k) + by(i,j,k) * by(i,j,k) + bz(i,j,k) * bz(i,j,k)) / dfloat(Nx*Ny*Nz)
        Energy = uEnergy + bEnergy
      enddo
    enddo
  ENDDO
  WRITE(10,*) time, uEnergy, bEnergy, Energy, uxEnergy, uyEnergy, uzEnergy, bxEnergy, byEnergy, bzEnergy!log(bEnergy)!, freq(t)
!ENDIF

ENDDO ! time

call dfftw_plan_dft_r2c_1d_(planf, M, freq, freq_k, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

do i = 1,M/2+1,1
  freq_k(i) = 2.0d0*freq_k(i)/dfloat(M)
  abs_freq_k(i) = abs(freq_k(i))
  sort(i) = i
  !write(30,*) 2.0d0*pi*dfloat(i)/(dfloat(int(1.0d0/dt))),abs(freq_k(i))
enddo

call exchange (abs_freq_k, sort, M/2+1)

write(30,*) 2.0d0*pi*dfloat(sort(1))/(f*dfloat(int(1.0d0/dt)))

contains
!===================================================================================

SUBROUTINE derive(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, B0x, B0y, B0z, &
                  dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
                  dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)

implicit none

INTEGER ( kind = 4 ) Nx, Ny, Nz, Nh, i, j, k
REAL ( kind = 8 ) Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm
REAL ( kind = 8 ) pi, time, kx, ky, kz
REAL ( kind = 8 ) xB0db, xbdB0, yB0db, ybdB0, zB0db, zbdB0
REAL ( kind = 8 ) B0dux, udB0x, B0duy, udB0y, B0duz, udB0z
REAL ( kind = 8 ) uxRe, uyRe, uzRe, bxRm, byRm, bzRm
REAL ( kind = 8 ) ux(Nx,Ny,Nz), ux_dum(Nx,Ny,Nz), dux_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) uy(Nx,Ny,Nz), uy_dum(Nx,Ny,Nz), duy_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) uz(Nx,Ny,Nz), uz_dum(Nx,Ny,Nz), duz_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) bx(Nx,Ny,Nz), bx_dum(Nx,Ny,Nz), dbx_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) by(Nx,Ny,Nz), by_dum(Nx,Ny,Nz), dby_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) bz(Nx,Ny,Nz), bz_dum(Nx,Ny,Nz), dbz_dt(Nx,Ny,Nz)
REAL ( kind = 8 ) B0x(Nx,Ny,Nz),B0y(Nx,Ny,Nz),B0z(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0x_dx(Nx,Ny,Nz),dB0x_dy(Nx,Ny,Nz),dB0x_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0y_dx(Nx,Ny,Nz),dB0y_dy(Nx,Ny,Nz),dB0y_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0z_dx(Nx,Ny,Nz),dB0z_dy(Nx,Ny,Nz),dB0z_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dux_dx(Nx,Ny,Nz), duy_dx(Nx,Ny,Nz), duz_dx(Nx,Ny,Nz)
REAL ( kind = 8 ) dux_dy(Nx,Ny,Nz), duy_dy(Nx,Ny,Nz), duz_dy(Nx,Ny,Nz)
REAL ( kind = 8 ) dux_dz(Nx,Ny,Nz), duy_dz(Nx,Ny,Nz), duz_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dbx_dy(Nx,Ny,Nz), dbx_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dby_dx(Nx,Ny,Nz), dby_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dbz_dx(Nx,Ny,Nz), dbz_dy(Nx,Ny,Nz)
REAL ( kind = 8 ) d2ux_dx2(Nx,Ny,Nz), d2ux_dy2(Nx,Ny,Nz), d2ux_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2uy_dx2(Nx,Ny,Nz), d2uy_dy2(Nx,Ny,Nz), d2uy_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2uz_dx2(Nx,Ny,Nz), d2uz_dy2(Nx,Ny,Nz), d2uz_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2bx_dx2(Nx,Ny,Nz), d2bx_dy2(Nx,Ny,Nz), d2bx_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2by_dx2(Nx,Ny,Nz), d2by_dy2(Nx,Ny,Nz), d2by_dz2(Nx,Ny,Nz)
REAL ( kind = 8 ) d2bz_dx2(Nx,Ny,Nz), d2bz_dy2(Nx,Ny,Nz), d2bz_dz2(Nx,Ny,Nz)
COMPLEX ( kind = 8 ) uxk(Nh,Ny,Nz), uyk(Nh,Ny,Nz), uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) bxk(Nh,Ny,Nz), byk(Nh,Ny,Nz), bzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikx_uxk(Nh,Ny,Nz), ikx_uyk(Nh,Ny,Nz), ikx_uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) iky_uxk(Nh,Ny,Nz), iky_uyk(Nh,Ny,Nz), iky_uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikz_uxk(Nh,Ny,Nz), ikz_uyk(Nh,Ny,Nz), ikz_uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) iky_bxk(Nh,Ny,Nz), ikz_bxk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikx_byk(Nh,Ny,Nz), ikz_byk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) ikx_bzk(Nh,Ny,Nz), iky_bzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_uxk(Nh,Ny,Nz), ky2_uxk(Nh,Ny,Nz), kz2_uxk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_uyk(Nh,Ny,Nz), ky2_uyk(Nh,Ny,Nz), kz2_uyk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_uzk(Nh,Ny,Nz), ky2_uzk(Nh,Ny,Nz), kz2_uzk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_bxk(Nh,Ny,Nz), ky2_bxk(Nh,Ny,Nz), kz2_bxk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_byk(Nh,Ny,Nz), ky2_byk(Nh,Ny,Nz), kz2_byk(Nh,Ny,Nz)
COMPLEX ( kind = 8 ) kx2_bzk(Nh,Ny,Nz), ky2_bzk(Nh,Ny,Nz), kz2_bzk(Nh,Ny,Nz)

COMMON/coeff/Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm

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
  do j = 1, Ny/2
    ky = 2.0d0 * pi * dfloat(j-1) / Ly
    do k = 1, Nz/2
      kz = 2.0d0 * pi * dfloat(k-1) / Lz
      ikx_uxk(i,j,k) = (0.0d0,1.0d0) * kx * uxk(i,j,k)
      ikx_uyk(i,j,k) = (0.0d0,1.0d0) * kx * uyk(i,j,k)
      ikx_uzk(i,j,k) = (0.0d0,1.0d0) * kx * uzk(i,j,k)
      iky_uxk(i,j,k) = (0.0d0,1.0d0) * ky * uxk(i,j,k)
      iky_uyk(i,j,k) = (0.0d0,1.0d0) * ky * uyk(i,j,k)
      iky_uzk(i,j,k) = (0.0d0,1.0d0) * ky * uzk(i,j,k)
      ikz_uxk(i,j,k) = (0.0d0,1.0d0) * kz * uxk(i,j,k)
      ikz_uyk(i,j,k) = (0.0d0,1.0d0) * kz * uyk(i,j,k)
      ikz_uzk(i,j,k) = (0.0d0,1.0d0) * kz * uzk(i,j,k)
      iky_bxk(i,j,k) = (0.0d0,1.0d0) * ky * bxk(i,j,k)
      ikz_bxk(i,j,k) = (0.0d0,1.0d0) * kz * bxk(i,j,k)
      ikx_byk(i,j,k) = (0.0d0,1.0d0) * kx * byk(i,j,k)
      ikz_byk(i,j,k) = (0.0d0,1.0d0) * kz * byk(i,j,k)
      ikx_bzk(i,j,k) = (0.0d0,1.0d0) * kx * bzk(i,j,k)
      iky_bzk(i,j,k) = (0.0d0,1.0d0) * ky * bzk(i,j,k)
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
      ikx_uxk(i,j,k) = (0.0d0,1.0d0) * kx * uxk(i,j,k)
      ikx_uyk(i,j,k) = (0.0d0,1.0d0) * kx * uyk(i,j,k)
      ikx_uzk(i,j,k) = (0.0d0,1.0d0) * kx * uzk(i,j,k)
      iky_uxk(i,j,k) = (0.0d0,1.0d0) * ky * uxk(i,j,k)
      iky_uyk(i,j,k) = (0.0d0,1.0d0) * ky * uyk(i,j,k)
      iky_uzk(i,j,k) = (0.0d0,1.0d0) * ky * uzk(i,j,k)
      ikz_uxk(i,j,k) = (0.0d0,1.0d0) * kz * uxk(i,j,k)
      ikz_uyk(i,j,k) = (0.0d0,1.0d0) * kz * uyk(i,j,k)
      ikz_uzk(i,j,k) = (0.0d0,1.0d0) * kz * uzk(i,j,k)
      iky_bxk(i,j,k) = (0.0d0,1.0d0) * ky * bxk(i,j,k)
      ikz_bxk(i,j,k) = (0.0d0,1.0d0) * kz * bxk(i,j,k)
      ikx_byk(i,j,k) = (0.0d0,1.0d0) * kx * byk(i,j,k)
      ikz_byk(i,j,k) = (0.0d0,1.0d0) * kz * byk(i,j,k)
      ikx_bzk(i,j,k) = (0.0d0,1.0d0) * kx * bzk(i,j,k)
      iky_bzk(i,j,k) = (0.0d0,1.0d0) * ky * bzk(i,j,k)
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
  do j = Ny/2+1, Ny
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    do k = 1, Nz/2
      kz = 2.0d0*pi*float(k-1)/Lz
      ikx_uxk(i,j,k) = (0.0d0,1.0d0) * kx * uxk(i,j,k)
      ikx_uyk(i,j,k) = (0.0d0,1.0d0) * kx * uyk(i,j,k)
      ikx_uzk(i,j,k) = (0.0d0,1.0d0) * kx * uzk(i,j,k)
      iky_uxk(i,j,k) = (0.0d0,1.0d0) * ky * uxk(i,j,k)
      iky_uyk(i,j,k) = (0.0d0,1.0d0) * ky * uyk(i,j,k)
      iky_uzk(i,j,k) = (0.0d0,1.0d0) * ky * uzk(i,j,k)
      ikz_uxk(i,j,k) = (0.0d0,1.0d0) * kz * uxk(i,j,k)
      ikz_uyk(i,j,k) = (0.0d0,1.0d0) * kz * uyk(i,j,k)
      ikz_uzk(i,j,k) = (0.0d0,1.0d0) * kz * uzk(i,j,k)
      iky_bxk(i,j,k) = (0.0d0,1.0d0) * ky * bxk(i,j,k)
      ikz_bxk(i,j,k) = (0.0d0,1.0d0) * kz * bxk(i,j,k)
      ikx_byk(i,j,k) = (0.0d0,1.0d0) * kx * byk(i,j,k)
      ikz_byk(i,j,k) = (0.0d0,1.0d0) * kz * byk(i,j,k)
      ikx_bzk(i,j,k) = (0.0d0,1.0d0) * kx * bzk(i,j,k)
      iky_bzk(i,j,k) = (0.0d0,1.0d0) * ky * bzk(i,j,k)
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
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      ikx_uxk(i,j,k) = (0.0d0,1.0d0) * kx * uxk(i,j,k)
      ikx_uyk(i,j,k) = (0.0d0,1.0d0) * kx * uyk(i,j,k)
      ikx_uzk(i,j,k) = (0.0d0,1.0d0) * kx * uzk(i,j,k)
      iky_uxk(i,j,k) = (0.0d0,1.0d0) * ky * uxk(i,j,k)
      iky_uyk(i,j,k) = (0.0d0,1.0d0) * ky * uyk(i,j,k)
      iky_uzk(i,j,k) = (0.0d0,1.0d0) * ky * uzk(i,j,k)
      ikz_uxk(i,j,k) = (0.0d0,1.0d0) * kz * uxk(i,j,k)
      ikz_uyk(i,j,k) = (0.0d0,1.0d0) * kz * uyk(i,j,k)
      ikz_uzk(i,j,k) = (0.0d0,1.0d0) * kz * uzk(i,j,k)
      iky_bxk(i,j,k) = (0.0d0,1.0d0) * ky * bxk(i,j,k)
      ikz_bxk(i,j,k) = (0.0d0,1.0d0) * kz * bxk(i,j,k)
      ikx_byk(i,j,k) = (0.0d0,1.0d0) * kx * byk(i,j,k)
      ikz_byk(i,j,k) = (0.0d0,1.0d0) * kz * byk(i,j,k)
      ikx_bzk(i,j,k) = (0.0d0,1.0d0) * kx * bzk(i,j,k)
      iky_bzk(i,j,k) = (0.0d0,1.0d0) * ky * bzk(i,j,k)
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

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_uxk, dux_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_uyk, duy_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_uzk, duz_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_uxk, dux_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_uyk, duy_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_uzk, duz_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_uxk, dux_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_uyk, duy_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_uzk, duz_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_bxk, dbx_dy, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_bxk, dbx_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_byk, dby_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikz_byk, dby_dz, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, ikx_bzk, dbz_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_3d_(planb, Nx, Ny, Nz, iky_bzk, dbz_dy, FFTW_ESTIMATE)
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
      dux_dx(i,j,k)   = dux_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      duy_dx(i,j,k)   = duy_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      duz_dx(i,j,k)   = duz_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      dux_dy(i,j,k)   = dux_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
      duy_dy(i,j,k)   = duy_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
      duz_dy(i,j,k)   = duz_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
      dux_dz(i,j,k)   = dux_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      duy_dz(i,j,k)   = duy_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      duz_dz(i,j,k)   = duz_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbx_dy(i,j,k)   = dbx_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbx_dz(i,j,k)   = dbx_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      dby_dx(i,j,k)   = dby_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      dby_dz(i,j,k)   = dby_dz(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbz_dx(i,j,k)   = dbz_dx(i,j,k)   / dfloat(Nx*Ny*Nz)
      dbz_dy(i,j,k)   = dbz_dy(i,j,k)   / dfloat(Nx*Ny*Nz)
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
      uxRe  = (d2ux_dx2(i,j,k) + d2ux_dy2(i,j,k) + d2ux_dz2(i,j,k)) / (rho0 * Re)
      uyRe  = (d2uy_dx2(i,j,k) + d2uy_dy2(i,j,k) + d2uy_dz2(i,j,k)) / (rho0 * Re)
      uzRe  = (d2uz_dx2(i,j,k) + d2uz_dy2(i,j,k) + d2uz_dz2(i,j,k)) / (rho0 * Re)
      bxRm  = (d2bx_dx2(i,j,k) + d2bx_dy2(i,j,k) + d2bx_dz2(i,j,k)) / Rm
      byRm  = (d2by_dx2(i,j,k) + d2by_dy2(i,j,k) + d2by_dz2(i,j,k)) / Rm
      bzRm  = (d2bz_dx2(i,j,k) + d2bz_dy2(i,j,k) + d2bz_dz2(i,j,k)) / Rm
      xB0db = B0y(i,j,k) * (dbx_dy(i,j,k)  - dby_dx(i,j,k))  + B0z(i,j,k) * (dbx_dz(i,j,k)  - dbz_dx(i,j,k))
      yB0db = B0x(i,j,k) * (dby_dx(i,j,k)  - dbx_dy(i,j,k))  + B0z(i,j,k) * (dby_dz(i,j,k)  - dbz_dy(i,j,k))
      zB0db = B0x(i,j,k) * (dbz_dx(i,j,k)  - dbx_dz(i,j,k))  + B0y(i,j,k) * (dbz_dy(i,j,k)  - dby_dz(i,j,k))
      xbdB0 = by(i,j,k)  * (dB0x_dy(i,j,k) - dB0y_dx(i,j,k)) + bz(i,j,k)  * (dB0x_dz(i,j,k) - dB0z_dx(i,j,k))
      ybdB0 = bx(i,j,k)  * (dB0y_dx(i,j,k) - dB0x_dy(i,j,k)) + bz(i,j,k)  * (dB0y_dz(i,j,k) - dB0z_dy(i,j,k))
      zbdB0 = bx(i,j,k)  * (dB0z_dx(i,j,k) - dB0x_dz(i,j,k)) + by(i,j,k)  * (dB0z_dy(i,j,k) - dB0y_dz(i,j,k))
      B0dux = B0x(i,j,k) * dux_dx(i,j,k)  + B0y(i,j,k) * dux_dy(i,j,k)  + B0z(i,j,k) * dux_dz(i,j,k)
      B0duy = B0x(i,j,k) * duy_dx(i,j,k)  + B0y(i,j,k) * duy_dy(i,j,k)  + B0z(i,j,k) * duy_dz(i,j,k)
      B0duz = B0x(i,j,k) * duz_dx(i,j,k)  + B0y(i,j,k) * duz_dy(i,j,k)  + B0z(i,j,k) * duz_dz(i,j,k)
      udB0x = ux(i,j,k)  * dB0x_dx(i,j,k) + uy(i,j,k)  * dB0x_dy(i,j,k) + uz(i,j,k)  * dB0x_dz(i,j,k)
      udB0y = ux(i,j,k)  * dB0y_dx(i,j,k) + uy(i,j,k)  * dB0y_dy(i,j,k) + uz(i,j,k)  * dB0y_dz(i,j,k)
      udB0z = ux(i,j,k)  * dB0z_dx(i,j,k) + uy(i,j,k)  * dB0z_dy(i,j,k) + uz(i,j,k)  * dB0z_dz(i,j,k)
      dux_dt(i,j,k) = uxRe - (xB0db + xbdB0) / rho0
      duy_dt(i,j,k) = uyRe - (yB0db + ybdB0) / rho0
      duz_dt(i,j,k) = uzRe - (zB0db + zbdB0) / rho0
      dbx_dt(i,j,k) = bxRm - B0dux  + udB0x
      dby_dt(i,j,k) = byRm - B0duy  + udB0y
      dbz_dt(i,j,k) = bzRm - B0duz  + udB0z
    enddo
  enddo
enddo

return
end subroutine derive

!====================================================================================

SUBROUTINE rk4(Nx, Ny, Nz, Nh, pi, time, ux, uy, uz, bx, by, bz, B0x, B0y, B0z, &
               dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
               dux_dt, duy_dt, duz_dt, dbx_dt, dby_dt, dbz_dt)
implicit none

INTEGER ( kind = 4 ) Nx, Ny, Nz, Nh, i, j, k
REAL ( kind = 8 ) pi, time, Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm
REAL ( kind = 8 ) B0x(Nx,Ny,Nz),B0y(Nx,Ny,Nz),B0z(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0x_dx(Nx,Ny,Nz),dB0x_dy(Nx,Ny,Nz),dB0x_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0y_dx(Nx,Ny,Nz),dB0y_dy(Nx,Ny,Nz),dB0y_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) dB0z_dx(Nx,Ny,Nz),dB0z_dy(Nx,Ny,Nz),dB0z_dz(Nx,Ny,Nz)
REAL ( kind = 8 ) ux(Nx,Ny,Nz),dux_dt(Nx,Ny,Nz),k1_ux(Nx,Ny,Nz),k2_ux(Nx,Ny,Nz),k3_ux(Nx,Ny,Nz),k4_ux(Nx,Ny,Nz),dum_ux(Nx,Ny,Nz)
REAL ( kind = 8 ) uy(Nx,Ny,Nz),duy_dt(Nx,Ny,Nz),k1_uy(Nx,Ny,Nz),k2_uy(Nx,Ny,Nz),k3_uy(Nx,Ny,Nz),k4_uy(Nx,Ny,Nz),dum_uy(Nx,Ny,Nz)
REAL ( kind = 8 ) uz(Nx,Ny,Nz),duz_dt(Nx,Ny,Nz),k1_uz(Nx,Ny,Nz),k2_uz(Nx,Ny,Nz),k3_uz(Nx,Ny,Nz),k4_uz(Nx,Ny,Nz),dum_uz(Nx,Ny,Nz)
REAL ( kind = 8 ) bx(Nx,Ny,Nz),dbx_dt(Nx,Ny,Nz),k1_bx(Nx,Ny,Nz),k2_bx(Nx,Ny,Nz),k3_bx(Nx,Ny,Nz),k4_bx(Nx,Ny,Nz),dum_bx(Nx,Ny,Nz)
REAL ( kind = 8 ) by(Nx,Ny,Nz),dby_dt(Nx,Ny,Nz),k1_by(Nx,Ny,Nz),k2_by(Nx,Ny,Nz),k3_by(Nx,Ny,Nz),k4_by(Nx,Ny,Nz),dum_by(Nx,Ny,Nz)
REAL ( kind = 8 ) bz(Nx,Ny,Nz),dbz_dt(Nx,Ny,Nz),k1_bz(Nx,Ny,Nz),k2_bz(Nx,Ny,Nz),k3_bz(Nx,Ny,Nz),k4_bz(Nx,Ny,Nz),dum_bz(Nx,Ny,Nz)

COMMON/coeff/Lx, Ly, Lz, dx, dy, dz, dt, rho0, B0, Re, Rm

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

CALL derive(Nx, Ny, Nz, Nh, pi, time+dt/2.0d0, dum_ux, dum_uy, dum_uz, dum_bx, dum_by, dum_bz, B0x, B0y, B0z, &
            dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
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

CALL derive(Nx, Ny, Nz, Nh, pi, time+dt/2.0d0, dum_ux, dum_uy, dum_uz, dum_bx, dum_by, dum_bz, B0x, B0y, B0z, &
            dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
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

CALL derive(Nx, Ny, Nz, Nh, pi, time+dt, dum_ux, dum_uy, dum_uz, dum_bx, dum_by, dum_bz, B0x, B0y, B0z, &
            dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz, dB0z_dx, dB0z_dy, dB0z_dz, &
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

!================== SORTING ===========================================

subroutine exchange (b, s, n)

integer, intent (in) ::  n
real*8, dimension (n), intent (inout) :: b
integer*8, dimension (n), intent (inout) :: s
real*8 :: c
integer*8 ::  i, j, d

do j = n-1, 1, -1
  do i = 1, j, 1
    if (b(i) < b(i+1)) then
      c=b(i)
      b(i)=b(i+1)
      b(i+1)=c
      d = s(i)
      s(i) = s(i+1)
      s(i+1) = d
    endif
  end do
end do

end subroutine exchange

!==================================================================================

end program linmhd3d

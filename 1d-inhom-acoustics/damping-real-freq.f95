! To run the program : gfortran -w -I/usr/local/include -L/usr/local/lib <filename.f95> -lfftw3 -lm; ./a.out
! It uses Runge-Kutta-4 Solver
program inhomogeneous_wave_equation
implicit none

include "fftw3.f"

integer ( kind = 4 ), parameter :: N = 1024
integer ( kind = 4 ), parameter :: Nh = N/2+1
integer ( kind = 4 ), parameter :: M = 10001
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

DOUBLE PRECISION :: x,u,b,du_dt,db_dt
DIMENSION x(N),u(N),b(N),du_dt(N),db_dt(N)

REAL ( kind = 8 ) freq(M), abs_freq_k(M/2+1)
COMPLEX ( kind = 8 ) freq_k(M/2+1)

integer ( kind = 8 ) :: planf, planb

integer ( kind = 4) i,t
REAL ( kind = 8 ) :: time, time_min, time_max, dt, rho0, B0, Re, Rm, dx, L, k, uEnergy, bEnergy!, Energy

integer*8 sort(M/2+1)

COMMON/comm/L, dx, dt, rho0, B0, Re, Rm

integer,parameter :: seed = 99999999
CALL srand(seed)

open(unit=10,file='Energy.dat',status='unknown')
open(unit=30,file='FFT.dat',status='unknown')

!===================== USER INPUTS ============================================

L = 10.0d0*pi
dx = L/dfloat(N)

time_min = 0.00d0
time_max = 100.0d0
dt = 0.010d0

k = 10.0d0

Re = 1e7
Rm = 1e7

rho0 = 1.0d0
B0 = 1.0d0

do i=1,N
  x(i) = dfloat(i-1) * dx
  u(i) = dsin(k*x(i))
  b(i) = dcos(k*x(i))
!  write (100,*) 0, x(i), u(i), b(i)
enddo
!write (100,*) 0, L, u(1), b(1)

!==================== MAIN PROGRAM ============================================

DO time = time_min, time_max, dt

!Energy = 0.0d0
!uEnergy = 0.0d0
bEnergy = 0.0d0

t = nint(time/dt) - int(time_min/dt)

!Calculating the time evolution in real-space...

CALL derive(N, Nh, pi, time, u, b, du_dt, db_dt)
CALL rk4(N, Nh, pi, time, u, b, du_dt, db_dt)

!PRINT*, "Printing Data"

!IF ( t /= 0 .and. MOD(t,100) == 0 ) THEN
!   DO i = 1, N
!    WRITE (t+100,*) t, x(i), u(i), b(i)
!  enddo
!  write (t+100,*) t, L, u(1), b(1)
!endif

!IF ( MOD(t,100) == 0 ) THEN
  DO i = 1, N
    !uEnergy = uEnergy + u(i) * u(i) / dfloat(2*N)
    bEnergy = bEnergy + b(i) * b(i) / dfloat(N)
    !Energy = uEnergy + bEnergy
  ENDDO

  freq(t+1) = b(N/2)
  WRITE(10,*) time, log(bEnergy), freq(t)
!ENDIF

enddo ! time

call dfftw_plan_dft_r2c_1d_(planf, M, freq, freq_k, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

do i = 1,M/2+1,1
freq_k(i) = 2.0d0*freq_k(i)/dfloat(M)
abs_freq_k(i) = abs(freq_k(i))
sort(i) = i
write(30,*) 2.0d0*pi*dfloat(i)/(dfloat(int(1.0d0/dt))),abs(freq_k(i))
enddo

call exchange (abs_freq_k, sort, M/2+1)

write(*,*) 2.0d0*pi*dfloat(sort(1))/(k*dfloat(int(1.0d0/dt)))

contains
!===================================================================================

SUBROUTINE derive(N, Nh, pi, time, u, b, du_dt, db_dt)

implicit none

INTEGER ( kind = 4 ) N, Nh, i
REAL ( kind = 8 ) pi, time, dt, rho0, B0, Re, Rm, dx, L, k
REAL ( kind = 8 ) u(N), u_dum(N), du_dx(N), d2u_dx2(N), du_dt(N)
REAL ( kind = 8 ) b(N), b_dum(N), db_dx(N), d2b_dx2(N), db_dt(N)
COMPLEX ( kind = 8 ) uk(Nh), ik_uk(Nh), k2_uk(Nh)
COMPLEX ( kind = 8 ) bk(Nh), ik_bk(Nh), k2_bk(Nh)

common/comm/L, dx, dt, rho0, B0, Re, Rm

do i = 1, N
  u_dum(i) = u(i)
  b_dum(i) = b(i)
enddo

call dfftw_plan_dft_r2c_1d_(planf, N, u_dum, uk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

call dfftw_plan_dft_r2c_1d_(planf, N, b_dum, bk, FFTW_ESTIMATE)
call dfftw_execute_(planf)
call dfftw_destroy_plan_(planf)

do i = 1, Nh
   k = 2.0d0 * pi * dfloat(i-1) / L
   ik_uk(i) = (0.0d0,1.0d0) * k * uk(i)
   ik_bk(i) = (0.0d0,1.0d0) * k * bk(i)
   k2_uk(i) = - k * k * uk(i)
   k2_bk(i) = - k * k * bk(i)
enddo

call dfftw_plan_dft_c2r_1d_(planb, N, ik_uk, du_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_1d_(planb, N, ik_bk, db_dx, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_1d_(planb, N, k2_uk, d2u_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

call dfftw_plan_dft_c2r_1d_(planb, N, k2_bk, d2b_dx2, FFTW_ESTIMATE)
call dfftw_execute_(planb)
call dfftw_destroy_plan_(planb)

do i = 1, N
  du_dx(i) = du_dx(i) / dfloat(N)
  db_dx(i) = db_dx(i) / dfloat(N)
  d2u_dx2(i) = d2u_dx2(i) / dfloat(N)
  d2b_dx2(i) = d2b_dx2(i) / dfloat(N)
enddo

DO i = 1, N
  du_dt(i) = - B0 * db_dx(i) / rho0 + d2u_dx2(i) / (rho0 * Re)
  db_dt(i) = - B0 * du_dx(i)        + d2b_dx2(i) / Rm
enddo

return
end

!====================================================================================

SUBROUTINE rk4(N, Nh, pi, time, u, b, du_dt, db_dt)
implicit none

INTEGER ( kind = 4 ) N, Nh, i
REAL ( kind = 8 ) pi, time, dt, rho0, B0, Re, Rm, dx, L
REAL ( kind = 8 ) u(N),du_dt(N),k1_u(N),k2_u(N),k3_u(N),k4_u(N),dum_u(N)
REAL ( kind = 8 ) b(N),db_dt(N),k1_b(N),k2_b(N),k3_b(N),k4_b(N),dum_b(N)

common/comm/L, dx, dt, rho0, B0, Re, Rm

!PRINT*, "RKIdx = 1"

do i = 1, N
  k1_u(i) = du_dt(i)
  dum_u(i) = u(i) + k1_u(i) * dt / 2.0d0
enddo

do i = 1, N
  k1_b(i) = db_dt(i)
  dum_b(i) = b(i) + k1_b(i) * dt / 2.0d0
enddo

CALL derive(N, Nh, pi, time+dt/2.0d0, dum_u, dum_b, du_dt, db_dt)

!PRINT*, "RKIdx = 2"

do i = 1, N
  k2_u(i) = du_dt(i)
  dum_u(i) = u(i) + k2_u(i) * dt/2.0d0
enddo

do i = 1, N
  k2_b(i) = db_dt(i)
  dum_b(i) = b(i) + k2_b(i) * dt/2.0d0
enddo

CALL derive(N, Nh, pi, time+dt/2.0d0, dum_u, dum_b, du_dt, db_dt)

!PRINT*, "RKIdx = 3"

do i = 1, N
  k3_u(i) = du_dt(i)
  dum_u(i) = u(i) + k3_u(i) * dt
enddo

do i = 1, N
  k3_b(i) = db_dt(i)
  dum_b(i) = b(i) + k3_b(i) * dt
enddo

CALL derive(N, Nh, pi, time+dt, dum_u, dum_b, du_dt, db_dt)

!PRINT*, "RKIdx = 4"

do i = 1, N
  k4_u(i) = du_dt(i)
  u(i) = u(i) + dt/6.0d0*(k1_u(i) + 2.0d0*k2_u(i) + 2.0d0*k3_u(i) + k4_u(i))
enddo

do i = 1, N
  k4_b(i) = db_dt(i)
  b(i) = b(i) + dt/6.0d0*(k1_b(i) + 2.0d0*k2_b(i) + 2.0d0*k3_b(i) + k4_b(i))
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

end program inhomogeneous_wave_equation

! To run the program : gfortran -w -I/usr/local/include -L/usr/local/lib <filename.f95> -lfftw3 -lm; ./a.out
! It uses Adam Bashforth Solver
program Berezin_Karpman
implicit none

include "fftw3.f"

integer, parameter :: N = 1024
integer, parameter :: Nh = N/2+1
double precision, parameter :: pi=3.14159265358979323846d0

double precision :: x,rho,vx,vy,vz,By,Bz
dimension x(N),rho(N),vx(N),vy(N),vz(N),By(N),Bz(N)
double precision :: dvx,dvy,dvz,dBy,dBz
dimension dvx(N),dvy(N),dvz(N),dBy(N),dBz(N)
! double precision :: ddvx,ddvy,ddvz,ddBy,ddBz
! dimension ddvx(N),ddvy(N),ddvz(N),ddBy(N),ddBz(N)
double precision :: rhovx,B2,curly,curlz
dimension rhovx(N),B2(N),curly(N),curlz(N)
double precision :: drhovx,dB2,dcurly,dcurlz
dimension drhovx(N),dB2(N),dcurly(N),dcurlz(N)
double precision :: sTdBy,sTdBz,dsTdBy,dsTdBz
dimension sTdBy(N),sTdBz(N),dsTdBy(N),dsTdBz(N)
double complex :: vxk,vyk,vzk,Byk,Bzk
dimension vxk(Nh),vyk(Nh),vzk(Nh),Byk(Nh),Bzk(Nh)
double complex :: dvxk,dvyk,dvzk,dByk,dBzk
dimension dvxk(Nh),dvyk(Nh),dvzk(Nh),dByk(Nh),dBzk(Nh)
! double complex :: ddvxk,ddvyk,ddvzk,ddByk,ddBzk
! dimension ddvxk(Nh),ddvyk(Nh),ddvzk(Nh),ddByk(Nh),ddBzk(Nh)
double complex :: rhovxk,B2k,curlyk,curlzk
dimension rhovxk(Nh),B2k(Nh),curlyk(Nh),curlzk(Nh)
double complex :: drhovxk,dB2k,dcurlyk,dcurlzk
dimension drhovxk(Nh),dB2k(Nh),dcurlyk(Nh),dcurlzk(Nh)
double complex :: sTdByk,sTdBzk,dsTdByk,dsTdBzk
dimension sTdByk(Nh),sTdBzk(Nh),dsTdByk(Nh),dsTdBzk(Nh)
integer*8 :: plan
real*8 :: L

integer :: i,t
real*8 :: dx,rho0,vx0,vy0,vz0,By0,Bz0,k,time,time_min,time_max,dt,nu
real*8 :: sinTheta,cosTheta
real*8, dimension (N) :: force_rho_n,force_vx_n,force_vy_n,force_vz_n,force_By_n,force_Bz_n
real*8, dimension (N) :: force_rho_o,force_vx_o,force_vy_o,force_vz_o,force_By_o,force_Bz_o
external :: rho0,vx0,vy0,vz0,By0,Bz0

L = 2.0d0*pi
dx = L/dfloat(N)

time_min = 0.00d0
time_max = 0.001d0
dt = 0.00010d0

nu = 0.01d0

sinTheta = 0.50d0
cosTheta = dsqrt(1.0d0 - sinTheta*sinTheta)

do i=1,N,1
  x(i)  = dfloat(i-1) * dx
  rho(i)= rho0(x(i))
  vx(i) = vx0(x(i))
  vy(i) = vy0(x(i))
  vz(i) = vz0(x(i))
  By(i) = By0(x(i))
  Bz(i) = Bz0(x(i))
  write (100,*) 0, x(i), rho(i), vx(i), vy(i), vz(i), By(i), Bz(i)
enddo

do time = time_min, time_max, dt

t = nint(time/dt) - int(time_min/dt)

do i=1,N,1
  rhovx(i) = (1.0d0 + rho(i)) * vx(i)
  B2(i)    = ( cosTheta + Bz(i) ) * ( cosTheta + Bz(i) ) + By(i) * By(i)
  curly(i) = vx(i)*By(i) - vy(i)*sinTheta
  curlz(i) = vz(i)*sinTheta - vx(i)*(cosTheta+Bz(i))
enddo

call dfftw_plan_dft_r2c_1d(plan,N,vx,vxk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,vx,vxk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,vy,vyk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,vy,vyk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,vz,vzk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,vz,vzk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,By,Byk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,By,Byk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,Bz,Bzk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,Bz,Bzk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,rhovx,rhovxk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,rhovx,rhovxk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,B2,B2k,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,B2,B2k)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,curly,curlyk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,curly,curlyk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,curlz,curlzk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,curlz,curlzk)
call dfftw_destroy_plan(plan)

!do i=1,Nh
!  if (i .ge. N/3+1) then     ! De - Aliazing Technique...
!    uk(i) = (0.0d0,0.0d0)
!    vk(i) = (0.0d0,0.0d0)
!  endif
!enddo

do i=1,Nh,1
  k = 2.0d0 * pi * dfloat(i-1) / L
  dvxk(i) = (0.d0,1.d0) * k * vxk(i)
  dvyk(i) = (0.d0,1.d0) * k * vyk(i)
  dvzk(i) = (0.d0,1.d0) * k * vzk(i)
  dByk(i) = (0.d0,1.d0) * k * Byk(i)
  dBzk(i) = (0.d0,1.d0) * k * Bzk(i)
  drhovxk(i) = (0.d0,1.d0) * k * rhovxk(i)
  dB2k(i)    = (0.d0,1.d0) * k * B2k(i)
  dcurlyk(i) = (0.d0,1.d0) * k * curlyk(i)
  dcurlzk(i) = (0.d0,1.d0) * k * curlzk(i)
  ! ddvxk(i) = - k * k * vxk(i)
  ! ddvyk(i) = - k * k * vyk(i)
  ! ddvzk(i) = - k * k * vzk(i)
  ! ddByk(i) = - k * k * Byk(i)
  ! ddBzk(i) = - k * k * Bzk(i)
enddo

call dfftw_plan_dft_c2r_1d(plan,N,dvxk,dvx,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dvxk,dvx)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dvyk,dvy,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dvyk,dvy)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dvzk,dvz,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dvzk,dvz)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dByk,dBy,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dByk,dBy)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dBzk,dBz,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dBzk,dBz)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,drhovxk,drhovx,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,drhovxk,drhovx)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dB2k,dB2,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dB2k,dB2)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dcurlyk,dcurly,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dcurlyk,dcurly)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dcurlzk,dcurlz,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dcurlzk,dcurlz)
call dfftw_destroy_plan(plan)

! call dfftw_plan_dft_c2r_1d(plan,N,ddvxk,ddvx,FFTW_ESTIMATE)
! call dfftw_execute_dft_c2r(plan,ddvxk,ddvx)
! call dfftw_destroy_plan(plan)
!
! call dfftw_plan_dft_c2r_1d(plan,N,ddvyk,ddvy,FFTW_ESTIMATE)
! call dfftw_execute_dft_c2r(plan,ddvyk,ddvy)
! call dfftw_destroy_plan(plan)
!
! call dfftw_plan_dft_c2r_1d(plan,N,ddvzk,ddvz,FFTW_ESTIMATE)
! call dfftw_execute_dft_c2r(plan,ddvzk,ddvz)
! call dfftw_destroy_plan(plan)
!
! call dfftw_plan_dft_c2r_1d(plan,N,ddByk,ddBy,FFTW_ESTIMATE)
! call dfftw_execute_dft_c2r(plan,ddByk,ddBy)
! call dfftw_destroy_plan(plan)
!
! call dfftw_plan_dft_c2r_1d(plan,N,ddBzk,ddBz,FFTW_ESTIMATE)
! call dfftw_execute_dft_c2r(plan,ddBzk,ddBz)
! call dfftw_destroy_plan(plan)

do i=1,N,1
  dvx(i) = dvx(i) / dfloat(N)
  dvy(i) = dvy(i) / dfloat(N)
  dvz(i) = dvz(i) / dfloat(N)
  dBy(i) = dBy(i) / dfloat(N)
  dBz(i) = dBz(i) / dfloat(N)
  drhovx(i) = drhovx(i) / dfloat(N)
  dB2(i)    = dB2(i)    / dfloat(N)
  dcurly(i) = dcurly(i) / dfloat(N)
  dcurlz(i) = dcurlz(i) / dfloat(N)
  ! ddvx(i) = ddvx(i) / dfloat(N)
  ! ddvy(i) = ddvy(i) / dfloat(N)
  ! ddvz(i) = ddvz(i) / dfloat(N)
  ! ddBy(i) = ddBy(i) / dfloat(N)
  ! ddBz(i) = ddBz(i) / dfloat(N)
enddo

do i=1,N,1
  sTdBy(i) = sinTheta * dBy(i) / (1.0d0 + rho(i))
  sTdBz(i) = sinTheta * dBz(i) / (1.0d0 + rho(i))
enddo

call dfftw_plan_dft_r2c_1d(plan,N,sTdBy,sTdByk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,sTdBy,sTdByk)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_1d(plan,N,sTdBz,sTdBzk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,sTdBz,sTdBzk)
call dfftw_destroy_plan(plan)

do i=1,Nh,1
  k = 2.0d0 * pi * dfloat(i-1) / L
  dsTdByk(i) = (0.d0,1.d0) * k * sTdByk(i)
  dsTdBzk(i) = (0.d0,1.d0) * k * sTdBzk(i)
enddo

call dfftw_plan_dft_c2r_1d(plan,N,dsTdByk,dsTdBy,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dsTdByk,dsTdBy)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan,N,dsTdBzk,dsTdBz,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dsTdBzk,dsTdBz)
call dfftw_destroy_plan(plan)

do i=1,N,1
  dsTdBy(i) = dsTdBy(i) / dfloat(N)
  dsTdBz(i) = dsTdBz(i) / dfloat(N)
enddo

do i = 1,N,1
  force_rho_n(i)= - drhovx(i)
  force_vx_n(i) = - vx(i) * dvx(i) - dB2(i) / ( 2.0d0 * (1.0d0 + rho(i)) )
  force_vy_n(i) = - vx(i) * dvy(i) + sinTheta * dBy(i) / (1.0d0 + rho(i))
  force_vz_n(i) = - vx(i) * dvz(i) + sinTheta * dBz(i) / (1.0d0 + rho(i))
  force_By_n(i) = - dcurly(i) + dsTdBz(i)
  force_Bz_n(i) =   dcurlz(i) + dsTdBy(i)
enddo

if (t == 0) then ! Taylor
  do i = 1,N,1
    rho(i)= rho(i) + dt * force_rho_n(i)
    vx(i) = vx(i)  + dt * force_vx_n(i)
    vy(i) = vy(i)  + dt * force_vy_n(i)
    vz(i) = vz(i)  + dt * force_vz_n(i)
    By(i) = By(i)  + dt * force_By_n(i)
    Bz(i) = Bz(i)  + dt * force_Bz_n(i)
  enddo
else ! Adams - Bashforth
  do i = 1,N,1
    rho(i)= rho(i)+ dt * ( (3.0/2.0) * force_rho_n(i)- (1.0/2.0) * force_rho_o(i) )
    vx(i) = vx(i) + dt * ( (3.0/2.0) * force_vx_n(i) - (1.0/2.0) * force_vx_o(i) )
    vy(i) = vy(i) + dt * ( (3.0/2.0) * force_vy_n(i) - (1.0/2.0) * force_vy_o(i) )
    vz(i) = vz(i) + dt * ( (3.0/2.0) * force_vz_n(i) - (1.0/2.0) * force_vz_o(i) )
    By(i) = By(i) + dt * ( (3.0/2.0) * force_By_n(i) - (1.0/2.0) * force_By_o(i) )
    Bz(i) = Bz(i) + dt * ( (3.0/2.0) * force_Bz_n(i) - (1.0/2.0) * force_Bz_o(i) )
  enddo
endif

! if (t .ge. 1000 .and. mod(t,500) == 0) then
  do i = 1,N,1
    write (t+100,*) t, x(i), rho(i), vx(i), vy(i), vz(i), By(i), Bz(i)
  enddo
! endif

do i = 1,N,1
  force_rho_o(i)= force_rho_n(i)
  force_vx_o(i) = force_vx_n(i)
  force_vy_o(i) = force_vy_n(i)
  force_vz_o(i) = force_vz_n(i)
  force_By_o(i) = force_By_n(i)
  force_Bz_o(i) = force_Bz_n(i)
enddo

enddo ! time

end program Berezin_Karpman

!==============================================================

function rho0(x)
implicit none
real*8 rho0,pi,x,L
pi = 3.14159265358979323846d0
L = 2.0d0*pi
rho0 = 1.0d0 !sin(1.0d0*x)
return
end

!================================================================

function vx0(x)
implicit none
real*8 vx0,pi,x,L
pi = 3.14159265358979323846d0
L = 2.0d0*pi
vx0 = 0.0d0 !sin(1.0d0*x)
return
end

!==============================================================

function vy0(x)
implicit none
real*8 vy0,pi,x,L
pi = 3.14159265358979323846d0
L = 2.0d0*pi
vy0 = 0.0d0 !sin(1.0d0*x)
return
end

!==============================================================

function vz0(x)
implicit none
real*8 vz0,pi,x,L
pi = 3.14159265358979323846d0
L = 2.0d0*pi
vz0 = 0.0d0 !sin(1.0d0*x)
return
end

!==============================================================

function By0(x)
implicit none
real*8 By0,pi,x,L
pi = 3.14159265358979323846d0
L = 2.0d0*pi
By0 = sin(1.0d0*x)
return
end

!==============================================================

function Bz0(x)
implicit none
real*8 Bz0,pi,x,L
pi = 3.14159265358979323846d0
L = 2.0d0*pi
Bz0 = sin(1.0d0*x)
return
end

!==============================================================

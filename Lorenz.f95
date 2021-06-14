program LorenzChimera
! Chimera states and intermittency in an ensemble of nonlocally coupled Lorenz systems
! I. A. Shepelev, G. I. Strelkova, and V. S. Anishchenko
! Chaos 28, 063119 (2018); https://doi.org/10.1063/1.5020009
! Important reference: Coherent traveling waves in nonlocally coupled chaotic systems,
! Volodymyr Dziubak,Yuri Maistrenko, and Eckehard Scholl
! PRE, 87, 032907 (2013)
implicit none

integer ( kind = 4 ), parameter :: Particle = 300
integer ( kind = 4 ), parameter :: N = 3 * Particle

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real*8 time,dt,xi,rho,beta,sigma,r
real*8 w(N),wp(N)
integer*8 i,j,time_max,time_start,time_write,timestep,t,h,th,P,pIndex

common/comm/xi,rho,beta,sigma,r,P

integer,parameter :: seed = 99999999
CALL srand(seed)

!================ File Names ==========================================

open(unit=10,file='Initial_Condition.dat',status='unknown')

!================ INITIAL CONDITION ===================================

 time_max   = 50
 time_write = 0
 timestep   = 1000
 time       = 0.0d0
 dt         = 0.001d0

 xi    = 10.0d0
 rho   = 28.0d0
 beta  = 8.0d0/3.0d0

! Choose parameters as region C of Figure 2 of Chaos, 28, 063119 (2018)
 sigma = 0.70d0
 r     = 0.090d0
 P     = int(r*N)

do i = 1, N, 3
  w(i)   = 1.0d0 * dcos(2.0d0*pi*dfloat(i)/dfloat(N) + 4.0d0*pi*(rand()-0.5))
  w(i+1) = 1.0d0 * dsin(2.0d0*pi*dfloat(i)/dfloat(N) + 4.0d0*pi*(rand()-0.5))
  w(i+2) = 1.0d0 * dsin(2.0d0*pi*dfloat(i)/dfloat(N) + 4.0d0*pi*(rand()-0.5))
  write(10,*) i,w(i),w(i+1),w(i+2)
enddo

  close(10)
!================ MAIN PROGRAM ========================================

do t = 1, time_max, 1
  print*, t
  do h = 1, timestep, 1
    th = t * timestep + h
    time = time + dt

    call derive(N,Particle,time,w,wp)
    call RK4(N,Particle,time,dt,w,wp)

    do i = 1, N, 3
      pIndex = (i-1)/3 + 1
      if (t .ge. time_write .and. mod(t,1) == 0 .and. i .le. 850) then ! Choose the i-limit depending on the stack limit of your machine.
        write(pIndex+100,*) time, pIndex, w(i), w(i+1), w(i+2)
      endif
    enddo

  enddo !h
enddo !t

contains

!================ DIFFERENTIAL EQUATION ===============================

subroutine derive(N,Particle,time,w,wp)
implicit none
real*8, intent (in) :: w(N)
real*8, intent (out) :: wp(N)
real*8 time,dt,xi,rho,beta,sigma,r
real*8 w_bar(N)
integer*4 N,Particle
integer*8 t,i,k,k0,k1,k2,P
common/comm/xi,rho,beta,sigma,r,P

do i = 1, N, 3
  w_bar(i)   = 0.0d0
  w_bar(i+1) = 0.0d0
  w_bar(i+2) = 0.0d0
enddo

do i = 1, N, 3
  do k = i - 3*P, i + 3*P, 3
    if (k < 1) then
      k0 = k + N
      k1 = k + 1 + N
      k2 = k + 2 + N
    elseif (k > N) then
      k0 = k - N
      k1 = k + 1 - N
      k2 = k + 2 - N
    else
      k0 = k
      k1 = k + 1
      k2 = k + 2
    endif
    w_bar(i)   = w_bar(i)   + w(k0)
    w_bar(i+1) = w_bar(i+1) + w(k1)
    w_bar(i+2) = w_bar(i+2) + w(k2)
  enddo
enddo


do i = 1, N, 3
  wp(i)   = xi * (w(i+1) - w(i)) - sigma * w(i) + sigma/(2.0d0*P) * w_bar(i)
  wp(i+1) = w(i) * (rho - w(i+2)) - w(i+1) - sigma * w(i+1) + sigma/(2.0d0*P) * w_bar(i+1)
  wp(i+2) = w(i) * w(i+1) - beta * w(i+2)
enddo

return
end subroutine derive

!=================== RK 4 subroutine ===================================

subroutine RK4(N,Particle,time,dt,w,wp)
implicit none
real*8, intent (inout) :: w(N)
real*8 time,dt,xi,rho,beta,sigma,r
real*8 wp(N),w_dum(N),w_new(N,4)
integer*4 N,Particle
integer*8 t,i,P
common/comm/xi,rho,beta,sigma,r,P

do i=1,N
  w_dum(i) = w(i)
enddo

  call derive(N,Particle,time,w_dum,wp)

do i=1,N
  w_new(i,1)=dt*wp(i)
  w_dum(i)=w(i)+w_new(i,1)/2.0d0
enddo

  call derive(N,Particle,time+dt/2.0,w_dum,wp)

do i=1,N
  w_new(i,2)=dt*wp(i)
  w_dum(i)=w(i)+w_new(i,2)/2.0d0
enddo

  call derive(N,Particle,time+dt/2.0,w_dum,wp)

do i=1,N
  w_new(i,3)=dt*wp(i)
  w_dum(i)=w(i)+w_new(i,3)
enddo

  call derive(N,Particle,time+dt,w_dum,wp)

do i=1,N
  w_new(i,4)=dt*wp(i)
  w(i)=w(i)+1.0d0/6.0d0*(w_new(i,1)+2.0d0*(w_new(i,2)+w_new(i,3))+w_new(i,4))
enddo

return
end subroutine RK4

!================= End Program =========================================

end program LorenzChimera

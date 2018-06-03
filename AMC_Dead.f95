program Chimera
! Rupak Mukherjee, Abhijit Sen, arXiv: 1710.09608, Chaos: 28, 053109 (2018).
implicit none

include "fftw3.f"

integer ( kind = 4 ), parameter :: Particle = 201
integer ( kind = 4 ), parameter :: N = 2 * Particle
integer ( kind = 4 ), parameter :: M = 10000

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real*8 time,dt,C1,C2,kappa,p,w_real_bar,w_img_bar,R,RA,RI,Order,phi
real*8 w(N),wp(N),w_img_sum(0:N), w_real_sum(-1:N)
real*8  w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A
integer*8 i,j,time_max,time_start,time_write,timestep,t,h,th,Dead
real*8 plan ! FFTW

double precision :: omega,abs_omega_k
dimension omega(M),abs_omega_k(M/2+1)
double complex :: omegak
dimension omegak(M/2+1)
integer*8 sort(M/2+1)

common/comm/kappa,C1,C2,Dead

!================ File Names ==========================================

open(unit=10,file='Initial_Condition.dat',status='unknown')
open(unit=20,file='chimera.dat',status='unknown')
open(unit=30,file='FFT.dat',status='unknown')
open(unit=40,file='Parameter.dat',status='unknown')

!================ INITIAL CONDITION ===================================

 time_max = 4500
 time_start = 4400
 time_write = 4495 
 timestep = 100
 time = 0.0d0
 dt = 0.01d0
 
 kappa = 0.7    
 C1 = -1.0      
 C2 = 2.0d0
 p = 0.2   
 Dead = int(p*N)
 R = 0.0d0
 RA = 0.0d0
 RI = 0.0d0
 
do i=1,N,1
  if (mod(i,2) .eq. 0) then
  w(i) = 1.0d0 * dcos(2.0d0*pi*dfloat(i)/dfloat(N))
  w(i-1) = 1.0d0 * dsin(2.0d0*pi*dfloat(i)/dfloat(N))
  write(10,*) i,w(i),w(i-1)
  endif
enddo

 close(10)
!================ MAIN PROGRAM ========================================
     
do t=1,time_max,1
  do h=1,timestep,1
  th = t*timestep+h
  time=time+dt
  
  call derive(N,Particle,time,w,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A) 
  call RK4(N,Particle,time,dt,w,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)
  
    if (t .gt. time_start .and. mod(t,1) == 0 .and. t .le. time_max) then 
    write(20,*) time, w_real_bar, w_real_bar_A, w_real_bar_I, dsqrt((w_real_bar)**2+(w_img_bar)**2),&
                dsqrt((w_real_bar_A)**2+(w_img_bar_A)**2), dsqrt((w_real_bar_I)**2+(w_img_bar_I)**2)
    R = R + dsqrt((w_real_bar)**2+(w_img_bar)**2)/dfloat((time_max-time_start)*timestep)
    RA = RA + dsqrt((w_real_bar_A)**2+(w_img_bar_A)**2)/dfloat((time_max-time_start)*timestep) 
    RI = RI + dsqrt((w_real_bar_I)**2+(w_img_bar_I)**2)/dfloat((time_max-time_start)*timestep)
    j = th - time_start * timestep - timestep
    omega(j) = w_real_bar
    endif
    
    if (t .ge. time_write .and. mod(t,1) == 0 .and. t .le. time_max) then
      do i=1,N,1
        if (mod(i,2) .eq. 0) then
        Order = dsqrt(w(i)**2 + w(i+1)**2)
	      if (w(i) .ge. 0.0d0 .and. w(i-1) .ge. 0.0d0) then
	      phi = datan(abs(w(i-1)/w(i)))
	      elseif (w(i) .lt. 0.0d0 .and. w(i-1) .ge. 0.0d0) then
	      phi = pi - datan(abs(w(i-1)/w(i)))
	      elseif (w(i) .lt. 0.0d0 .and. w(i-1) .lt. 0.0d0) then
	      phi = pi + datan(abs(w(i-1)/w(i)))
	      elseif (w(i) .ge. 0.0d0 .and. w(i-1) .lt. 0.0d0) then
	      phi = 2.0d0*pi - datan(abs(w(i-1)/w(i)))
	      endif
        write(th,*)time,i/2,w(i),w(i-1),Order,phi
        endif 
      enddo ! i
    endif
     
    if (t .ge. time_write .and. mod(t,1) == 0 .and. t .le. time_max) then
    close (th)
    endif
    
  enddo !h
enddo !t

  call dfftw_plan_dft_r2c_1d(plan,M,omega,omegak,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan,omega,omegak)
  call dfftw_destroy_plan(plan)

do i = 1,M/2+1,1
omegak(i) = 2.0d0*omegak(i)/dfloat(M)
abs_omega_k(i) = abs(omegak(i))
sort(i) = i
write(30,*) 2.0d0*pi*dfloat(i)/(dfloat(timestep)),abs(omegak(i))
enddo

call exchange (abs_omega_k, sort, M/2+1) 

write(40,*) kappa, C1, R, RA, RI, 2.0d0*pi*dfloat(sort(1))/(dfloat(timestep))

 close (20)
 close (30)
 close (40)

contains

!================ DIFFERENTIAL EQUATION ===============================

subroutine derive(N,Particle,time,w,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)
implicit none
real*8, intent (in) :: w(N)
real*8, intent (out) :: wp(N)
real*8 time,dt,C1,C2,kappa,w_real_bar,w_img_bar
real*8  r(N),w_img_sum(0:N), w_real_sum(-1:N)
real*8  w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A
integer*4 N,Particle
integer*8 t,i,Dead
common/comm/kappa,C1,C2,Dead
     
w_real_sum(-1) = 0.0d0
w_img_sum(0) = 0.0d0
w_img_bar_I = 0.0d0
w_img_bar_A = 0.0d0
w_real_bar_I = 0.0d0
w_real_bar_A = 0.0d0

do i=1,N,1
  if (mod(i,2) .eq. 0) then 
  w_img_sum(i)=w_img_sum(i-2)+w(i)
    if (i .le. Dead) then
    w_img_bar_I = w_img_bar_I + w(i)
    elseif (i .gt. Dead) then
    w_img_bar_A = w_img_bar_A + w(i)
    endif
  elseif (mod(i,2) .eq. 1) then
  w_real_sum(i)=w_real_sum(i-2)+w(i)
    if (i .le. Dead) then
    w_real_bar_I = w_real_bar_I + w(i)
    elseif (i .gt. Dead) then
    w_real_bar_A = w_real_bar_A + w(i)
    endif
  endif
enddo

w_real_bar= w_real_sum(N-1)/Particle
w_img_bar= w_img_sum(N)/Particle
w_real_bar_I = w_real_bar_I/(Dead/2)
w_real_bar_A = w_real_bar_A/((N-Dead)/2)
w_img_bar_I = w_img_bar_I/(Dead/2)
w_img_bar_A = w_img_bar_A/((N-Dead)/2)

do i=1,N, 2

  if (i .le. Dead) then
  ! Inactive/Dead
  r(i) = -1.0d0 
  else
  ! Active
  r(i) = +1.0d0 
  endif
  
wp(i)   = r(i)*w(i) - (w(i)*w(i) + w(i+1)*w(i+1)) * (w(i) - C2*w(i+1)) &
          + kappa * ( (w_real_bar - w(i)) - C1 * (w_img_bar - w(i+1)) )
          
wp(i+1) = r(i)*w(i+1) - (w(i)*w(i) + w(i+1)*w(i+1)) * (w(i+1) + C2*w(i)) &
          + kappa * ( C1 * (w_real_bar - w(i)) + (w_img_bar - w(i+1)) )
enddo
   
return
end subroutine derive

!=================== RK 4 subroutine ===================================

subroutine RK4(N,Particle,time,dt,w,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)	
implicit none
real*8, intent (inout) :: w(N)
real*8 time,dt,C1,C2,kappa,w_real_bar,w_img_bar
real*8 wp(N),w_dum(N),w_new(N,4)
real*8  w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A
integer*4 N,Particle
integer*8 t,i,Dead
common/comm/kappa,C1,C2,Dead

do i=1,N
w_dum(i) = w(i)
enddo
	
call derive(N,Particle,time,w_dum,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)
do i=1,N
w_new(i,1)=dt*wp(i)
w_dum(i)=w(i)+w_new(i,1)/2.0d0	
enddo	
	
call derive(N,Particle,time+dt/2.0,w_dum,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)
do i=1,N
w_new(i,2)=dt*wp(i)
w_dum(i)=w(i)+w_new(i,2)/2.0d0
enddo

call derive(N,Particle,time+dt/2.0,w_dum,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)
do i=1,N
w_new(i,3)=dt*wp(i)
w_dum(i)=w(i)+w_new(i,3)
enddo	

call derive(N,Particle,time+dt,w_dum,wp,w_real_bar,w_img_bar,w_img_bar_I,w_img_bar_A,w_real_bar_I,w_real_bar_A)
do i=1,N
w_new(i,4)=dt*wp(i)
w(i)=w(i)+1.0d0/6.0d0*(w_new(i,1)+2.0d0*(w_new(i,2)+w_new(i,3))+w_new(i,4))
enddo

return
end subroutine RK4

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

!================= End Program =========================================
    
end program Chimera    	  

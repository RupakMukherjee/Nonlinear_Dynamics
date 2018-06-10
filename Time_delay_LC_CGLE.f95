program LCCGLE
! Solution of Locally Coupled Complex Ginzburg Landau Equation with Time delay.
implicit none

integer ( kind = 4 ), parameter :: Particle = 500
integer ( kind = 4 ), parameter :: N = 2 * Particle
integer ( kind = 4 ), parameter :: timestep = 100
integer ( kind = 4 ), parameter :: tau = 1*timestep + 1
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real*8 time,dt,C1,C2,kappa,omega0,R,phi,ka
real*8 w(N),wp(N),w_s(tau,N),w_tau(N)
integer*4 i,time_max,t,h,th,s

common/comm/kappa,C1,C2,omega0

open(unit=10,file='Initial_Condition.dat',status='unknown')
open(unit=20,file='chimera.dat',status='unknown')

!================ INITIAL CONDITION ===================================

 time_max = 100 ! 4500
 time = 0.0d0
 dt = 0.01d0
 
 kappa = 10.0d0
 C1 = 2.0d0   
 C2 = -1.0d0
 omega0 = 10.0d0

 ka = 4.0d0*pi/dfloat(Particle)
 
do i=1,N,2
  if (mod(i,2) == 1) then
  w(i) = 1.0d0 * dsin(ka*dfloat(i/2))
  w(i+1) = 1.0d0
  write(10,*) i/2,w(i),w(i+1)
  endif
enddo

  close(10)

do s = 1,tau,1
  do i = 1,N,1
  w_s(s,i) = w(i)
  enddo
enddo  

s = 1

!================ MAIN PROGRAM ========================================
     
do t = 1, time_max, 1
  do h = 1, timestep, 1
  th = t * timestep + h
  time = time + dt
  
    if (mod(th,tau) == 0) then
    s = 1
    else
    s = s + 1
    endif
    
    do i = 1,N
    w_tau(i) = w_s(s,i)
    enddo
    
  call derive (N,time,w,w_tau,wp) 
  call RK4 (N,time,dt,w,w_tau,wp)
  
    do i = 1,N
    w_s(s,i) = w(i)
    enddo
  
    if (t .gt. 95 .and. mod(t,1) == 0) then
      do i=1,N,1
        if (mod(i,2) == 1) then
        R = sqrt(w(i)**2 + w(i+1)**2)
	      if (w(i) .ge. 0.0d0 .and. w(i+1) .ge. 0.0d0) then
	      phi = datan(abs(w(i+1)/w(i)))
	      elseif (w(i) .lt. 0.0d0 .and. w(i+1) .ge. 0.0d0) then
	      phi = pi - datan(abs(w(i+1)/w(i)))
	      elseif (w(i) .lt. 0.0d0 .and. w(i+1) .lt. 0.0d0) then
	      phi = pi + datan(abs(w(i+1)/w(i)))
	      elseif (w(i) .ge. 0.0d0 .and. w(i+1) .lt. 0.0d0) then
	      phi = 2.0d0*pi - datan(abs(w(i+1)/w(i)))
	      endif
          if (mod(t,1) == 0) then
          write(th,*)time,i/2,w(i),w(i+1),R,phi
          endif
        write(20,*)time,i/2,w(i),w(i+1),R,phi
        endif 
      enddo
    endif
     
  enddo !h
  if (mod(t,1) == 0) then
  close(th)  
  endif
enddo !t

  close(20)
  
contains

!================ DIFFERENTIAL EQUATION ===============================

subroutine derive(N,time,w_new,w_tau_old,wp)
implicit none
real*8, intent (in) :: w_new(N),w_tau_old(N)
real*8, intent (out) :: wp(N)
real*8 time,dt,C1,C2,kappa,omega0,w2
real*8  w(N),w_tau(-1:N+2)
integer*4 N,th,i
common/comm/kappa,C1,C2,omega0

w_tau(-1)  = w_tau_old(N-1)
w_tau(0)   = w_tau_old(N)
w_tau(N+1) = w_tau_old(1)
w_tau(N+2) = w_tau_old(2)

do i = 1, N, 2
w(i) = w_new(i)
w(i+1) = w_new(i+1)
w2 = w(i)*w(i) + w(i+1)*w(i+1)

  if (i == 1) then
    w_tau(i+2) = w_tau_old(i+2)
    w_tau(i+3) = w_tau_old(i+3)
  elseif (i == N-1) then
    w_tau(i-2) = w_tau_old(i-2)
    w_tau(i-1) = w_tau_old(i-1)
  elseif (i /= 1 .and. i /= N-1) then
    w_tau(i-2) = w_tau_old(i-2)
    w_tau(i-1) = w_tau_old(i-1)
    w_tau(i+2) = w_tau_old(i+2)
    w_tau(i+3) = w_tau_old(i+3)
  endif

wp(i)   = w(i) - omega0 * w(i+1) - w2 * ( w(i) - C2*w(i+1) ) &
               + kappa * (w_tau(i+2)+w_tau(i-2)-2.0d0*w(i)) - kappa * C1 * (w_tau(i+3)+w_tau(i-1)-2.0d0*w(i+1))
       
wp(i+1) = w(i+1) + omega0 * w(i) - w2 * ( w(i+1) + C2*w(i) ) &
               + kappa * (w_tau(i+3)+w_tau(i-1)-2.0d0*w(i+1)) + kappa * C1 * (w_tau(i+2)+w_tau(i-2)-2.0d0*w(i)) 

enddo
   
return
end subroutine derive

!=================== RK 4 subroutine ===================================

subroutine RK4(N,time,dt,w,w_tau,wp)	
implicit none
real*8, intent (inout) :: w(N)
real*8 time,dt,C1,C2,kappa,omega0
real*8 wp(N),w_dum(N),w_new(N,4),w_tau(N)
integer*4 N,th,i
common/comm/kappa,C1,C2,omega0

do i=1,N
w_dum(i) = w(i)
enddo
	
call derive(N,time,w_dum,w_tau,wp)
do i=1,N
w_new(i,1)=dt*wp(i)
w_dum(i)=w(i)+w_new(i,1)/2.0d0	
enddo	
	
call derive(N,time+dt/2.0,w_dum,w_tau,wp)
do i=1,N
w_new(i,2)=dt*wp(i)
w_dum(i)=w(i)+w_new(i,2)/2.0d0
enddo

call derive(N,time+dt/2.0,w_dum,w_tau,wp)
do i=1,N
w_new(i,3)=dt*wp(i)
w_dum(i)=w(i)+w_new(i,3)
enddo	

call derive(N,time+dt,w_dum,w_tau,wp)
do i=1,N
w_new(i,4)=dt*wp(i)
w(i)=w(i)+1.0d0/6.0d0*(w_new(i,1)+2.0d0*(w_new(i,2)+w_new(i,3))+w_new(i,4))
enddo

return
end subroutine RK4

!================= End Program =========================================
    
end program LCCGLE   	  

program weno

implicit none

integer, parameter :: N = 128
double precision, parameter :: pi=3.14159265358979323846d0

double precision :: dx,gamma_1,gamma_2,gamma_3,epslion,gamma_1R,gamma_2R,gamma_3R,time,time_min,time_max,dt,L

double precision :: x, u,uhp_1,uhp_2,uhp_3,beta_1,beta_2,beta_3,w1_bar,w2_bar,w3_bar,w_d,w1,w2,w3,uhp,u_bar,u_1,u_2,uhR,u_bar_dum
double precision ::uhR_1,uhR_2,uhR_3,beta_1R,beta_2R,beta_3R,w1R_bar,w2R_bar,w3R_bar,wR_d,w1R,w2R,w3R,RHS
double precision :: a, du_dt, d2u_dt2, du_dx, d2u_dx2,a_bar,a_bar_dum,a_1,a_2
dimension x(N),u(-1:N+2),u_bar(-1:N+2),u_bar_dum(-1:N+2),u_1(-1:N+2),u_2(-1:N+2),uhp(0:N),uhR(0:N)
dimension uhp_1(N),uhp_2(N),uhp_3(N),beta_1(N),beta_2(N),beta_3(N)
dimension w1_bar(N),w2_bar(N),w3_bar(N),w_d(N),w1(N),w2(N),w3(N)
dimension uhR_1(N),uhR_2(N),uhR_3(N),beta_1R(N),beta_2R(N),beta_3R(N)
dimension w1R_bar(N),w2R_bar(N),w3R_bar(N),wR_d(N),w1R(N),w2R(N),w3R(N),RHS(N)
dimension a(N), du_dt(N), d2u_dt2(N), du_dx(-1:N+2), d2u_dx2(N),a_bar(N),a_bar_dum(N),a_1(N),a_2(N)

integer :: i
integer*8 t

!====================================================

L = 2.0d0*pi
dx = L/dfloat(N-1)
   
time_min = 0.0d0
time_max = 1.50d0
dt = 0.0010d0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do i = 1,N
  x(i) = dfloat(i-1)*dx
  u(i) = dsin(x(i)) 
enddo

u(-1) = u(N-2)
u(0) = u(N-1)
u(N+1) = u(2)
u(N+2) = u(3)

gamma_1 = 1.0d0/16.0d0
gamma_2 = 5.0d0/8.0d0
gamma_3 = 5.0d0/16.0d0

epslion = 1.0d0/10.0d0**6

do i = 1,N
  uhp_1(i) = (3.0d0/8.0d0)*u(i-2) - (5.0d0/4.0d0)*u(i-1) + (15.0d0/8.0d0)*u(i)
  uhp_2(i) = (-1.0d0/8.0d0)*u(i-1) + (3.0d0/4.0d0)*u(i) + (3.0d0/8.0d0)*u(i+1)
  uhp_3(i) = (3.0d0/8.0d0)*u(i) + (3.0d0/4.0d0)*u(i+1) - (1.0d0/8.0d0)*u(i+2)

  beta_1(i) = (1.0d0/3.0d0)*(4.0d0*u(i-2)*u(i-2) - 19.0d0*u(i-2)*u(i-1) + 25.0d0*u(i-1)*u(i-1) &
                           + 11.0d0*u(i-2)*u(i) - 31.0d0*u(i-1)*u(i) + 10.0d0*u(i)*u(i))
  beta_2(i) = (1.0d0/3.0d0)*(4.0d0*u(i-1)*u(i-1) - 13.0d0*u(i-1)*u(i) + 13.0d0*u(i)*u(i) &
                           + 5.0d0*u(i-1)*u(i+1) - 13.0d0*u(i)*u(i+1) + 4.0d0*u(i+1)*u(i+1))
  beta_3(i) = (1.0d0/3.0d0)*(10.0d0*u(i)*u(i) - 31.0d0*u(i)*u(i+1) + 25.0d0*u(i+1)*u(i+1) &
                           + 11.0d0*u(i)*u(i+2) - 19.0d0*u(i+1)*u(i+2) + 4.0d0*u(i+2)*u(i+2))

  w1_bar(i) = (gamma_1)/(epslion+beta_1(i))**2
  w2_bar(i) = (gamma_2)/(epslion+beta_2(i))**2
  w3_bar(i) = (gamma_3)/(epslion+beta_3(i))**2

  w_d(i) = w1_bar(i) + w2_bar(i) + w3_bar(i)

  w1(i) = w1_bar(i)/w_d(i)
  w2(i) = w2_bar(i)/w_d(i)
  w3(i) = w3_bar(i)/w_d(i)

  uhp(i) = w1(i)*uhp_1(i) + w2(i)*uhp_2(i) + w3(i)*uhp_3(i)
enddo

uhp(0) = uhp(N)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do i = 1,N
  x(i) = dfloat(i-1)*dx  
  u_bar(i) = dsin(x(i)) !(uhp(i) + uhp(i-1)) / dx
  a(i) = dcos(x(i))
  du_dt(i) = a(i)
  write(100,*) x(i), u_bar(i)
enddo

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

  call weno_reconstruction (N,dx,u_bar,du_dx)
  call weno_reconstruction (N,dx,du_dx,d2u_dx2)
  
  do i = 1,N
    d2u_dt2(i) = d2u_dx2(i)
    u_bar_dum(i) = u_bar(i)
    u_1(i) = u_bar(i) + dt * du_dt(i)
    a_bar_dum(i) = a(i)
    a_1(i) = a(i) + dt * d2u_dt2(i)
    du_dt(i) = a_1(i)
  enddo
  
  call weno_reconstruction (N,dx,u_1,du_dx)
  call weno_reconstruction (N,dx,du_dx,d2u_dx2)
    
  do i = 1,N
    d2u_dt2(i) = d2u_dx2(i)
    u_2(i) = 3.0d0*u_bar_dum(i)/4.0d0 + u_1(i)/4.0d0 + dt * du_dt(i)/4.0d0
    a_2(i) = 3.0d0*a_bar_dum(i)/4.0d0 + a_1(i)/4.0d0 + dt * d2u_dt2(i)/4.0d0    
    du_dt(i) = a_2(i)
  enddo
  
  call weno_reconstruction (N,dx,u_2,du_dx)
  call weno_reconstruction (N,dx,du_dx,d2u_dx2)
  
  do i = 1,N
    d2u_dt2(i) = d2u_dx2(i)
    u_bar(i) = u_bar_dum(i)/3.0d0 + 2.0d0*u_2(i)/3.0d0 + dt * 2.0d0 * du_dt(i)/3.0d0
    a_bar(i) = a_bar_dum(i)/3.0d0 + 2.0d0*a_2(i)/3.0d0 + dt * 2.0d0 * d2u_dt2(i)/3.0d0    
    du_dt(i) = a_bar(i)
    write(t+101,*) x(i), u_bar(i)
  enddo
  
  close (t+101)
  
enddo ! time

end program weno

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine weno_reconstruction (N,dx,u_bar,RHS)
implicit none

double precision :: dx,gamma_1,gamma_2,gamma_3,epslion,gamma_1R,gamma_2R,gamma_3R,time,time_min,time_max,dt,L

double precision :: x, u,uhp_1,uhp_2,uhp_3,beta_1,beta_2,beta_3,w1_bar,w2_bar,w3_bar,w_d,w1,w2,w3,uhp,u_bar,uhR
double precision :: uhR_1,uhR_2,uhR_3,beta_1R,beta_2R,beta_3R,w1R_bar,w2R_bar,w3R_bar,wR_d,w1R,w2R,w3R,RHS
dimension x(N),u(-1:N+2),u_bar(-1:N+2),uhp(0:N),uhR(0:N)
dimension uhp_1(N),uhp_2(N),uhp_3(N),beta_1(N),beta_2(N),beta_3(N)
dimension w1_bar(N),w2_bar(N),w3_bar(N),w_d(N),w1(N),w2(N),w3(N)
dimension uhR_1(N),uhR_2(N),uhR_3(N),beta_1R(N),beta_2R(N),beta_3R(N)
dimension w1R_bar(N),w2R_bar(N),w3R_bar(N),wR_d(N),w1R(N),w2R(N),w3R(N),RHS(N)

integer :: i,N
integer*8 t

epslion = 1.0d0/10.0d0**6

gamma_1R = 1.0d0/10.0d0
gamma_2R = 3.0d0/5.0d0
gamma_3R = 3.0d0/10.0d0

  u_bar(-1) = u_bar(N-2)
  u_bar(0) = u_bar(N-1)
  u_bar(N+1) = u_bar(2)
  u_bar(N+2) = u_bar(3)

  do i = 1,N
    uhR_1(i) = (1.0d0/3.0d0)*u_bar(i-2) - (7.0d0/6.0d0)*u_bar(i-1) + (11.0d0/6.0d0)*u_bar(i)
    uhR_2(i) = (-1.0d0/6.0d0)*u_bar(i-1) + (5.0d0/6.0d0)*u_bar(i) + (1.0d0/3.0d0)*u_bar(i+1)
    uhR_3(i) = (1.0d0/3.0d0)*u_bar(i) + (5.0d0/6.0d0)*u_bar(i+1) - (1.0d0/6.0d0)*u_bar(i+2)

    beta_1R(i) = (13.0d0/12.0d0)*(u_bar(i-2) - 2.0d0*u_bar(i-1) + u_bar(i))**2 &
                +(1.0d0/4.0d0)*(u_bar(i-2) - 4.0d0*u_bar(i-1) + 3.0d0*u_bar(i))**2
    beta_2R(i) = (13.0d0/12.0d0)*(u_bar(i-1) - 2.0d0*u_bar(i) + u_bar(i+1))**2 &
                +(1.0d0/4.0d0)*(u_bar(i-1) - u_bar(i+1) )**2
    beta_3R(i) = (13.0d0/12.0d0)*(u_bar(i) - 2.0d0*u_bar(i+1) + u_bar(i+2))**2 &
                +(1.0d0/4.0d0)*(3.0d0*u_bar(i) - 4.0d0*u_bar(i+1) + u_bar(i+2))**2

    w1R_bar(i) = (gamma_1R)/(epslion+beta_1R(i))**2
    w2R_bar(i) = (gamma_2R)/(epslion+beta_2R(i))**2
    w3R_bar(i) = (gamma_3R)/(epslion+beta_3R(i))**2

    wR_d(i) = w1R_bar(i) + w2R_bar(i) + w3R_bar(i)
      
    w1R(i) = w1R_bar(i)/wR_d(i)
    w2R(i) = w2R_bar(i)/wR_d(i)
    w3R(i) = w3R_bar(i)/wR_d(i)
    
    uhR(i) = w1R(i)*uhR_1(i) + w2R(i)*uhR_2(i) + w3R(i)*uhR_3(i)
  enddo

  uhR(0) = uhR(N-1)

  do i = 1,N
    RHS(i) = - (uhR(i) - uhR(i-1)) / dx
  enddo

end subroutine weno_reconstruction

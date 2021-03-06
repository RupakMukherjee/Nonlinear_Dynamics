PROGRAM SIR
!
!
! The code simulates the following coupled ordinary nonlinear differential equations:
! \begin{equnarray}
! \frac{dS}{dt} = - R_0*S*I
! \frac{dI}{dt} = R_0*S*I - I
! \frac{dR}{dt} = I
! \end{equnarray}
! where, $R_0 = f(I) = N_0 + \frac{\epsa*I^{1+\alpha}}{1+\epsa * \delta * I^{1+\alpha}}$ and 'S' denotes 'Suspected', 'I' denotes 'Infected' and 'R' denotes 'Recovered' fraction of a population.
! Reference: 
! a. Equation 3 of http://dx.doi.org/10.1155/2014/12062
! b. Equation 4.7 of Mathematical Biology, Lecture notes for MATH 4333, Jeffrey R. Chasnov, Hong Kong University of Science and Technology.
! To run the code: gfortran <program>.f95; ./a.out
! To Plot: gnuplot
! p "Time_Evolution.dat" u 1:2 w lp, "Time_Evolution.dat" u 1:3 w lp, "Time_Evolution.dat" u 1:4 w lp, "Time_Evolution.dat" u 1:($2+$3+$4) w lp 
!
! 
implicit none
integer, parameter :: n = 3
integer i
real*8 x,y(n),yp(n),xst,xmax,h

open(unit=10,file='Time_Evolution.dat',status='unknown')

y(1) = 0.80d0 ! S = Suspected
y(2) = 0.05d0 ! I = Infected
y(3) = 0.00d0 ! R = Recovered

xst = 0.0d0
xmax = 50.0d0
h = 0.0010d0

do x=xst,xmax,h
  
  call derivative(x,y,yp,n)
  call rk4(x,h,y,yp,n)
    
  write(10,*) x+h,(y(i),i=1,n)
  
  call flush (10)

enddo

end program SIR

!================ DIFFERENTIAL EQUATION =============================================

subroutine derivative(x,y,yp,n)
REAL*8 y(n),yp(n),x,R0,N0,alpha,epsa,delta
integer n,i

  alpha = 0.50d0
  epsa = 0.030d0
  delta = 0.000010d0
  N0 = 2.0d0
  R0 = N0 + ( epsa*(y(2))**(1.0d0+alpha) ) / ( 1.0d0+epsa*delta*((y(2))**(1.0d0+alpha)) )

  yp(1) = -R0*y(1)*y(2)         ! dS/dt
  yp(2) = R0*y(1)*y(2) - y(2)   ! dI/dt
  yp(3) = y(2)                  ! dR/dt
  
  return
end 

!==================== RK 4 SUBROUTINE ===============================================

subroutine rk4(x,h,y,yp,n)
implicit none
real*8 x,h,y(n),yp(n),k1(n),k2(n),k3(n),k4(n),dum(n)
integer n,i

 do i = 1,n,1
  k1(i) = yp(i)
  dum(i) = y(i) + k1(i)*h/2.0d0
 enddo
   call derivative(x+h/2.0d0,dum,yp,n)
 do i = 1,n,1
  k2(i) = yp(i)
  dum(i) = y(i) + k2(i)*h/2.0d0
 enddo
   call derivative(x+h/2.0d0,dum,yp,n)
 do i = 1,n,1
  k3(i) = yp(i)
  dum(i) = y(i) + k3(i)*h
 enddo
   call derivative(x+h,dum,yp,n)
 do i = 1,n,1
  k4(i) = yp(i)
  y(i) = y(i) + h/6.0d0*(k1(i) + 2.0d0*k2(i) + 2.0d0*k3(i) + k4(i))
 enddo
 
 return
end

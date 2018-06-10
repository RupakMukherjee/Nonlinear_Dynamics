Program spiral
! Complex Ginzburg-Landau Equation (Pseudo-Spectral Solution)
! https://www.youtube.com/watch?v=DeGvax-jalc&t=16s
use omp_lib
implicit none

integer ( kind = 4 ), parameter :: Nx = 1536
integer ( kind = 4 ), parameter :: Ny = 1536
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,C1,C3,E
real ( kind = 8 ) x(Nx),y(Ny),Ax(Nx,Ny),Ay(Nx,Ny),Ax_new(Nx,Ny),Ay_new(Nx,Ny),R(Nx,Ny),Theta(Nx,Ny)
real ( kind = 8 ) NLx_new(Nx,Ny),NLx_old(Nx,Ny),NLy_new(Nx,Ny),NLy_old(Nx,Ny)
integer ( kind = 8 ) plan_forward,plan_backward ! FFTW
integer ( kind = 4 ) iret,thread_num,num_thread,proc_num ! OMP
real ( kind = 8 ) t1,t2
 character (len=90) :: filename

common/comm/iret,thread_num,Lx,Ly,C1,C3,dt

integer,parameter :: seed = 99999999
call srand(seed)

open(unit=5,file='System_information.dat',status='unknown')
open(unit=10,file='Initial_Condition.dat',status='unknown')
open(unit=20,file='Energy.dat',status='unknown')

!===================== USER INPUTS ============================================		

! Define Number of Threads.
proc_num = omp_get_num_procs()
thread_num = 26
call omp_set_num_threads (thread_num)

write(5,*) "Number of Processors=", proc_num, "Number of Threads=", thread_num  
  call flush (5)

Lx = 384.0*pi
Ly = 384.0*pi
dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)
  C1 = 2.0d0
  C3 = 0.850d0 !0.50d0

time_min = 0.00d0
time_max = 1000.0d0 
dt = 0.0010d0    

!===================== INITIAL TIME DATA ===============================================

do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0+real(j-1)*dy
    Ax(i,j) = dsin(2.0d0*pi*x(i)/Lx)
    Ay(i,j) = dcos(2.0d0*pi*y(j)/Ly)
    !Ax(i,j) = 2.0d0*rand()-1.0d0
    !Ay(i,j) = 2.0d0*rand()-1.0d0
    R(i,j) = dsqrt(Ax(i,j)*Ax(i,j)+Ay(i,j)*Ay(i,j))
    if (Ax(i,j) .ge. 0.0d0 .and. Ay(i,j) .ge. 0.0d0) then
	  Theta(i,j) = datan(abs(Ay(i,j)/Ax(i,j)))
	elseif (Ax(i,j) .lt. 0.0d0 .and. Ay(i,j) .ge. 0.0d0) then
	  Theta(i,j) = pi - datan(abs(Ay(i,j)/Ax(i,j)))
	elseif (Ax(i,j) .lt. 0.0d0 .and. Ay(i,j) .lt. 0.0d0) then
	  Theta(i,j) = pi + datan(abs(Ay(i,j)/Ax(i,j)))
	elseif (Ax(i,j) .ge. 0.0d0 .and. Ay(i,j) .lt. 0.0d0) then
	  Theta(i,j) = 2.0d0*pi - datan(abs(Ay(i,j)/Ax(i,j)))
	endif
    write(10,*) x(i),y(j),Ax(i,j),Ay(i,j),R(i,j),Theta(i,j)
    NLx_old(i,j) = 0.0d0
    NLy_old(i,j) = 0.0d0
  end do
end do

  close (10)

!======================= MAIN PROGRAM =====================================================

t1 = omp_get_wtime()

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

!if (mod(t,100) == 0) then
!print*, t/100
!endif

!=================================================================================
  call derive (Nx,Ny,Nh,time,Ax,Ay,Ax_new,Ay_new,NLx_new,NLx_old,NLy_new,NLy_old) 
  call ab (Nx,Ny,Nh,time,Ax,Ay,Ax_new,Ay_new,NLx_new,NLx_old,NLy_new,NLy_old)
!=================================================================================

!$OMP PARALLEL SHARED(Ax,Ay,Ax_new,Ay_new,NLx_new,NLy_new,NLx_old,NLy_old) PRIVATE(i,j)
!$OMP DO 
do i = 1,Nx
  do j = 1,Ny
    Ax(i,j) = Ax_new(i,j)
    Ay(i,j) = Ay_new(i,j)
    NLx_old(i,j) = NLx_new(i,j)
    NLy_old(i,j) = NLy_new(i,j)
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL 

E = 0.0d0

if (mod(dfloat(t),0.1d0/dt) == 0.0) then
do i = 1,Nx
  do j = 1,Ny
    E = E + (Ax(i,j)**2+Ay(i,j)**2)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo  
endif

if (mod(dfloat(t),0.10d0/dt) == 0.0) then
do i = 1,Nx
  do j = 1,Ny
    R(i,j) = dsqrt(Ax(i,j)*Ax(i,j)+Ay(i,j)*Ay(i,j))
      if (Ax(i,j) .ge. 0.0d0 .and. Ay(i,j) .ge. 0.0d0) then
	    Theta(i,j) = datan(abs(Ay(i,j)/Ax(i,j)))
	  elseif (Ax(i,j) .lt. 0.0d0 .and. Ay(i,j) .ge. 0.0d0) then
	    Theta(i,j) = pi - datan(abs(Ay(i,j)/Ax(i,j)))
	  elseif (Ax(i,j) .lt. 0.0d0 .and. Ay(i,j) .lt. 0.0d0) then
	    Theta(i,j) = pi + datan(abs(Ay(i,j)/Ax(i,j)))
	  elseif (Ax(i,j) .ge. 0.0d0 .and. Ay(i,j) .lt. 0.0d0) then
	    Theta(i,j) = 2.0d0*pi - datan(abs(Ay(i,j)/Ax(i,j)))
      endif
    write (filename, '( "/scratch/Rupak/Chimera/Sin/C3_0p85/fort.", I8.8 )' ) t+100
    open(unit=t+100,file=filename,status='unknown')
    write(t+100,*) x(i),y(j),Ax(i,j),Ay(i,j),R(i,j),Theta(i,j)
  enddo
enddo 
endif
 
if (mod(dfloat(t),0.1d0/dt) == 0.0) then
  write(20,*) time, E
  call flush(20)
endif

if (mod(dfloat(t),0.10d0/dt) == 0.0) then
  close (t+100)
endif

enddo ! time

t2 = omp_get_wtime()

write(5,*) "Time taken for the run =",(t2 - t1)/(60.0d0*60.0d0),"Hours"
  close(5)

contains

!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive(Nx,Ny,Nh,time,Ax,Ay,Ax_new,Ay_new,NLx_new,NLx_old,NLy_new,NLy_old)
implicit none
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0
integer ( kind = 4 ) Nx,Ny,Nh,i,j,iret,thread_num
real ( kind = 8 ) time,dt,Lx,Ly,C1,C3,kx,ky
real ( kind = 8 ) Ax(Nx,Ny),Ax_new(Nx,Ny),Ay(Nx,Ny),Ay_new(Nx,Ny)
real ( kind = 8 ) NLx_new(Nx,Ny), NLx_old(Nx,Ny)
real ( kind = 8 ) NLy_new(Nx,Ny), NLy_old(Nx,Ny)
real ( kind = 8 ) Ax_dum(Nx,Ny),Ay_dum(Nx,Ny),kx2_Ax(Nx,Ny),kx2_Ay(Nx,Ny),ky2_Ax(Nx,Ny),ky2_Ay(Nx,Ny),A2(Nx,Ny)
complex ( kind = 8 ) Axk(Nh,Ny),Ayk(Nh,Ny),kx2_Axk(Nh,Ny),kx2_Ayk(Nh,Ny),ky2_Axk(Nh,Ny),ky2_Ayk(Nh,Ny)
integer ( kind = 8 ) plan_forward,plan_backward

common/comm/iret,thread_num,Lx,Ly,C1,C3,dt

!$OMP PARALLEL SHARED(A2,Ax,Ay,Ax_dum,Ay_dum) PRIVATE(i,j)
!$OMP DO 
do i = 1,Nx
  do j = 1,Ny
    A2(i,j) = Ax(i,j)*Ax(i,j) + Ay(i,j)*Ay(i,j)
    Ax_dum(i,j) = Ax(i,j)
    Ay_dum(i,j) = Ay(i,j)
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL 

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Ax_dum, Axk, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Ay_dum, Ayk, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)
  
!$OMP PARALLEL SHARED(Lx,Ly,Axk,Ayk,kx2_Axk,kx2_Ayk,ky2_Axk,ky2_Ayk) PRIVATE(i,j,kx,ky)
!$OMP DO 
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    kx2_Axk(i,j) = - kx*kx*Axk(i,j)
    kx2_Ayk(i,j) = - kx*kx*Ayk(i,j)
    ky2_Axk(i,j) = - ky*ky*Axk(i,j)
    ky2_Ayk(i,j) = - ky*ky*Ayk(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    kx2_Axk(i,j) = - kx*kx*Axk(i,j)
    kx2_Ayk(i,j) = - kx*kx*Ayk(i,j)
    ky2_Axk(i,j) = - ky*ky*Axk(i,j)
    ky2_Ayk(i,j) = - ky*ky*Ayk(i,j)
  enddo
enddo 
!$OMP END DO
!$OMP END PARALLEL 

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, kx2_Axk, kx2_Ax, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, kx2_Ayk, kx2_Ay, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ky2_Axk, ky2_Ax, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ky2_Ayk, ky2_Ay, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
!$OMP PARALLEL SHARED(kx2_Ax,kx2_Ay,ky2_Ax,ky2_Ay) PRIVATE(i,j)
!$OMP DO 
do i = 1,Nx
  do j = 1,Ny
    kx2_Ax(i,j) = kx2_Ax(i,j)/(dfloat(Nx)*dfloat(Ny))
    kx2_Ay(i,j) = kx2_Ay(i,j)/(dfloat(Nx)*dfloat(Ny))
    ky2_Ax(i,j) = ky2_Ax(i,j)/(dfloat(Nx)*dfloat(Ny))
    ky2_Ay(i,j) = ky2_Ay(i,j)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL 

!$OMP PARALLEL SHARED(C1,C3,A2,Ax,Ay,kx2_Ax,kx2_Ay,ky2_Ax,ky2_Ay,NLx_new,NLy_new) PRIVATE(i,j)
!$OMP DO 
do i = 1,Nx
  do j = 1,Ny  
    NLx_new(i,j) = Ax(i,j) + kx2_Ax(i,j) + ky2_Ax(i,j) - C1*(kx2_Ay(i,j) + ky2_Ay(i,j)) - A2(i,j)*Ax(i,j) - C3*A2(i,j)*Ay(i,j)
    NLy_new(i,j) = Ay(i,j) + kx2_Ay(i,j) + ky2_Ay(i,j) + C1*(kx2_Ax(i,j) + ky2_Ax(i,j)) - A2(i,j)*Ay(i,j) + C3*A2(i,j)*Ax(i,j)
  enddo
enddo   
!$OMP END DO
!$OMP END PARALLEL 

return

end subroutine derive

!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab (Nx,Ny,Nh,time,Ax,Ay,Ax_new,Ay_new,NLx_new,NLx_old,NLy_new,NLy_old)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh,i,j,iret,thread_num
real ( kind = 8 ) time,dt,Lx,Ly,C1,C3
real ( kind = 8 ) Ax(Nx,Ny),Ax_new(Nx,Ny),Ay(Nx,Ny),Ay_new(Nx,Ny)
real ( kind = 8 ) NLx_new(Nx,Ny), NLx_old(Nx,Ny)
real ( kind = 8 ) NLy_new(Nx,Ny), NLy_old(Nx,Ny)

common/comm/iret,thread_num,Lx,Ly,C1,C3,dt

!$OMP PARALLEL SHARED(Ax,Ay,Ax_new,Ay_new,NLx_new,NLy_new,NLx_old,NLy_old) PRIVATE(i,j)
!$OMP DO 
do i = 1,Nx
  do j = 1,Ny
    Ax_new(i,j) = Ax(i,j) + ( (3.0d0/2.0d0)*NLx_new(i,j) - (1.0d0/2.0d0)*NLx_old(i,j) )*dt
    Ay_new(i,j) = Ay(i,j) + ( (3.0d0/2.0d0)*NLy_new(i,j) - (1.0d0/2.0d0)*NLy_old(i,j) )*dt
  end do
end do
!$OMP END DO
!$OMP END PARALLEL 

return

end subroutine ab

!====================================================================================

end program spiral  


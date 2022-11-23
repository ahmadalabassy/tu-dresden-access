! Program for solving a quadratic system of linear equations
!***************************************************************************
!***************************************************************************

program quad_system

use global_variables										! size of the system n, coefficients matrix M(n,n), right handside vector r

implicit none

real, allocatable:: MM(:,:),M_tria(:,:)						!coefficients matrix
real, allocatable:: rr(:),xx(:),r_tria(:)					!right side and the unknowns x vectors 
real start,finish,time,t_max,precond_time,tria_time,sol_time!variables for computing runtime
integer             i,j                             		!loop indices
integer				p,cholesky_limit						!size of systems
logical             det_exist,sym,pos_def_det    			!termination variables for cholecsky

! Setting the size of the size of the system
!***************************************************************************
p=60
cholesky_limit=10

! Opening the input file and reading the size of the system
!***************************************************************************
call input(p)
open(10, file='system.txt')
read(10,*)													! size of system already stored in p

! Allocation of the matrix and vector sizes
!***************************************************************************
allocate(MM(p,p))
allocate(rr(p))
allocate(xx(p))

! Reading the coefficients matrix and the vector of the right side
!***************************************************************************
do i=1,p
		read(10,*) (MM(i,j), j=1,p)
end do

do i=1,p
		read(10,*) rr(i)
end do

! Closing the input file
!***************************************************************************
close(10)

! Opening the input files of Gaussian and Cholesky runtime results
!**************************************************************************!
open(10,file='data_cholesky.txt') !will be created as empty fle if it doesn't exist
open(11,file='data_gaussian.txt')
open(12,file='preconditions_cholesky.txt')
open(13,file='triangular_matrix_cholesky.txt')
open(14,file='solution_cholesky.txt')
open(15,file='triangular_matrix_gaussian.txt')
open(16,file='solution_gaussian.txt')

! Creating systems of linear equations of size 1 to n then solving them
!***************************************************************************
do n=1,p

allocate(M(n,n))
allocate(r(n))
allocate(M_tria(n,n))
allocate(r_tria(n))
allocate(x(n))

write(*,*) 'n=', n

if (n.eq.p) then
	M=MM
	r=rr
else
	do i=1,n
		r(i)=rr(i)
		do j=1,n
			M(i,j)=MM(i,j)
		end do
	end do
end if

!limiting maximum size of system to 12 for Cholesky decomposition
!***************************************************************************
if (n.le.cholesky_limit) then
	
	! Calling the subroutine for the solution of the system of equations by Cholesky decomposition
	!***************************************************************************
	write(*,*) 'cholesky: '
	call cpu_time(start)
	call cholesky_precondition(sym,pos_def_det)
	call cpu_time(precond_time)
	if((sym.eqv..true.).and.(pos_def_det.eqv..true.)) then
		call cholesky_triangular_matrix(M_tria)
		call cpu_time(tria_time)
		call cholesky_solution(M_tria)
	end if
	call cpu_time(finish)

	! Output of Cholesky elimination
	!***************************************************************************
	write(*,*) 'symmetry criterion: ',sym
	write(*,*) 'positive definiteness for determinant criterion: ',pos_def_det

	if((sym.eqv..true.).and.(pos_def_det.eqv..true.)) then
			write(*,*) 'solution=',x
			time=finish-start
			sol_time=finish-tria_time
			tria_time=tria_time-precond_time
			precond_time=precond_time-start
			write(*,*) 'Total computing runtime =', time,'seconds'
			write(*,*) 'Preconditions check runtime =', precond_time,'seconds'
			write(*,*) 'Triangular matrix computing runtime =', tria_time,' seconds'
			write(*,*) 'Solution computing runtime =', sol_time,' seconds'
	else
			write(*,*) 'Preconditions for Cholesky decomposition are not met'
	end if

	! Writing data in cholesky input files
	!**************************************************************************!
	write(10,*) n,time
	write(12,*) n,precond_time
	write(13,*) n,tria_time
	write(14,*) n,sol_time
	if(n.eq.cholesky_limit) t_max=time
	
end if

! Calling the subroutine for the solution of the system of equations by Gaussian elimination
!***************************************************************************
write(*,*) 'gaussian: '
call cpu_time(start)
do i=1,10000
	call gauss_triangular_matrix(det_exist,M_tria,r_tria)
end do
call cpu_time(tria_time)
do i=1,10000
	call gauss_solution(M_tria,r_tria)
end do
call cpu_time(finish)

! Output of Gaussian elimination
!***************************************************************************
if(det_exist.eqv..true.) then
		write(*,*) 'solution=',x
else
		write(*,*) 'no unique solution'
end if

time=finish-start
sol_time=finish-tria_time
tria_time=tria_time-start

write(*,*) 'Total computing runtime =', time,'/10000 seconds'
write(*,*) 'Triangular matrix computing runtime =', tria_time,'/10000 seconds'
write(*,*) 'Solution computing runtime =', sol_time,'/10000 seconds'
write(*,*)

! Writing data in gaussian input files
!**************************************************************************!
write(11,*) n,time
write(15,*) n,tria_time
write(16,*) n,sol_time

deallocate(M)
deallocate(M_tria)
deallocate(r)
deallocate(r_tria)
deallocate(x)

end do

! Closing the input files of Gaussian and Cholesky runtime results
!**************************************************************************!
close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)

! Plotting data
!**************************************************************************!
call gnuplot_2D(p,real(int(t_max)+1))

end
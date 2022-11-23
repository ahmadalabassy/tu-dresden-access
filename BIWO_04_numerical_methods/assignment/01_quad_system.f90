! Program for solving a quadratic system of linear equations
!***************************************************************************
!***************************************************************************

program quad_system

use global_variables							! size of the system n, coefficients matrix M(n,n), right handside vector r

implicit none

real start,finish								!variables for computing runtime
integer             i,j                         ! loop indices
logical             det_exist,sym,pos_def_det   ! termination variable

! Opening the input file and reading the size of the system
!***************************************************************************
open(10, file='system.txt')
read(10,*) n

! Allocation of the matrix and vector sizes
!***************************************************************************
allocate(M(n,n))
allocate(r(n))
allocate(x(n))

! Reading the coefficients matrix and the vector of the right side
!***************************************************************************
do i=1,n
		read(10,*) (M(i,j), j=1,n)
end do

do i=1,n
		read(10,*) r(i)
end do

! Closing the input file
!***************************************************************************
close(10)

! Calling the subroutine for the solution of the system of equations by Cholesky decomposition
!***************************************************************************
call cpu_time(start)
!do i=1,10000
	call cholesky(sym,pos_def_det)
!end do
call cpu_time(finish)

! Output of Cholesky decomposition
!***************************************************************************
write(*,*) 'cholesky: '
write(*,*) 'symmetry criterion: ',sym
write(*,*) 'positive definiteness for determinant criterion: ',pos_def_det

if((sym.eqv..true.).and.(pos_def_det.eqv..true.)) then
		write(*,*) 'solution=',x
		write(*,*) 'computing runtime =', finish-start,'seconds'!,'/1000'
else
		write(*,*) 'Preconditions for Cholesky decomposition are not met'
end if

! Calling the subroutine for the solution of the system of equations by Gaussian elimination
!***************************************************************************
call cpu_time(start)
do i=1,10000
	call gauss(det_exist)
end do
call cpu_time(finish)

! Output of Gaussian elimination with a prompt as a "pause"
!***************************************************************************
write(*,*) 'gaussian: '
if(det_exist.eqv..true.) then
		write(*,*) 'solution=',x
else
		write(*,*) 'no unique solution'
end if
write(*,*) 'computing runtime =', finish-start,'/10000 seconds'
read(*,*)

end
! Subroutine for solving a quadratic system of linear equations M*x=r
! by using the cholesky decomposition
!**************************************************************************!
!**************************************************************************!

subroutine cholesky_solution(C)

use global_variables

implicit none

real C(n,n)							! the lower trianglular matrix in cholesky decomposition													
real   y(n)                 		! vector of the right side corresponding to matrix A 
real temp
integer i,k               	    	! loop indexes

! solving vector y
!***************************************************************************
do i=1,n
	temp=0.0
	do k=1,i-1
		temp = temp+C(i,k)*y(k)
	end do
	y(i)=(r(i)-temp)/C(i,i)
end do

! Solving the system of equations in the triangular form C_t*x=y
!**************************************************************************!
do i=n,1,-1
		temp=0.0
		do k=i+1,n
			temp=temp+C(k,i)*x(k)
		end do
		x(i)=(y(i)-temp)/C(i,i)
end do

end subroutine
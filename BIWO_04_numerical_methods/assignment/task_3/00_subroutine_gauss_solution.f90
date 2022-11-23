! Subroutine for solving a quadratic system of linear equations M*x=r
! by using the Gaussian algorithm with partial pivoting
!**************************************************************************!
!**************************************************************************!

subroutine gauss_solution(A,rr)

use global_variables

implicit none

real   A(n,n)                ! matrix M in triangular form
real   rr(n)                 ! vector of the right side corresponding to matrix A 
real   d	                 ! temporary variables
integer i,j                  ! loop indexes

! Solving the system of equations in the triangular form A*x=rr
!**************************************************************************!

x(n)=rr(n)/A(n,n)
do i=n-1,1,-1
		d=0.0
		do j=n,i+1,-1
			d=d+A(i,j)*x(j)
		end do
		x(i)=(rr(i)-d)/A(i,i)
end do

end subroutine







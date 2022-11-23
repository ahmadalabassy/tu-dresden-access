! Subroutine for solving a quadratic system of linear equations M*x=r
! by using the cholesky decomposition
!**************************************************************************!
!**************************************************************************!

subroutine cholesky_triangular_matrix(C)

use global_variables

implicit none

real C(n,n)										! the lower trianglular matrix in cholesky decomposition													
real temp
integer i,j,k               	    			! loop indexes

! Assigning initial values to variables
!***************************************************************************
C=0

! Creating the lower triangular matrix C
!***************************************************************************
do i=1,n
	do j=1,n
		temp=0.0
		if(i.eq.j) then							!elements of the main diagonal
			do k=1,i-1
				temp = temp+C(i,k)**2
			end do
			C(i,j)=sqrt(M(i,i)-temp)
		else if (i.gt.j) then					!lower traingle elements
			do k=1,j-1
				temp = temp+C(i,k)*C(j,k)
			end do
			C(i,j)=(M(i,j)-temp)/C(j,j)
		end if
	end do
end do

end subroutine
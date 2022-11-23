! Subroutine for solving a quadratic system of linear equations M*x=r
! by using the cholesky decomposition
!**************************************************************************!
!**************************************************************************!

subroutine cholesky_precondition(sym_check,pos_def_dcheck)

use global_variables

implicit none
											
real det, laplace_det
integer i,j	        							! loop indexes
logical sym_check,pos_def_dcheck				! termination variables

! Check of matrix symmetry criterion
!**************************************************************************!
sym_check=.true.
do i=1,n
	do j=i+1,n									!less computation time we just need to check one half
			if(M(i,j).ne.M(j,i)) then
					sym_check = .false.
					exit
			end if
	end do
end do

! Check of positive definiteness for the determinant criterion
!**************************************************************************!
pos_def_dcheck=.true.
do i=1,n
		det = laplace_det(M(1:i,1:i),i)			!assigning values for minor matrix and calculating determinant using laplace expansion
		if (det.le.0) then						!check if determinant is less than or equal zero
			pos_def_dcheck = .false.
			exit
		end if
end do

end subroutine
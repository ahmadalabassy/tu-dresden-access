! Subroutine for solving a quadratic system of linear equations M*x=r
! by using the Gaussian algorithm with partial pivoting
!**************************************************************************!
!**************************************************************************!

subroutine gauss(det_exist)

use global_variables

implicit none

real   A(n,n)                ! matrix M in triangular form
real   rr(n)                 ! vector of the right side corresponding to matrix A 
real   d, temp               ! temporary variables
integer i,j,k                ! loop indexes
integer piv                  ! pivot row 
logical det_exist            ! termination variable

! Default values of the variables
!**************************************************************************!
A=M
rr=r
det_exist=.true.


! Transformation of the matrix A in triangular form
!**************************************************************************!
do k=1,n-1
		temp=A(k,k)                                         ! "memorizing" the diagonal element and
        piv=k                                               ! taking row k as pivot row
        do i=k+1,n          
				if(Abs(A(i,k)).gt.Abs(temp))   then         ! determining the largest element of
						temp=A(i,k)							! remaining column
                        piv=i                               ! corresponging row = pivot row
				end if
		end do

        if(piv.ne.k) then                                   ! if pivot row |= row k, then
				do j=k,n								    ! swap pivot row and row k
						temp=A(piv,j)						! 
                        A(piv,j)=A(k,j)                     ! - of matrix A
                        A(k,j)=temp                         ! 
                end do                                      !
                temp=rr(piv)                                ! - of the vector of the right side rr
                rr(piv)=rr(k)                               !
                rr(k)=temp                                  !
		end if

        if(A(k,k).eq.0) then                                ! if all elements of the remaining column = 0,
				det_exist=.false.							! then the det(A) is = 0
                return                                      ! return to the calling (sub-)program
		end if

        do j=k+1,n                                          ! subtracting the k-th row (with
				d=A(j,k)/A(k,k)								! factor d) from the remaining rows
                do i=k,n                                    ! - of matrix A (= diagonalization)
						A(j,i)=A(j,i)-d*A(k,i)				!
                end do                                      ! - of the vector of the right side rr
				rr(j)=rr(j)-d*rr(k)
		end do
end do

if(A(n,n).eq.0) then                                        ! if diagonal element(n,n) = 0,
        det_exist=.false.                                   ! then det(A) is = 0
        return                                              ! return to the calling (sub-)program
end if

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







! Function for calculating the determinant 
!**************************************************************************!
!**************************************************************************!

function elimination_det(A,n) RESULT(det)

implicit none

real    A(n,n)      ! matrix for which the determinant is calculated !call by reference.. any changes here reflect on the matrix A in the main programme
real	A_tria(n,n) ! matrix A in triangular form
real    det         ! determinant of matrix A i.e. A_tria
real    temp		! temporary variables
integer n           ! size of matrix A
integer i,j,k       ! loop indeces
integer piv			! pivot row
logical det_exist   ! termination variable

! Default values for the variables
!**************************************************************************!
det_exist=.true.
det=1.0
A_tria=A

! Conversion of matrix A in triangular form (without pivoting)
!**************************************************************************!
do k=1,n-1
		temp=A_tria(k,k)                                  					! "memorizing" the diagonal element and
        piv=k                                               				! taking row k as pivot row
        do i=k+1,n          
				if(Abs(A_tria(i,k)).gt.Abs(temp))   then    				! determining the largest element of
						temp=A_tria(i,k)									! remaining column
                        piv=i                               				! corresponging row = pivot row
				end if
		end do

        if(piv.ne.k) then                                   				! if pivot row |= row k, then
				do j=k,n								    				! swap pivot row and row k
						temp=A_tria(piv,j)									! 
                        A_tria(piv,j)=A_tria(k,j)           				!
                        A_tria(k,j)=temp                    				! 
                end do 
				det_exist=.true.											! change of sign of the determinant
				det=-det													! due to row swapping
		end if

		if(A_tria(k,k).eq.0) then											! all elements of the column k are null,
                det=0.0														! it means that det(A)=0
                return														! return to the calling (sub-)program
		end if

        do i=k+1,n															! subtracting the k-th row (with
                temp=A_tria(i,k)/A_tria(k,k)								! factor temp from the remaining rows
                do j=k,n													! (= diagonalization)
						A_tria(i,j)=A_tria(i,j)-temp*A_tria(k,j)
				end do
		end do
end do

! Calculation of the determinant as a product of the diagonal elements
!**************************************************************************!
do i=1,n
		det=det*A_tria(i,i)
end do
end function

        


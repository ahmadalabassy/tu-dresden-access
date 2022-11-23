! Subroutine for solving a quadratic system of linear equations M*x=r
! by using the cholesky decomposition
!**************************************************************************!
!**************************************************************************!

subroutine cholesky(sym_check,pos_def_dcheck)

use global_variables

implicit none

real, allocatable:: C(:,:)							! the lower trianglular matrix in cholesky decomposition													
real   y(n)                 						! vector of the right side corresponding to matrix A 
real det, laplace_det, temp
integer i,j,k               	    				! loop indexes
logical sym_check,pos_def_dcheck   ! termination variables

! Check of matrix symmetry criterion
!**************************************************************************!
sym_check=.true.
do i=1,n
	do j=i+1,n										!less computation time we just need to check one half
			if(M(i,j).ne.M(j,i)) then
					sym_check = .false.
					exit
			end if
	end do
end do

! Check of positive definiteness for the determinant criterion
!**************************************************************************!
pos_def_dcheck=.true.
do k=1,n
		det = laplace_det(M(1:k,1:k),k)			!assigning values for minor matrix and calculating determinant using laplace expansion
		if (det.le.0) then						!check if determinant is less than or equal zero
			pos_def_dcheck = .false.
			exit
		end if
end do

! Check if conditions are met
!***************************************************************************
if((sym_check.eqv..false.).or.(pos_def_dcheck.eqv..false.)) return

! Assigning initial values to variables
!***************************************************************************
allocate(C(n,n))
C=0
y=0

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
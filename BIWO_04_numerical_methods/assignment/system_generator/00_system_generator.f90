!Program to generate positive definite matrices
!**************************************************************************!
!**************************************************************************!
program system

implicit none

real, allocatable:: A(:,:),Tria_lower(:,:),Tria_upper(:,:)	!A:output matrix of matrix multiplication C(Lower Triangula)
integer n													!size of the system
integer i,j													!loope indices

! Opening the input file and reading the size of the system
!************************************************************************
open(10,file='system.txt')

! Assigning initial values to variables
!***************************************************************************
n=10

! Generating the system file
!***************************************************************************

	! Allocation of the matrix and vector sizes
	!***************************************************************************
	allocate (Tria_lower(n,n))
	allocate (Tria_upper(n,n))
	allocate (A(n,n))

	! Generating the lower triangular matrix
	!***************************************************************************
	Tria_lower=0
	do i=1,n
		do j=1,n
			if(i.ge.j) then
				call random_number(Tria_lower(i,j))
				Tria_lower(i,j)=int(10*Tria_lower(i,j))
				if(i.eq.j.and.Tria_lower(i,j).eq.0.0) Tria_lower(i,j)=i
			end if
		end do
	end do

	! Generating the upper triangular matrix
	!***************************************************************************
	do i=1,n
		do j=1,n
			Tria_upper(i,j)=Tria_lower(j,i)
		end do
	end do

	! Generating the output positive definite matrix
	!***************************************************************************
	A=matmul(Tria_lower,Tria_upper)

	! Generating the output positive definite matrix
	!***************************************************************************
	write (10,*) n
	do i=1,n
		write(10,*) A(i,1:n)
	end do
	do i=1,n
		write(10,*) 1
	end do
	deallocate (A)
	deallocate (Tria_lower)
	deallocate (Tria_upper)

end
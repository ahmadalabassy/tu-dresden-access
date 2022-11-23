! Function to calculate determinant of a matrix using Laplace expansion
!**************************************************************************!
!**************************************************************************!

recursive real function laplace_det(A,n) RESULT(det)

implicit none

integer n													!size of the system
real A(n,n), A_sub(n-1,n-1)									!coefficients matrix and the sub cofactors matrix
integer i,j,k 												!loop indices

if(n.eq.1) then												!base case when n is 1
	det = A(n,n)
else if (n.eq.2) then										!base case when n is 2
	det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
else
	det= 0.0												!initial value
	do k=1,n												!calculating determinants of elements of the first row
        do i=1,n-1 											!i = row index
			do j=1,n-1  									!j = column index
				if(j.lt.k) then
					A_sub(i,j)=A(i+1,j)
				else if(j.ge.k) then
					A_sub(i,j)=A(i+1,j+1)
				end if
			end do
        end do
		det = det+((-1)**(1+k))*A(1,k)*laplace_det(A_sub,n-1) !recursive calculation of determinant
	end do 
end if

end function
! Module for declaration of global variables
!**************************************************************************!
!**************************************************************************!
module global_variables

implicit none

integer n								 ! size of the system
real, allocatable:: M(:,:)				 ! coefficients square matrix
real, allocatable:: r(:),x(:)            ! right side and the unknowns x vectors

end module global_variables
! Module for declaration of global variables
!**************************************************************************!
!**************************************************************************!
module global_variables

implicit none

integer n								 ! size of the system
real, allocatable:: M(:,:)				 ! coefficients square matrix
real, allocatable:: r(:),x(:)            ! vector of the right side

end module global_variables
!
! This is a simple test module to test
! running F90 code within the interface
!

subroutine EcoSIM_DataTest() bind(C)

  use, intrinsic :: iso_c_binding

  implicit none

  write(*,*) "Okay calling the function works."

end subroutine EcoSIM_DataTest

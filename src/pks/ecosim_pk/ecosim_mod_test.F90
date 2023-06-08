!
! This is a simple test module to test
! running F90 code within the interface
!

module ecosim_datatest_mod
  implicit none
contains
  subroutine ecosim_datatest() bind(C)

    use, intrinsic :: iso_c_binding

    write(*,*) "Okay calling the function works."

  end subroutine EcoSIM_DataTest

end module ecosim_datatest_mod

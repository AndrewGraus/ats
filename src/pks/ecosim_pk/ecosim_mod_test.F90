!
! This is a simple test module to test
! running F90 code within the interface
!

module ecosim_datatest_mod
  implicit none
contains
  subroutine ecosim_datatest(col) bind(C)

    use, intrinsic :: iso_c_binding
    integer (c_int), value, intent(in) :: col

    write(*,*) "Okay calling the ncol function works."
    write(*,*) "num col is: ", col

  end subroutine ecosim_datatest

end module ecosim_datatest_mod

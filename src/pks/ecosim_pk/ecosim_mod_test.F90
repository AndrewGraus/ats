!
! This is a simple test module to test
! running F90 code within the interface
!

module ecosim_datatest_mod
  use BGCContainers_module, only : BGCProperties, BGCSizes

  implicit none

contains
  subroutine ecosim_datatest(col, props) bind(C)

    use, intrinsic :: iso_c_binding
    integer (c_int), value, intent(in) :: col
    type(BGCProperties) :: props
    integer :: i
    integer :: len

    write(*,*) "Okay calling the props function works."
    write(*,*) "num col is: ", col
    write(*,*) "printing the BGCProps data: "

    len = props%volume%size

    write(*,*) "the length is: ", len
    !do i = 1, props%volume%size
    !  write(*,*) props%volume%data(i)
    !end do

    write(*,*) "the properties are finished"

  end subroutine ecosim_datatest

end module ecosim_datatest_mod

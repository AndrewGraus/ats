!
! This is a simple test module to test
! running F90 code within the interface
!

module ecosim_datatest_mod
  use BGCContainers_module, only : BGCSizes, &
   BGCProperties
  implicit none

  public :: &
      ecosim_datatest

  private :: &
      SetBGCSizes

contains
  subroutine ecosim_datatest(col, props, sizes) bind(C)

    use, intrinsic :: iso_c_binding
    use BGCContainers_module, only : BGCSizes, &
     BGCProperties
    integer (c_int), value, intent(in) :: col
    type(BGCProperties), intent(in) :: props
    type(BGCSizes), intent(out) :: sizes
    integer :: i
    integer :: len

    write(*,*) "calling set sizes"
    write(*,*) "does sizes exist?"
    write(*,*) sizes%num_components

    call SetBGCSizes(sizes)

    write(*,*) "Okay calling the props function works."
    write(*,*) "num col is: ", col
    write(*,*) "printing the BGCProps data: "

    len = props%volume%size

    write(*,*) "the length is: ", len
    do i = 1, len
      write(*,*) props%volume%data(i)
    end do

    write(*,*) "the properties are finished"

  end subroutine ecosim_datatest

  subroutine SetBGCSizes(sizes)

    use BGCContainers_module, only : BGCSizes

    implicit none

    type (BGCSizes), intent(out) :: sizes

    sizes%num_components = 1
    sizes%ncells_per_col_ = 100

  end subroutine SetBGCSizes

end module ecosim_datatest_mod

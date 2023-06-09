!
! This is a simple test module to test
! running F90 code within the interface
!

module ecosim_datatest_mod
  use BGCContainers_module, only : BGCProperties, BGCSizes

  implicit none

interface
  subroutine AllocateBGCProperties(sizes, props) bind(C, name='AllocateBGCProperties')
    use BGCContainers_module, only : BGCSizes, BGCProperties
    implicit none
    type(BGCSizes) :: sizes
    type(BGCProperties) :: props
  end subroutine
end interface
interface
  subroutine FreeBGCProperties(props) bind(C, name='FreeBGCProperties')
    use BGCContainers_module, only : BGCProperties
    implicit none
    type(BGCProperties) :: props
  end subroutine
end interface

contains
  subroutine ecosim_datatest(col, props) bind(C)

    use, intrinsic :: iso_c_binding
    integer (c_int), value, intent(in) :: col
    type(BGCProperties) :: props

    write(*,*) "Okay calling the props function works."
    write(*,*) "num col is: ", col
    write(*,*) "printing the BGCProps data: "
    do i = 1, props%size
      write(*,*) props%volume%data(i)
    end do

    write(*,*)

  end subroutine ecosim_datatest

end module ecosim_datatest_mod

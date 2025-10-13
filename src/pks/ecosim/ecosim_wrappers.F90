!
!There needs to be a wrapper for the eocsim f90 driver
!as there are differences between how gfortran and intel compilers
!handle mangling conventions.
!
! Copied from the alquimia wrapper:
! **************************************************************************** !
!
! PFloTran Alquimia Inteface Wrappers
!
! Author: Benjamin Andre
!
! Different fortran compilers use different name mangling conventions
! for fortran modules:
!
!    gfortran : ___modulename_MOD_procedurename
!
!    intel : _modulename_mp_procedurename_
!
!    as a consequence we can't put the alquimia interface into a
!    module and call it directly from C/C++. Instead we use
!    some simple wrapper functions.
!
! Notes:
!
!  * Function call signatures are dictated by the alquimia API!
!
!  * alquimia data structures defined in AlquimiaContainers_module
!    (alquimia_containers.F90) are dictated by the alquimia API.
!
! **************************************************************************** !

subroutine EcoSIM_DataTest() bind(c)
    use, intrinsic :: iso_c_binding
    !Simple test for wrapper functionality
    implicit none

    write(*,*) "in the data test"

  end subroutine EcoSIM_DataTest

! **************************************************************************** !

subroutine EcoSIM_Setup(properties, state, sizes, num_iterations,&
                        num_columns, ncells_per_col_) bind(C)

  use, intrinsic :: iso_c_binding

  use BGCContainers_module
  use ATSCPLMod, only : ATS2EcoSIMData, Init_EcoSIM, EcoSIM2ATSData

  implicit none

  ! function parameters
  !character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (BGCSizes), intent(out) :: sizes
  type (BGCState), intent(inout) :: state
  !type (BGCAuxiliaryData), intent(inout) :: aux_data
  type (BGCProperties), intent(in) :: properties
  integer, intent(inout) :: num_columns
  integer, intent(inout) :: num_iterations
  integer, intent(inout) :: ncells_per_col_

  call ATS2EcoSIMData(num_columns, state, properties, sizes)

  call Init_EcoSIM(sizes)

  call EcoSIM2ATSData(num_columns, state, sizes)

end subroutine EcoSIM_Setup

! **************************************************************************** !

subroutine EcoSIM_Shutdown() bind(C)

  !For now this does nothing, but it should clear all
  !the data structures
  use, intrinsic :: iso_c_binding

  use BGCContainers_module

  implicit none

  ! function parameters
  !character(kind=c_char), dimension(*), intent(in) :: input_filename
  !type (BGCSizes), intent(out) :: sizes
  !type (BGCState), intent(in) :: state
  !type (BGCAuxiliaryData), intent(in) :: aux_data
  !type (BGCProperties), intent(in) :: properties
  !integer :: num_columns, jz, js
  !integer, intent(in) :: num_iterations

end subroutine EcoSIM_Shutdown

! **************************************************************************** !

subroutine EcoSIM_Advance( &
     delta_t, &
     properties, &
     state, &
     sizes, &
     num_iterations, &
     num_columns) bind(C)

  use, intrinsic :: iso_c_binding
  use BGCContainers_module
  use ATSCPLMod, only : Run_EcoSIM_one_step, ATS2EcoSIMData, EcoSIM2ATSData

  implicit none

  ! function parameters
  real (c_double), value, intent(in) :: delta_t
  type (BGCProperties), intent(in) :: properties
  type (BGCState), intent(inout) :: state
  type (BGCSizes), intent(out) :: sizes
  !type (BGCEngineStatus), intent(out) :: status
  integer, intent(in) :: num_columns
  integer, intent(in) :: num_iterations

  call ATS2EcoSIMData(num_columns, state, properties, sizes)

  call Run_EcoSIM_one_step(sizes)

  call EcoSIM2ATSData(num_columns, state, sizes)

end subroutine EcoSIM_Advance

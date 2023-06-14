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

subroutine EcoSIM_Setup(props, state, sizes, &
                          num_iterations, ncol) bind(C)

  use, intrinsic :: iso_c_binding

  use BGCContainers_module
  use ATSCPLMod, only : ATS2EcoSIMData, Init_EcoSIM, EcoSIM2ATSData

  implicit none

  ! function parameters
  !character(kind=c_char), dimension(*), intent(in) :: input_filename
  type (BGCSizes), intent(out) :: sizes
  type (BGCState), intent(inout) :: state
  !type (BGCAuxiliaryData), intent(inout) :: aux_data
  type (BGCProperties), intent(in) :: props
  integer :: ncol, jz, js
  integer, intent(in) :: num_iterations

  write(*,*) "starting driver transfer ATS2EcoSIMData"

  call ATS2EcoSIMData(ncol, state, props, sizes)

  write(*,*) "starting driver Init_EcoSIM"

  !call Init_EcoSIM(jz,js,ncol)

  write(*,*) "starting driver transfer EcoSIM2ATSData"

  !call EcoSIM2ATSData()

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
  !type (BGCProperties), intent(in) :: props
  !integer :: ncol, jz, js
  !integer, intent(in) :: num_iterations

end subroutine EcoSIM_Shutdown

! **************************************************************************** !

subroutine EcoSIM_Advance( &
     delta_t, &
     props, &
     state, &
     aux_data, &
     num_iterations, &
     ncol) bind(C)

  use, intrinsic :: iso_c_binding
  use BGCContainers_module
  use ATSCPLMod, only : Run_EcoSIM_one_step, ATS2EcoSIMData, EcoSIM2ATSData

  implicit none

  ! function parameters
  real (c_double), value, intent(in) :: delta_t
  type (BGCProperties), intent(in) :: props
  type (BGCState), intent(inout) :: state
  type (BGCAuxiliaryData), intent(inout) :: aux_data
  !type (BGCEngineStatus), intent(out) :: status
  integer :: ncol
  integer, intent(in) :: num_iterations

  call ATS2EcoSIMData(ncol, state, aux_data, props)

  call Run_EcoSIM_one_step()

  call EcoSIM2ATSData()

end subroutine EcoSIM_Advance

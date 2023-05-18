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

subroutine EcoSIM_Setup(input_filename, hands_off, &
     pft_engine_state, sizes, functionality, status) bind(C)

  use, intrinsic :: iso_c_binding

  use BGCContainers_module
  use ATSCPLMod, only : ATS2EcoSIMData, Init_EcoSIM

  implicit none

  ! function parameters
  character(kind=c_char), dimension(*), intent(in) :: input_filename
  logical (c_bool), value, intent(in) :: hands_off
  type (c_ptr), intent(out) :: pft_engine_state
  type (BGCSizes), intent(out) :: sizes
  type (BGCEngineFunctionality), intent(out) :: functionality
  type (BGCEngineStatus), intent(out) :: status

  call ATS2EcoSIMData(filter_col,data_1d,var_1d,data_2d,var_2d,data_3d,var_3d)

  call Init_EcoSIM(jz,js,ncol)

  call EcoSIM2ATSData()

end subroutine EcoSIM_Setup

! **************************************************************************** !

subroutine EcoSIM_Advance( &
     pft_engine_state, &
     delta_t, &
     properties, &
     state, &
     aux_data, &
     status) bind(C)

  use, intrinsic :: iso_c_binding

  use BGCContainers_module
  use ATSCPLMod, only : Run_EcoSIM_one_step

  implicit none

  ! function parameters
  type (c_ptr), intent(inout) :: pft_engine_state
  real (c_double), value, intent(in) :: delta_t
  type (BGCProperties), intent(in) :: properties
  type (BGCState), intent(inout) :: state
  type (BGCAuxiliaryData), intent(inout) :: aux_data
  type (BGCEngineStatus), intent(out) :: status

  call ATS2EcoSIMData(filter_col,data_1d,var_1d,data_2d,var_2d,data_3d,var_3d)

  call Run_EcoSIM_one_step()

  call EcoSIM2ATSData()
  
end subroutine EcoSIM_Advance

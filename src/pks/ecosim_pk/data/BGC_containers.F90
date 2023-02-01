
!
! Alquimia Copyright (c) 2013-2016, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any
! required approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Alquimia is available under a BSD license. See LICENSE.txt for more
! information.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Technology Transfer and Intellectual Property
! Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
!
! NOTICE.  This software was developed under funding from the U.S. Department
! of Energy.  As such, the U.S. Government has been granted for itself and
! others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
! license in the Software to reproduce, prepare derivative works, and perform
! publicly and display publicly.  Beginning five (5) years after the date
! permission to assert copyright is obtained from the U.S. Department of Energy,
! and subject to any subsequent five (5) year renewals, the U.S. Government is
! granted for itself and others acting on its behalf a paid-up, nonexclusive,
! irrevocable, worldwide license in the Software to reproduce, prepare derivative
! works, distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Authors: Benjamin Andre <bandre@lbl.gov>
!

! **************************************************************************** !
!
! Alquimia Containers module
!
! Author: Benjamin Andre
!
! WARNINGS:
!
!   * The alquimia data structures defined in the this are dictated by
!     the alquimia API! Do NOT change them unless you make
!     corresponding changes to the c containers (and doc).
!
!   * The order of the data members matters! If num_primary is the
!     first member and num_minerals is the third, then they must be in
!     those positions on both the c and fortran side of the interface!
!
!   * The names of the containers and their members should be the same
!     for both c and fortran. The language interface doesn't require
!     this, but it makes it easier for the reader to understand what is
!     going on.
!
! **************************************************************************** !

module AlquimiaContainers_module

  use, intrinsic :: iso_c_binding

  implicit none

  integer (c_int), parameter :: kBGCMaxStringLength = 512
  integer (c_int), parameter :: kBGCMaxWordLength = 32

  integer (c_int), parameter :: kBGCNoError = 0
  integer (c_int), parameter :: kBGCErrorInvalidEngine = 1
  integer (c_int), parameter :: kBGCErrorUnknownConstraintName = 2
  integer (c_int), parameter :: kBGCErrorUnsupportedFunctionality = 3
  integer (c_int), parameter :: kBGCErrorEngineIntegrity = 4577

  character (13), parameter :: kBGCStringTotalAqueous = 'total_aqueous'
  character (12), parameter :: kBGCStringTotalSorbed = 'total_sorbed'
  character (25), parameter :: kBGCStringTotalAqueousPlusSorbed = 'total_aqueous_plus_sorbed'
  character (4), parameter :: kBGCStringFree = 'free'
  character (2), parameter :: kBGCStringPH = 'pH'
  character (7), parameter :: kBGCStringMineral = 'mineral'
  character (3), parameter :: kBGCStringGas = 'gas'
  character (6), parameter :: kBGCStringCharge = 'charge'

  type, public, bind(c) :: BGCVectorDouble
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type BGCVectorDouble

  type, public, bind(c) :: BGCVectorInt
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type BGCVectorInt

  type, public, bind(c) :: BGCVectorString
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type BGCVectorString

  type, public, bind(c) :: BGCSizes
     integer (c_int) :: num_primary
     integer (c_int) :: num_sorbed
     integer (c_int) :: num_minerals
     integer (c_int) :: num_aqueous_complexes
     integer (c_int) :: num_aqueous_kinetics
     integer (c_int) :: num_surface_sites
     integer (c_int) :: num_ion_exchange_sites
     integer (c_int) :: num_isotherm_species
     integer (c_int) :: num_aux_integers
     integer (c_int) :: num_aux_doubles
  end type BGCSizes

  type, public, bind(c) :: BGCState
     ! I think I have to write the data as vector doubles
     type (BGCVectorDouble) :: fluid_density
     type (BGCVectorDouble) :: gas_density
     type (BGCVectorDouble) :: ice_density
     type (BGCVectorDouble) :: porosity
     type (BGCVectorDouble) :: water_content
     type (BGCVectorDouble) :: temperature
     type (BGCVectorDouble) :: total_mobile
  end type BGCState

  type, public, bind(c) :: BGCProperties
     type (BGCVectorDouble) :: liquid_saturation
     type (BGCVectorDouble) :: gas_saturation
     type (BGCVectorDouble) :: ice_saturation
     type (BGCVectorDouble) :: elevation
     type (BGCVectorDouble) :: relative_permeability
     type (BGCVectorDouble) :: conductivity
     type (BGCVectorDouble) :: volume
  end type BGCProperties

  type, public, bind(c) :: BGCAuxiliaryData
     type (BGCVectorInt) :: aux_ints
     type (BGCVectorDouble) :: aux_doubles
  end type BGCAuxiliaryData

end module AlquimiaContainers_module

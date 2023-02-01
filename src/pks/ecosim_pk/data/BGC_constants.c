/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California,
** through Lawrence Berkeley National Laboratory (subject to receipt of any
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
**
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software,
** please contact Berkeley Lab's Technology Transfer and Intellectual Property
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
**
** NOTICE.  This software was developed under funding from the U.S. Department
** of Energy.  As such, the U.S. Government has been granted for itself and
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
** license in the Software to reproduce, prepare derivative works, and perform
** publicly and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of Energy,
** and subject to any subsequent five (5) year renewals, the U.S. Government is
** granted for itself and others acting on its behalf a paid-up, nonexclusive,
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly,
** and to permit others to do so.
**
** Authors: Benjamin Andre <bandre@lbl.gov>
*/


#include "BGC_constants.h"

/* String lengths */
const int kBGCMaxStringLength = 512;
const int kBGCMaxWordLength = 32;

/* Geochemistry Engine Strings */
const char* kBGCStringPFloTran = "PFloTran";
const char* kBGCStringCrunchFlow = "CrunchFlow";
const char* kBGCStringTotal = "total_aqueous";
const char* kBGCStringTotalSorbed = "total_sorbed";
const char* kBGCStringTotalAqueousPlusSorbed = "total_aqueous_plus_sorbed";
const char* kBGCStringFree = "free";
const char* kBGCStringPH = "pH";
const char* kBGCStringMineral = "mineral";
const char* kBGCStringGas = "gas";
const char* kBGCStringCharge = "charge";

/* Error Codes */
const int kBGCNoError = 0;
const int kBGCErrorInvalidEngine = 1;
const int kBGCErrorUnknownConstraintName = 2;
const int kBGCErrorUnsupportedFunctionality = 3;
const int kBGCErrorEngineIntegrity = 4577;

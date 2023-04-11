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

/*This code is a modification of the constants code contained within alquimia
as of right now it's only needed to hold the size of the string, it could be
eliminated eventually and this code could be removed along with the
corresponding code file*/

#ifndef BGC_CONSTANTS_H_
#define BGC_CONSTANTS_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* String lengths */
extern const int kBGCMaxStringLength;
extern const int kBGCMaxWordLength;

/* Geochemistry Engine Strings */
extern const char* kBGCStringPFloTran;
extern const char* kBGCStringCrunchFlow;
extern const char* kBGCStringTotalAqueous;
extern const char* kBGCStringTotalSorbed;
extern const char* kBGCStringTotalAqueousPlusSorbed;
extern const char* kBGCStringFree;
extern const char* kBGCStringPH;
extern const char* kBGCStringMineral;
extern const char* kBGCStringGas;
extern const char* kBGCStringCharge;

/* Error Codes */
extern const int kBGCNoError;
extern const int kBGCErrorInvalidEngine;
extern const int kBGCErrorUnknownConstraintName;
extern const int kBGCErrorUnsupportedFunctionality;
extern const int kBGCErrorEngineIntegrity;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif     /* BGC_CONSTANTS_H_ */

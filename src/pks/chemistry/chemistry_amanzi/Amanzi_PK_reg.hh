/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Amanzi_PK registration
*/

#include "Amanzi_PK_ATS.hh"

namespace Amanzi {
namespace AmanziChemistry {

RegisteredPKFactory<Amanzi_PK_ATS> Amanzi_PK_ATS::reg_("chemistry amanzi");

}  // namespace AmanziChemistry
}  // namespace Amanzi

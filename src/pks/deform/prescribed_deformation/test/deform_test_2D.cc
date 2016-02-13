#include <cstdlib>
#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_MSTK.hh"

#include "../deform.hh"

/* **************************************************************** */
TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Deform;

  // first message
  std::cout << "Test: DeformMesh using a simple quadrilateral mesh" << endl;

  // sequential/parallel communicator object
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // read parameter list from input XML file
  Teuchos::Ptr<Teuchos::ParameterList> plist;
  string xml_fname = "test/simple_mesh_deform_test_2D.xml";
  updateParametersFromXmlFile(xml_fname, plist);

  // intantiate the mesh: get the region list
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  
  // intantiate the mesh: build the geometric model
  const int dim2D = 2 ;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(dim2D, region_list, (Epetra_MpiComm *)comm));
  
  // intantiate the mesh: build the geometric model
  int nx(100), ny(100) ;
  Teuchos::RCP<Mesh> mesh = 
    Teuchos::rcp(new Mesh_MSTK( 0., 0., 1., 1., nx, ny, comm, gm));

  // create and initialize the deform mesh class
  DeformMesh deform_test(*plist,mesh);
  deform_test.bell_shaped_profile();

  // say goodbye and exit
  deform_test.print_goodbye();
}

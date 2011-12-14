#ifndef _VIS_HPP_
#define _VIS_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Mesh.hh"
#include "hdf5mpi_mesh.hh"

namespace Amanzi {

  class Vis {

  public:

    Vis(Teuchos::ParameterList& plist, Epetra_MpiComm *comm); 
    Vis(); // created with this constructor this object will not create any output 
    ~Vis();
   
    void create_files(const Amanzi::AmanziMesh::Mesh& mesh);
    void create_timestep(const double& time, const int& cycle);
    void finalize_timestep();
    const bool dump_requested(int cycle);
    void write_vector(const Epetra_MultiVector& vec,
                      const std::vector<std::string>& names);
    void write_vector(const Epetra_Vector& vec, std::string name);
    const bool is_disabled();

  private:    
    void read_parameters(Teuchos::ParameterList& plist);
    
    std::string filebasename; 
    Teuchos::ParameterList plist;
    
    int interval;
    int start;
    int end;
    Teuchos::Array<int> steps;

    Amanzi::HDF5_MPI *viz_output; 

    // disable visualization dumps alltogether
    bool disabled;

    // the Epetra communicator
    Epetra_MpiComm *comm;
  };

}
#endif

#ifndef PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_

#include "CompositeVector.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf,
                        const Teuchos::Ptr<CompositeVector>& sub);

void
CopySubsurfaceToSurface(const CompositeVector& sub,
                        const Teuchos::Ptr<CompositeVector>& surf);


} // namespace


#endif


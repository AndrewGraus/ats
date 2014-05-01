/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   PK factory for self-registering PKs.

   See a more thorough factory discussion in $ATS_DIR/src/factory/factory.hh.

   Simplest usage:

   // pk_implementation.hh
   #include "pk.hh"
   class DerivedPK : public PK {
     DerivedPK(Teuchos::ParameterList& plist,
               const Teuchos::RCP<TreeVector>& solution);
     ...
   private:
     static RegisteredPKFactory<PK,DerivedPK> factory_; // my factory entry
     ...
   };

   ------------------------------------------------------------------------- */

#ifndef ATS_PK_FACTORY_HH_
#define ATS_PK_FACTORY_HH_

#include <iostream>
#include <map>
#include <string>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "TreeVector.hh"
#include "PK.hh"

namespace Amanzi {

class PKFactory {

public:
  typedef std::map<std::string,
                   PK* (*)(const Teuchos::RCP<Teuchos::ParameterList>&,
                           Teuchos::ParameterList&,
                           const Teuchos::RCP<TreeVector>&)> map_type;

  static Teuchos::RCP<PK> CreatePK(const Teuchos::RCP<Teuchos::ParameterList>& plist,
          Teuchos::ParameterList& FElist,
          const Teuchos::RCP<TreeVector>& soln) {
    std::string s = plist->get<std::string>("PK type");
    map_type::iterator iter = GetMap()->find(s);
    if (iter == GetMap()->end()) {
      std::stringstream errmsg;
      errmsg << "PK Factory: cannot find PK type \"" << s << "\"";
      for (map_type::iterator iter=GetMap()->begin();
           iter!=GetMap()->end(); ++iter) {
        errmsg << std::endl << "  option: " << iter->first;
      }
      Errors::Message message(errmsg.str());
      Exceptions::amanzi_throw(message);
    }
    return Teuchos::rcp(iter->second(plist, FElist, soln));
  }

protected:
  static map_type* GetMap() {
    if (!map_) {
      map_ = new map_type;
    }
    return map_;
  }

private:
  static map_type* map_;
};


template<typename T> PK* CreateT(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& soln) {
  return new T(plist, FElist, soln);
}


template<typename T>
class RegisteredPKFactory : public PKFactory {
public:
  // Constructor for the registered factory.  Needs some error checking in
  // case a name s is already in the map? (i.e. two implementations trying to
  // call themselves the same thing) --etc
  RegisteredPKFactory(const std::string& s) {
    GetMap()->insert(std::pair<std::string, PK* (*)(const Teuchos::RCP<Teuchos::ParameterList>&,
            Teuchos::ParameterList&,
            const Teuchos::RCP<TreeVector>&)>(s, &CreateT<T>));
  }
};

} // namespace

#endif

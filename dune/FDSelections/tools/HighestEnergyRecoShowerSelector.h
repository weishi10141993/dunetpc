#ifndef HIGHESTENERGYSHOWERSELECTOR_H_SEEN
#define HIGHESTENERGYSHOWERSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"
//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class HighestEnergyRecoShowerSelector : public RecoShowerSelector{
    public:
      explicit HighestEnergyRecoShowerSelector(fhicl::ParameterSet const& ps) 
        :
        fShowerModuleLabel(ps.get< std::string> ("ModuleLabels.ShowerModuleLabel")),
        fShowerEnergyAlg(ps.get<fhicl::ParameterSet>("ShowerEnergyAlg")) {};

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;
      std::string fShowerModuleLabel;
      //shower::ShowerEnergyAlg fShowerEnergyAlg;

  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestEnergyRecoShowerSelector)
#endif

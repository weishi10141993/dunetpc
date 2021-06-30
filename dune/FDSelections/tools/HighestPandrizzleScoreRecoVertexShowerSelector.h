#ifndef HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTOR_H_SEEN
#define HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "dune/FDSelections/pandrizzle/PandrizzleAlg.h"

//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class HighestPandrizzleScoreRecoVertexShowerSelector : public RecoShowerSelector{
    public:
      explicit HighestPandrizzleScoreRecoVertexShowerSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;

      std::string fShowerModuleLabel;
      std::string fPFParticleModuleLabel;
      FDSelection::PandrizzleAlg fPandrizzleAlg;

      // Pandrizzle Stuff
      double fRecoNuVtxX;
      double fRecoNuVtxY;
      double fRecoNuVtxZ;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector)
#endif

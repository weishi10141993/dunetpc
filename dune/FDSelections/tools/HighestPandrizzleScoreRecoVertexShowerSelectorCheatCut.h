#ifndef HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTORCHEATCUT_H_SEEN
#define HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTORCHEATCUT_H_SEEN

//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "dune/FDSelections/pandrizzle/PandrizzleAlg.h"

//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut : public RecoShowerSelector{
    public:
      explicit HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;

      std::string fHitsModuleLabel;
      std::string fClusterModuleLabel;
      std::string fShowerModuleLabel;
      std::string fPFParticleModuleLabel;

      double fShowerCompletenessCut;
      double fShowerPurityCut;

      FDSelection::PandrizzleAlg fPandrizzleAlg;

      // Pandrizzle Stuff
      double fRecoNuVtxX;
      double fRecoNuVtxY;
      double fRecoNuVtxZ;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut)
#endif

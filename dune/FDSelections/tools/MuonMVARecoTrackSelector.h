#ifndef MUONMVATRACKSELECTOR_H_SEEN
#define MUONMVATRACKSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
//LARSOFT
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class MuonMVARecoTrackSelector : public RecoTrackSelector{
    public:
      explicit MuonMVARecoTrackSelector(fhicl::ParameterSet const& ps) 
        :
        fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
        fPIDModuleLabel(ps.get< std::string> ("ModuleLabels.PIDModuleLabel")) {};
    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;
      std::string fTrackModuleLabel;
      std::string fPIDModuleLabel;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::MuonMVARecoTrackSelector)
#endif

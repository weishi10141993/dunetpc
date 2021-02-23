#ifndef HIGHESTCTPMUONSCOREVERTEXTRACKSELECTOR_H_SEEN
#define HIGHESTCTPMUONSCOREVERTEXTRACKSELECTOR_H_SEEN
//STL
#include <iostream>
#include <limits>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
//LArSoft
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
//DUNE
#include "dune/TrackPID/algorithms/CTPHelper.h"
#include "dune/TrackPID/products/CTPResult.h"
//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class HighestCTPMuonScoreRecoVertexTrackSelector : public RecoTrackSelector{
    public:
      explicit HighestCTPMuonScoreRecoVertexTrackSelector(fhicl::ParameterSet const& ps) 
        :
        fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
        fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
        fConvTrackPID(ps.get<fhicl::ParameterSet>("ctpHelper")) {};


    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;
      std::string fTrackModuleLabel;
      std::string fPFParticleModuleLabel;
      ctp::CTPHelper fConvTrackPID;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestCTPMuonScoreRecoVertexTrackSelector)
#endif

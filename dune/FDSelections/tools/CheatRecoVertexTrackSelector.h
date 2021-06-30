#ifndef CHEATRECOVERTEXTRACKSELECTOR_H_SEEN
#define CHEATRECOVERTEXTRACKSELECTOR_H_SEEN
//STL
#include <iostream>
#include <limits>
//ROOT
#include "TMVA/Reader.h"
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
#include "dune/FDSelections/pandizzle/PandizzleAlg.h"

//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class CheatRecoVertexTrackSelector : public RecoTrackSelector{
    public:
      explicit CheatRecoVertexTrackSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;

      std::string fNuGenModuleLabel;
      std::string fTrackModuleLabel;
      std::string fPFParticleModuleLabel;
      std::string fPandizzleWeightFileName;

      FDSelection::PandizzleAlg fPandizzleAlg;
      TMVA::Reader fPandizzleReader;

      float fTMVAPFPMichelNHits;
      float fTMVAPFPMichelElectronMVA;
      float fTMVAPFPMichelRecoEnergyPlane2;
      float fTMVAPFPTrackDeflecAngleSD;
      float fTMVAPFPTrackLength;
      float fTMVAPFPTrackEvalRatio;
      float fTMVAPFPTrackConcentration;
      float fTMVAPFPTrackCoreHaloRatio;
      float fTMVAPFPTrackConicalness;
      float fTMVAPFPTrackdEdxStart;
      float fTMVAPFPTrackdEdxEnd;
      float fTMVAPFPTrackdEdxEndRatio;
      float fTMVAPFPTrackPIDA;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::CheatRecoVertexTrackSelector)
#endif

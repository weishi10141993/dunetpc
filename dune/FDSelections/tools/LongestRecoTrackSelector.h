#ifndef LONGESTTRACKSELECTOR_H_SEEN
#define LONGESTTRACKSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class LongestRecoTrackSelector : public RecoTrackSelector{
    public:
      explicit LongestRecoTrackSelector(fhicl::ParameterSet const& ps) 
        :
        fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")) {};
    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;
      std::string fTrackModuleLabel;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::LongestRecoTrackSelector)
#endif

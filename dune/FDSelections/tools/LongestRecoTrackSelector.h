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
      explicit LongestRecoTrackSelector(fhicl::ParameterSet const& ps) {};
      void Test();
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::LongestRecoTrackSelector)
#endif

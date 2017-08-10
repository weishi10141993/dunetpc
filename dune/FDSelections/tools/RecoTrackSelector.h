#ifndef TRACKSELECTOR_H_SEEN
#define TRACKSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//LARSOFT

//Base class for selecting reconstructed tracks in the DUNE far detector

namespace FDSelectionTools{
  class RecoTrackSelector {
    public:
      virtual ~RecoTrackSelector() noexcept = default;
      void Test() { std::cout<< "it runs!"<<std::endl; };
  };
}
#endif

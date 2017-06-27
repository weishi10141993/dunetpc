////////////////////////////////////////////////////////////////////////
// Class:       NumuCutSelection
// Plugin Type: analyzer (art v2_07_03)
// File:        NumuCutSelection_module.cc
//
// Generated at Tue Jun 27 06:07:56 2017 by Dominic Brailsford using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace FDSelections {
  class NumuCutSelection;
}


class FDSelections::NumuCutSelection : public art::EDAnalyzer {
public:
  explicit NumuCutSelection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NumuCutSelection(NumuCutSelection const &) = delete;
  NumuCutSelection(NumuCutSelection &&) = delete;
  NumuCutSelection & operator = (NumuCutSelection const &) = delete;
  NumuCutSelection & operator = (NumuCutSelection &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;

private:

  // Declare member data here.

};


FDSelections::NumuCutSelection::NumuCutSelection(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void FDSelections::NumuCutSelection::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  std::cout<<"Running"<<std::endl;
  std::cout<<"Running"<<std::endl;

}

void FDSelections::NumuCutSelection::beginJob()
{
  // Implementation of optional member function here.
}

void FDSelections::NumuCutSelection::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void FDSelections::NumuCutSelection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FDSelections::NumuCutSelection)

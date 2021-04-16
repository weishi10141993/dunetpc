///////////////////////////////////////////////
// PandrizzleAlg.cxx
//
// Reco and true PID stuff up
// D Brailsford & M Wallbank, June 2017
///////////////////////////////////////////////

#include "PandrizzleAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


FDSelection::PandrizzleAlg::PandrizzleAlg(const fhicl::ParameterSet& pset) //:
  //fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg"))
{
}


void FDSelection::PandrizzleAlg::Run(const art::Event& evt) {
}

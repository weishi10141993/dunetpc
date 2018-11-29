///////////////////////////////////////////////
// PandizzleAlg.cxx
//
// Reco and true PID stuff up
// D Brailsford & M Wallbank, June 2017
///////////////////////////////////////////////

#include "PandizzleAlg.h"

FDSelection::PandizzleAlg::PandizzleAlg(const fhicl::ParameterSet& pset) {
  fTrackModuleLabel  = pset.get<std::string>("TrackModuleLabel");
  fShowerModuleLabel = pset.get<std::string>("ShowerModuleLabel");
  fPIDModuleLabel    = pset.get<std::string>("PIDModuleLabel");
  fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
  InitialiseTrees();
}

void FDSelection::PandizzleAlg::InitialiseTrees() {
  fSignalTree = tfs->make<TTree>("DizzleSigTree","Pandizzle Signal Tree");
  fBackgroundTree = tfs->make<TTree>("DizzleBgTree","Pandizzle Background Tree");
  //I am lazy.  Sue me.
  std::map<std::string, TTree*> treeMap;
  treeMap["signal"] = fSignalTree;
  treeMap["background"] = fBackgroundTree;
  for (std::map<std::string, TTree*>::iterator mapIt = treeMap.begin(); mapIt != treeMap.end(); mapIt++){
    TTree *tree = mapIt->second;
    tree->Branch("Run",         &fRun);
    tree->Branch("SubRun",      &fSubRun);
    tree->Branch("Event",      &fEvent);
  }
}

void FDSelection::PandizzleAlg::Run(const art::Event& evt) {

  //Grab the run and subrun info before doing anything
  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();

  //anab::MVAPIDResult res;
  //res.evalRatio = 90.;

  //mf::LogWarning("PandizzleAlg") << "No tracks";

  return;

}



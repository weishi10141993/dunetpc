#ifndef PANDIZZLEALG_H_SEEN
#define PANDIZZLEALG_H_SEEN

///////////////////////////////////////////////
// PandizzleAlg.h
//
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"

//constexpr int kMaxObjects = 999;

namespace FDSelection {
  class PandizzleAlg;
}

class FDSelection::PandizzleAlg {
 public:

  PandizzleAlg(const fhicl::ParameterSet& pset);

  void Run(const art::Event& evt);

 private:

  /*
  /// Returns the true track ID associated with this hit (if more than one, returns the one with highest energy)
  int ParticleID(const art::Ptr<recob::Hit>& hit);

  /// Returns the true particles associated with this object
  std::map<int,double> TrueParticles(const std::vector<art::Ptr<recob::Hit> >& hits);

  /// Returns the true particle most likely associated with this object
  int TrueParticle(const std::vector<art::Ptr<recob::Hit> >& hits);
  */

  /// Initialise the tree
  void InitialiseTrees();

  // module labels
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPIDModuleLabel;
  std::string fPFParticleModuleLabel;

  // tree
  TTree* fSignalTree;
  TTree *fBackgroundTree;
  int fRun;
  int fSubRun;
  int fEvent;
  /*
  int fNObjects;
  bool fTrack[kMaxObjects], fShower[kMaxObjects];
  int fTruePDG[kMaxObjects], fRecoPDG[kMaxObjects];
  double fElectronMVA[kMaxObjects], fMuonMVA[kMaxObjects], fPhotonMVA[kMaxObjects], fProtonMVA[kMaxObjects], fPionMVA[kMaxObjects];
  double fObjPurity[kMaxObjects];
  bool fPrimary[kMaxObjects];
  double fTrueEnergy[kMaxObjects], fRecoEnergy[kMaxObjects];
  double fTrueX[kMaxObjects], fTrueY[kMaxObjects], fTrueZ[kMaxObjects];
  double fTrueEndX[kMaxObjects], fTrueEndY[kMaxObjects], fTrueEndZ[kMaxObjects];
  double fRecoX[kMaxObjects], fRecoY[kMaxObjects], fRecoZ[kMaxObjects];
  double fRecoEndX[kMaxObjects], fRecoEndY[kMaxObjects], fRecoEndZ[kMaxObjects];
  double fRecoLength[kMaxObjects];
  int fRecoPoints[kMaxObjects];
  */

  // services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<art::TFileService> tfs;

};

#endif

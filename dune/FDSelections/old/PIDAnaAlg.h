#ifndef PIDANAALG_H_SEEN
#define PIDANAALG_H_SEEN

///////////////////////////////////////////////
// PIDAnaAlg.h
//
// Reco and true PID stuff up
// D Brailsford & M Wallbank, June 2017
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
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

#include "PandizzleAlg.h"

constexpr int kMaxObjects = 999;

namespace FDSelection {
  class PIDAnaAlg;
}

class FDSelection::PIDAnaAlg {
 public:

  PIDAnaAlg(const fhicl::ParameterSet& pset);

  void Run(const art::Event& evt);

 private:

  /// Returns the true track ID associated with this hit (if more than one, returns the one with highest energy)
  int ParticleID(const art::Ptr<recob::Hit>& hit);

  /// Returns the true particles associated with this object
  std::map<int,double> TrueParticles(const std::vector<art::Ptr<recob::Hit> >& hits);

  /// Returns the true particle most likely associated with this object
  int TrueParticle(const std::vector<art::Ptr<recob::Hit> >& hits);

  /// Initialise the tree
  void InitialiseTree();

  // module labels
  std::string fTrackModuleLabel, fShowerModuleLabel, fPIDModuleLabel;

  // tree
  TTree* fTree;
  int fRun;
  int fSubRun;
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

  // services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<art::TFileService> tfs;

};

#endif

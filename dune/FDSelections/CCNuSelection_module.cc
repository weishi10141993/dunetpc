////////////////////////////////////////////////////////////////////////
// Class:       CCNuSelection
// Plugin Type: analyzer (art v2_07_03)
// File:        CCNuSelection_module.cc
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

//STL
#include <iostream>
//ROOT
#include "TTree.h"
#include "TMVA/Reader.h"

//ART
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/make_tool.h" 
//LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"


//DUNE
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

//Custom
//#include "PIDAnaAlg.h"
#include "PandizzleAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"
#include "tools/RecoShowerSelector.h"


constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

namespace FDSelection {
  class CCNuSelection;
}


class FDSelection::CCNuSelection : public art::EDAnalyzer {
public:
  explicit CCNuSelection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CCNuSelection(CCNuSelection const &) = delete;
  CCNuSelection(CCNuSelection &&) = delete;
  CCNuSelection & operator = (CCNuSelection const &) = delete;
  CCNuSelection & operator = (CCNuSelection &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginSubRun(art::SubRun const & sr) override;
  void endSubRun(art::SubRun const & sr) override;
  void endJob() override;

private:

  //Delcare private functions
  void Reset();      //Resets all tree vars
  void GetEventInfo(art::Event const & evt); //Grab event-level info that is necessary/handy for selecting events but doesn't really fall under the 'selection' banner
  void GetTruthInfo(art::Event const & evt);  //Grab the truth info from the art record
  void FillVertexInformation(art::Event const & evt); //Grab the vertex parameters 
  void RunTrackSelection(art::Event const & evt);  //Run the selection and dump relevant info to the truee
  double CalculateTrackCharge(art::Ptr<recob::Track> const track, std::vector< art::Ptr< recob::Hit> > const track_hits); //Calculate hit charge on track as measured by the collection plane
  bool IsTrackContained(art::Ptr<recob::Track> const track, std::vector< art::Ptr<recob::Hit > > const track_hits, art::Event const & evt); // check if the track is contained in the detector
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt);
  art::Ptr<recob::Track> GetTrackMatchedPFParticle(art::Ptr<recob::PFParticle> const pfparticle, art::Event const & evt);

  void FillChildPFPInformation(art::Ptr<recob::Track> const track, art::Event const & evt); //Grab the vertex parameters associated with this track

  //Shower stuff
  void RunShowerSelection(art::Event const & evt);  //Run the selection and dump relevant info to the truee
  double CalculateShowerCharge(art::Ptr<recob::Shower> const track, std::vector< art::Ptr< recob::Hit> > const track_hits); //Calculate hit charge on track as measured by the collection plane
  bool IsPFParticlePrimaryDaughter(art::Ptr<recob::PFParticle> const & pfparticle, std::vector<art::Ptr<recob::PFParticle > > const & pfparticles);



  // Declare member data here.

  //Algs
  //PIDAnaAlg fPIDAnaAlg;
  PandizzleAlg fPandizzleAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  shower::ShowerEnergyAlg fShowerEnergyAlg;


  //Tools
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;
  std::unique_ptr<FDSelectionTools::RecoShowerSelector> fRecoShowerSelector;


  TTree *fTree; //The selection tree
  //Generic stuff
  int fRun;
  int fSubRun;
  int fEvent;
  int fIsMC;
  //Detector stuff
  double fT0;
  //Neutrino stuff
  int fNuPdg; //Interaction PDG
  int fBeamPdg; //PDG at point of creation
  int fNC;    // 1=is NC, 0=otherwise
  int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  int fTargetZ; //Atomic number of scattering target
  double fQ2; 
  double fENu; 
  double fW; //X-Sec params
  double fX;
  double fY;
  double fNuMomX; //Neutrino momentums
  double fNuMomY;
  double fNuMomZ;
  double fNuMomT;
  double fNuX; //Interaction positions
  double fNuY;
  double fNuZ;
  double fNuT;
  //Outgoing particle count
  int fNPiP; //Number of pi+
  int fNPim; //Number of pi-
  int fNPi0; //Number of pi0
  int fNPhotons; //Number of photons
  int fNProtons; //Number of protons
  int fNNeutrons; //Number of neutrinos
  int fNOther; //Number of other particles
  //Outgoing lepton stuff
  int fLepPDG;
  double fMomLepX;
  double fMomLepY;
  double fMomLepZ;
  double fMomLepT;
  double fLepEndX;
  double fLepEndY;
  double fLepEndZ;
  double fLepEndT;
  double fLepNuAngle;
  //Selection stuff
  //true bits
  int fSelTrackTruePDG;
  int fSelTrackTruePrimary;
  double fSelTrackTrueMomX;
  double fSelTrackTrueMomY;
  double fSelTrackTrueMomZ;
  double fSelTrackTrueMomT;
  double fSelTrackTrueStartX;
  double fSelTrackTrueStartY;
  double fSelTrackTrueStartZ;
  double fSelTrackTrueStartT;
  double fSelTrackTrueEndX;
  double fSelTrackTrueEndY;
  double fSelTrackTrueEndZ;
  double fSelTrackTrueEndT;
  //reco bits
  int fSelTrackRecoNHits;
  double fSelTrackRecoCompleteness;
  double fSelTrackRecoHitPurity;
  double fSelTrackRecoMomX;
  double fSelTrackRecoMomY;
  double fSelTrackRecoMomZ;
  double fSelTrackRecoMomT;
  double fSelTrackRecoStartX;
  double fSelTrackRecoStartY;
  double fSelTrackRecoStartZ;
  double fSelTrackRecoStartT;
  double fSelTrackRecoEndX;
  double fSelTrackRecoEndY;
  double fSelTrackRecoEndZ;
  double fSelTrackRecoEndT;
  double fSelTrackRecoUpstreamX;
  double fSelTrackRecoUpstreamY;
  double fSelTrackRecoUpstreamZ;
  double fSelTrackRecoUpstreamT;
  double fSelTrackRecoDownstreamX;
  double fSelTrackRecoDownstreamY;
  double fSelTrackRecoDownstreamZ;
  double fSelTrackRecoDownstreamT;
  double fSelTrackRecoEndClosestToVertexX;
  double fSelTrackRecoEndClosestToVertexY;
  double fSelTrackRecoEndClosestToVertexZ;
  double fSelTrackRecoLength;
  bool   fSelTrackRecoContained;
  double fSelTrackRecoCharge;
  double fSelTrackRecoMomMCS;
  double fSelTrackRecoMomContained;
  //int fSelTrackRecoNMatchedVertices;
  double fSelTrackRecoVertexX;
  double fSelTrackRecoVertexY;
  double fSelTrackRecoVertexZ;
  int fSelTrackRecoNChildPFP;
  int fSelTrackRecoNChildTrackPFP;
  int fSelTrackRecoNChildShowerPFP;
  //MVA bits
  double fSelTrackMVAElectron;
  double fSelTrackMVAPion;
  double fSelTrackMVAMuon;
  double fSelTrackMVAProton;
  double fSelTrackMVAPhoton;
  //Pandizzle var
  double fSelTrackPandizzleVar;
  //Event-level bits
  double fRecoEventCharge; //Total collected charge (as measured by the collection planes)
  //reco energy bits
  double fNumuRecoMomLep; //Reco lepton momentum
  double fNumuRecoEHad; //Reco hadronic energy
  double fNumuRecoENu; //Reco neutrino energy

  //Reco neutrino bits
  double fRecoNuVtxX;
  double fRecoNuVtxY;
  double fRecoNuVtxZ;
  int fRecoNuVtxNShowers;
  int fRecoNuVtxNTracks;
  int fRecoNuVtxNChildren;

  //Shower Selection stuff
  //true bits
  int fSelShowerTruePDG;
  int fSelShowerTruePrimary;
  double fSelShowerTrueMomX;
  double fSelShowerTrueMomY;
  double fSelShowerTrueMomZ;
  double fSelShowerTrueMomT;
  double fSelShowerTrueStartX;
  double fSelShowerTrueStartY;
  double fSelShowerTrueStartZ;
  double fSelShowerTrueStartT;
  double fSelShowerTrueEndX;
  double fSelShowerTrueEndY;
  double fSelShowerTrueEndZ;
  double fSelShowerTrueEndT;
  //reco bits
  int fSelShowerRecoNHits;
  double fSelShowerRecoCompleteness;
  double fSelShowerRecoHitPurity;
  double fSelShowerRecoMomX;
  double fSelShowerRecoMomY;
  double fSelShowerRecoMomZ;
  double fSelShowerRecoMomT;
  double fSelShowerRecoStartX;
  double fSelShowerRecoStartY;
  double fSelShowerRecoStartZ;
  double fSelShowerRecoStartT;
  double fSelShowerRecoEndX;
  double fSelShowerRecoEndY;
  double fSelShowerRecoEndZ;
  double fSelShowerRecoEndT;
  double fSelShowerRecoUpstreamX;
  double fSelShowerRecoUpstreamY;
  double fSelShowerRecoUpstreamZ;
  double fSelShowerRecoUpstreamT;
  double fSelShowerRecoDownstreamX;
  double fSelShowerRecoDownstreamY;
  double fSelShowerRecoDownstreamZ;
  double fSelShowerRecoDownstreamT;
  double fSelShowerRecoEndClosestToVertexX;
  double fSelShowerRecoEndClosestToVertexY;
  double fSelShowerRecoEndClosestToVertexZ;
  double fSelShowerRecoLength;
  bool   fSelShowerRecoContained;
  double fSelShowerRecoCharge;
  double fSelShowerRecoMomMCS;
  double fSelShowerRecoMomContained;
  int fSelShowerRecoNMatchedVertices;
  double fSelShowerRecoVertexX;
  double fSelShowerRecoVertexY;
  double fSelShowerRecoVertexZ;
  int fSelShowerRecoNChildPFP;
  int fSelShowerRecoNChildTrackPFP;
  int fSelShowerRecoNChildShowerPFP;
  double fSelShowerRecoDirX;
  double fSelShowerRecoDirY;
  double fSelShowerRecoDirZ;
  double fSelShowerRecodEdx;
  bool fSelShowerRecoIsPrimaryPFPDaughter;


  //MVA bits
  double fSelShowerMVAElectron;
  double fSelShowerMVAPion;
  double fSelShowerMVAMuon;
  double fSelShowerMVAProton;
  double fSelShowerMVAPhoton;
  //reco energy bits
  double fNueRecoMomLep; //Reco lepton momentum
  double fNueRecoEHad; //Reco hadronic energy
  double fNueRecoENu; //Reco neutrino energy



  //POT tree stuff
  TTree* fPOTTree;
  double fPOT;

  //Fhicl pset labels
  std::string fNuGenModuleLabel;
  std::string fLargeantModuleLabel;
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fVertexModuleLabel;
  std::string fPIDModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fHitsModuleLabel;
  std::string fPOTModuleLabel;
  std::string fNumuEnergyRecoModuleLabel;
  std::string fNueEnergyRecoModuleLabel;

 
  //Processing flags
  bool fUsePandoraVertex;


  //TMVAStuff
  TMVA::Reader fReader;
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


FDSelection::CCNuSelection::CCNuSelection(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)   ,
  //fPIDAnaAlg(pset.get<fhicl::ParameterSet>("ModuleLabels"))   ,
  fPandizzleAlg(pset) ,
  fCalorimetryAlg          (pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(pset.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))},
  fNuGenModuleLabel        (pset.get< std::string >("ModuleLabels.NuGenModuleLabel")),
  fLargeantModuleLabel     (pset.get< std::string >("ModuleLabels.LargeantModuleLabel")),
  fTrackModuleLabel        (pset.get< std::string >("ModuleLabels.TrackModuleLabel")),
  fShowerModuleLabel        (pset.get< std::string >("ModuleLabels.ShowerModuleLabel")),
  fPFParticleModuleLabel   (pset.get< std::string >("ModuleLabels.PFParticleModuleLabel")),
  fVertexModuleLabel   (pset.get< std::string >("ModuleLabels.VertexModuleLabel")),
  fPIDModuleLabel          (pset.get< std::string >("ModuleLabels.PIDModuleLabel")),
  fParticleIDModuleLabel          (pset.get< std::string >("ModuleLabels.ParticleIDModuleLabel")),
  fHitsModuleLabel         (pset.get< std::string >("ModuleLabels.HitsModuleLabel")),
  fPOTModuleLabel          (pset.get< std::string >("ModuleLabels.POTModuleLabel")),
  fNumuEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.NumuEnergyRecoModuleLabel")),
  fNueEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.NueEnergyRecoModuleLabel")),
  fUsePandoraVertex        (pset.get< bool >("UsePandoraVertex")),
  fReader("") {
  fReader.AddVariable("PFPMichelNHits",&fTMVAPFPMichelNHits);
  fReader.AddVariable("PFPMichelElectronMVA",&fTMVAPFPMichelElectronMVA);
  fReader.AddVariable("PFPMichelRecoEnergyPlane2",&fTMVAPFPMichelRecoEnergyPlane2);
  fReader.AddVariable("PFPTrackDeflecAngleSD",&fTMVAPFPTrackDeflecAngleSD);
  fReader.AddVariable("PFPTrackLength",&fTMVAPFPTrackLength);
  fReader.AddVariable("PFPTrackEvalRatio",&fTMVAPFPTrackEvalRatio);
  fReader.AddVariable("PFPTrackConcentration",&fTMVAPFPTrackConcentration);
  fReader.AddVariable("PFPTrackCoreHaloRatio",&fTMVAPFPTrackCoreHaloRatio);
  fReader.AddVariable("PFPTrackConicalness",&fTMVAPFPTrackConicalness);
  fReader.AddVariable("PFPTrackdEdxStart",&fTMVAPFPTrackdEdxStart);
  fReader.AddVariable("PFPTrackdEdxEnd",&fTMVAPFPTrackdEdxEnd);
  fReader.AddVariable("PFPTrackdEdxEndRatio",&fTMVAPFPTrackdEdxEndRatio);
  //fReader.AddVariable("PFPTrackPIDA",&fTMVAPFPTrackPIDA);
  fReader.BookMVA("BDTG","/dune/app/users/dbrailsf/oscillation/nu_mu/cutsel/trainings/pandizzle/dataset_pandizzle/weights/TMVAClassification_BDTG.weights.xml");
}

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
  //Get the generic stuff that can be pulled from the top of the record
  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();
  fIsMC = !(evt.isRealData());

  GetEventInfo(evt);
  if (fIsMC) GetTruthInfo(evt);
  FillVertexInformation(evt);
  RunTrackSelection(evt);
  RunShowerSelection(evt);


  //fPIDAnaAlg.Run(evt);
  //fPandizzleAlg.Run(evt);

  fTree->Fill();
  Reset(); //Reset at the end of the event
}

void FDSelection::CCNuSelection::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("ccnusel","CC nu selection");
    fTree->Branch("Run",&fRun);
    fTree->Branch("SubRun",&fSubRun);
    fTree->Branch("Event",&fEvent);
    fTree->Branch("IsMC",&fIsMC);
    fTree->Branch("T0",&fT0);
    fTree->Branch("NuPdg",&fNuPdg);
    fTree->Branch("BeamPdg",&fBeamPdg);
    fTree->Branch("NC",&fNC);
    fTree->Branch("Mode",&fMode);
    fTree->Branch("TargetZ",&fTargetZ);
    fTree->Branch("Q2",&fQ2);
    fTree->Branch("Enu",&fENu);
    fTree->Branch("W",&fW);
    fTree->Branch("X",&fX);
    fTree->Branch("Y",&fY);
    fTree->Branch("NuMomX",&fNuMomX);
    fTree->Branch("NuMomY",&fNuMomY);
    fTree->Branch("NuMomZ",&fNuMomZ);
    fTree->Branch("NuMomT",&fNuMomT);
    fTree->Branch("NuX",&fNuX);
    fTree->Branch("NuY",&fNuY);
    fTree->Branch("NuZ",&fNuZ);
    fTree->Branch("NuT",&fNuT);
    fTree->Branch("NPiP",&fNPiP);
    fTree->Branch("NPim",&fNPim);
    fTree->Branch("NPi0",&fNPi0);
    fTree->Branch("NPhotons",&fNPhotons);
    fTree->Branch("NProton",&fNProtons);
    fTree->Branch("NNeutrons",&fNNeutrons);
    fTree->Branch("NOther",&fNOther);
    fTree->Branch("LepPDG",&fLepPDG);
    fTree->Branch("MomLepX",&fMomLepX);
    fTree->Branch("MomLepY",&fMomLepY);
    fTree->Branch("MomLepZ",&fMomLepZ);
    fTree->Branch("MomLepT",&fMomLepT);
    fTree->Branch("LepEndX",&fLepEndX);
    fTree->Branch("LepEndY",&fLepEndY);
    fTree->Branch("LepEndZ",&fLepEndZ);
    fTree->Branch("LepEndT",&fLepEndT);
    fTree->Branch("LepNuAngle",&fLepNuAngle);
    fTree->Branch("SelTrackTruePDG",&fSelTrackTruePDG);
    fTree->Branch("SelTrackTruePrimary",&fSelTrackTruePrimary);
    fTree->Branch("SelTrackTrueMomX",&fSelTrackTrueMomX);
    fTree->Branch("SelTrackTrueMomY",&fSelTrackTrueMomY);
    fTree->Branch("SelTrackTrueMomZ",&fSelTrackTrueMomZ);
    fTree->Branch("SelTrackTrueMomT",&fSelTrackTrueMomT);
    fTree->Branch("SelTrackTrueStartX",&fSelTrackTrueStartX);
    fTree->Branch("SelTrackTrueStartY",&fSelTrackTrueStartY);
    fTree->Branch("SelTrackTrueStartZ",&fSelTrackTrueStartZ);
    fTree->Branch("SelTrackTrueStartT",&fSelTrackTrueStartT);
    fTree->Branch("SelTrackTrueEndX",&fSelTrackTrueEndX);
    fTree->Branch("SelTrackTrueEndY",&fSelTrackTrueEndY);
    fTree->Branch("SelTrackTrueEndZ",&fSelTrackTrueEndZ);
    fTree->Branch("SelTrackTrueEndT",&fSelTrackTrueEndT);
    fTree->Branch("SelTrackRecoNHits",&fSelTrackRecoNHits);
    fTree->Branch("SelTrackRecoCompleteness",&fSelTrackRecoCompleteness);
    fTree->Branch("SelTrackRecoHitPurity",&fSelTrackRecoHitPurity);
    fTree->Branch("SelTrackRecoMomX",&fSelTrackRecoMomX);
    fTree->Branch("SelTrackRecoMomY",&fSelTrackRecoMomY);
    fTree->Branch("SelTrackRecoMomZ",&fSelTrackRecoMomZ);
    fTree->Branch("SelTrackRecoMomT",&fSelTrackRecoMomT);
    fTree->Branch("SelTrackRecoStartX",&fSelTrackRecoStartX);
    fTree->Branch("SelTrackRecoStartY",&fSelTrackRecoStartY);
    fTree->Branch("SelTrackRecoStartZ",&fSelTrackRecoStartZ);
    fTree->Branch("SelTrackRecoStartT",&fSelTrackRecoStartT);
    fTree->Branch("SelTrackRecoEndX",&fSelTrackRecoEndX);
    fTree->Branch("SelTrackRecoEndY",&fSelTrackRecoEndY);
    fTree->Branch("SelTrackRecoEndZ",&fSelTrackRecoEndZ);
    fTree->Branch("SelTrackRecoEndT",&fSelTrackRecoEndT);
    fTree->Branch("SelTrackRecoUpstreamX",&fSelTrackRecoUpstreamX);
    fTree->Branch("SelTrackRecoUpstreamY",&fSelTrackRecoUpstreamY);
    fTree->Branch("SelTrackRecoUpstreamZ",&fSelTrackRecoUpstreamZ);
    fTree->Branch("SelTrackRecoUpstreamT",&fSelTrackRecoUpstreamT);
    fTree->Branch("SelTrackRecoDownstreamX",&fSelTrackRecoDownstreamX);
    fTree->Branch("SelTrackRecoDownstreamY",&fSelTrackRecoDownstreamY);
    fTree->Branch("SelTrackRecoDownstreamZ",&fSelTrackRecoDownstreamZ);
    fTree->Branch("SelTrackRecoDownstreamT",&fSelTrackRecoDownstreamT);
    fTree->Branch("SelTrackRecoEndClosestToVertexX",&fSelTrackRecoEndClosestToVertexX);
    fTree->Branch("SelTrackRecoEndClosestToVertexY",&fSelTrackRecoEndClosestToVertexY);
    fTree->Branch("SelTrackRecoEndClosestToVertexZ",&fSelTrackRecoEndClosestToVertexZ);
    fTree->Branch("SelTrackRecoLength",&fSelTrackRecoLength);
    fTree->Branch("SelTrackRecoContained",&fSelTrackRecoContained);
    fTree->Branch("SelTrackRecoCharge",&fSelTrackRecoCharge);
    fTree->Branch("SelTrackRecoMomMCS",&fSelTrackRecoMomMCS);
    fTree->Branch("SelTrackRecoMomContained",&fSelTrackRecoMomContained);
    //fTree->Branch("SelTrackRecoNMatchedVertices",&fSelTrackRecoNMatchedVertices);
    fTree->Branch("SelTrackRecoVertexX",&fSelTrackRecoVertexX);
    fTree->Branch("SelTrackRecoVertexY",&fSelTrackRecoVertexY);
    fTree->Branch("SelTrackRecoVertexZ",&fSelTrackRecoVertexZ);
    fTree->Branch("SelTrackRecoNChildPFP",&fSelTrackRecoNChildPFP);
    fTree->Branch("SelTrackRecoNChildTrackPFP",&fSelTrackRecoNChildTrackPFP);
    fTree->Branch("SelTrackRecoNChildShowerPFP",&fSelTrackRecoNChildShowerPFP);
    fTree->Branch("SelTrackMVAElectron",&fSelTrackMVAElectron);
    fTree->Branch("SelTrackMVAPion",&fSelTrackMVAPion);
    fTree->Branch("SelTrackMVAMuon",&fSelTrackMVAMuon);
    fTree->Branch("SelTrackMVAProton",&fSelTrackMVAProton);
    fTree->Branch("SelTrackMVAPhoton",&fSelTrackMVAPhoton);
    fTree->Branch("SelTrackPandizzleVar",&fSelTrackPandizzleVar);
    fTree->Branch("TMVAPFPMichelNHits",&fTMVAPFPMichelNHits);
    fTree->Branch("TMVAPFPMichelElectronMVA",&fTMVAPFPMichelElectronMVA);
    fTree->Branch("TMVAPFPMichelRecoEnergyPlane2",&fTMVAPFPMichelRecoEnergyPlane2);
    fTree->Branch("TMVAPFPTrackDeflecAngleSD",&fTMVAPFPTrackDeflecAngleSD);
    fTree->Branch("TMVAPFPTrackLength",&fTMVAPFPTrackLength);
    fTree->Branch("TMVAPFPTrackEvalRatio",&fTMVAPFPTrackEvalRatio);
    fTree->Branch("TMVAPFPTrackConcentration",&fTMVAPFPTrackConcentration);
    fTree->Branch("TMVAPFPTrackCoreHaloRatio",&fTMVAPFPTrackCoreHaloRatio);
    fTree->Branch("TMVAPFPTrackConicalness",&fTMVAPFPTrackConicalness);
    fTree->Branch("TMVAPFPTrackdEdxStart",&fTMVAPFPTrackdEdxStart);
    fTree->Branch("TMVAPFPTrackdEdxEnd",&fTMVAPFPTrackdEdxEnd);
    fTree->Branch("TMVAPFPTrackdEdxEndRatio",&fTMVAPFPTrackdEdxEndRatio);
    fTree->Branch("TMVAPFPTrackPIDA",&fTMVAPFPTrackPIDA);
    fTree->Branch("SelShowerTruePDG",&fSelShowerTruePDG);
    fTree->Branch("SelShowerTruePrimary",&fSelShowerTruePrimary);
    fTree->Branch("SelShowerTrueMomX",&fSelShowerTrueMomX);
    fTree->Branch("SelShowerTrueMomY",&fSelShowerTrueMomY);
    fTree->Branch("SelShowerTrueMomZ",&fSelShowerTrueMomZ);
    fTree->Branch("SelShowerTrueMomT",&fSelShowerTrueMomT);
    fTree->Branch("SelShowerTrueStartX",&fSelShowerTrueStartX);
    fTree->Branch("SelShowerTrueStartY",&fSelShowerTrueStartY);
    fTree->Branch("SelShowerTrueStartZ",&fSelShowerTrueStartZ);
    fTree->Branch("SelShowerTrueStartT",&fSelShowerTrueStartT);
    fTree->Branch("SelShowerTrueEndX",&fSelShowerTrueEndX);
    fTree->Branch("SelShowerTrueEndY",&fSelShowerTrueEndY);
    fTree->Branch("SelShowerTrueEndZ",&fSelShowerTrueEndZ);
    fTree->Branch("SelShowerTrueEndT",&fSelShowerTrueEndT);
    fTree->Branch("SelShowerRecoNHits",&fSelShowerRecoNHits);
    fTree->Branch("SelShowerRecoCompleteness",&fSelShowerRecoCompleteness);
    fTree->Branch("SelShowerRecoHitPurity",&fSelShowerRecoHitPurity);
    fTree->Branch("SelShowerRecoMomX",&fSelShowerRecoMomX);
    fTree->Branch("SelShowerRecoMomY",&fSelShowerRecoMomY);
    fTree->Branch("SelShowerRecoMomZ",&fSelShowerRecoMomZ);
    fTree->Branch("SelShowerRecoMomT",&fSelShowerRecoMomT);
    fTree->Branch("SelShowerRecoStartX",&fSelShowerRecoStartX);
    fTree->Branch("SelShowerRecoStartY",&fSelShowerRecoStartY);
    fTree->Branch("SelShowerRecoStartZ",&fSelShowerRecoStartZ);
    fTree->Branch("SelShowerRecoStartT",&fSelShowerRecoStartT);
    fTree->Branch("SelShowerRecoEndX",&fSelShowerRecoEndX);
    fTree->Branch("SelShowerRecoEndY",&fSelShowerRecoEndY);
    fTree->Branch("SelShowerRecoEndZ",&fSelShowerRecoEndZ);
    fTree->Branch("SelShowerRecoEndT",&fSelShowerRecoEndT);
    fTree->Branch("SelShowerRecoUpstreamX",&fSelShowerRecoUpstreamX);
    fTree->Branch("SelShowerRecoUpstreamY",&fSelShowerRecoUpstreamY);
    fTree->Branch("SelShowerRecoUpstreamZ",&fSelShowerRecoUpstreamZ);
    fTree->Branch("SelShowerRecoUpstreamT",&fSelShowerRecoUpstreamT);
    fTree->Branch("SelShowerRecoDownstreamX",&fSelShowerRecoDownstreamX);
    fTree->Branch("SelShowerRecoDownstreamY",&fSelShowerRecoDownstreamY);
    fTree->Branch("SelShowerRecoDownstreamZ",&fSelShowerRecoDownstreamZ);
    fTree->Branch("SelShowerRecoDownstreamT",&fSelShowerRecoDownstreamT);
    fTree->Branch("SelShowerRecoEndClosestToVertexX",&fSelShowerRecoEndClosestToVertexX);
    fTree->Branch("SelShowerRecoEndClosestToVertexY",&fSelShowerRecoEndClosestToVertexY);
    fTree->Branch("SelShowerRecoEndClosestToVertexZ",&fSelShowerRecoEndClosestToVertexZ);
    fTree->Branch("SelShowerRecoLength",&fSelShowerRecoLength);
    fTree->Branch("SelShowerRecoContained",&fSelShowerRecoContained);
    fTree->Branch("SelShowerRecoCharge",&fSelShowerRecoCharge);
    fTree->Branch("SelShowerRecoMomMCS",&fSelShowerRecoMomMCS);
    fTree->Branch("SelShowerRecoMomContained",&fSelShowerRecoMomContained);
    fTree->Branch("SelShowerRecoNMatchedVertices",&fSelShowerRecoNMatchedVertices);
    fTree->Branch("SelShowerRecoVertexX",&fSelShowerRecoVertexX);
    fTree->Branch("SelShowerRecoVertexY",&fSelShowerRecoVertexY);
    fTree->Branch("SelShowerRecoVertexZ",&fSelShowerRecoVertexZ);
    fTree->Branch("SelShowerRecoNChildPFP",&fSelShowerRecoNChildPFP);
    fTree->Branch("SelShowerRecoNChildTrackPFP",&fSelShowerRecoNChildTrackPFP);
    fTree->Branch("SelShowerRecoNChildShowerPFP",&fSelShowerRecoNChildShowerPFP);
    fTree->Branch("SelShowerRecoDirX",&fSelShowerRecoDirX);
    fTree->Branch("SelShowerRecoDirY",&fSelShowerRecoDirY);
    fTree->Branch("SelShowerRecoDirZ",&fSelShowerRecoDirZ);
    fTree->Branch("SelShowerRecodEdx",&fSelShowerRecodEdx);
    fTree->Branch("SelShowerRecoIsPrimaryPFPDaughter",&fSelShowerRecoIsPrimaryPFPDaughter);

    fTree->Branch("SelShowerMVAElectron",&fSelShowerMVAElectron);
    fTree->Branch("SelShowerMVAPion",&fSelShowerMVAPion);
    fTree->Branch("SelShowerMVAMuon",&fSelShowerMVAMuon);
    fTree->Branch("SelShowerMVAProton",&fSelShowerMVAProton);
    fTree->Branch("SelShowerMVAPhoton",&fSelShowerMVAPhoton);
    fTree->Branch("RecoNuVtxX",&fRecoNuVtxX);
    fTree->Branch("RecoNuVtxY",&fRecoNuVtxY);
    fTree->Branch("RecoNuVtxZ",&fRecoNuVtxZ);
    fTree->Branch("RecoNuVtxNShowers",&fRecoNuVtxNShowers);
    fTree->Branch("RecoNuVtxNTracks",&fRecoNuVtxNTracks);
    fTree->Branch("RecoNuVtxNChildren",&fRecoNuVtxNChildren);

    fTree->Branch("RecoEventCharge",&fRecoEventCharge);
    fTree->Branch("NumuRecoMomLep",&fNumuRecoMomLep);
    fTree->Branch("NumuRecoEHad",&fNumuRecoEHad);
    fTree->Branch("NumuRecoENu",&fNumuRecoENu);
    fTree->Branch("NueRecoMomLep",&fNueRecoMomLep);
    fTree->Branch("NueRecoEHad",&fNueRecoEHad);
    fTree->Branch("NueRecoENu",&fNueRecoENu);



    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("POT",&fPOT);
    fPOTTree->Branch("Run",&fRun);
    fPOTTree->Branch("SubRun",&fSubRun);


    Reset();  //Default value all variables now
}

void FDSelection::CCNuSelection::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}
void FDSelection::CCNuSelection::endSubRun(const art::SubRun& sr){
  //Need the run and subrun
  fRun = sr.run();
  fSubRun = sr.subRun();
  //Need the POT (obvs)
  art::Handle< sumdata::POTSummary > potListHandle;

  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    fPOT = potListHandle->totpot;
  else
    fPOT = 0.;
  if (fPOTTree) fPOTTree->Fill();
}

void FDSelection::CCNuSelection::endJob()
{
  // Implementation of optional member function here.
}

void FDSelection::CCNuSelection::Reset()
{
  //Generic stuff
  fRun = kDefInt;
  fSubRun = kDefInt;
  fEvent = kDefInt;
  fIsMC = kDefInt;
  //Detector stuff
  fT0 = kDefDoub;
  //Neutrino stuff
  fNuPdg = kDefInt; 
  fBeamPdg = kDefInt; 
  fNC = kDefInt;    
  fMode = kDefInt; 
  fTargetZ = kDefInt;
  fQ2 = kDefDoub; 
  fENu = kDefDoub; 
  fW = kDefDoub; 
  fX = kDefDoub;
  fY = kDefDoub;
  fNuMomX = kDefDoub; 
  fNuMomY = kDefDoub;
  fNuMomZ = kDefDoub;
  fNuMomT = kDefDoub;
  fNuX = kDefDoub; 
  fNuY = kDefDoub;
  fNuZ = kDefDoub;
  fNuT = kDefDoub;
  fNPiP = 0;
  fNPim = 0;
  fNPi0 = 0;
  fNPhotons = 0;
  fNProtons = 0;
  fNNeutrons = 0;
  fNOther = 0;
  //Outgoing lepton stuff
  fLepPDG = kDefInt;
  fMomLepX = kDefDoub;
  fMomLepY = kDefDoub;
  fMomLepZ = kDefDoub;
  fMomLepT = kDefDoub;
  fLepEndX = kDefDoub;
  fLepEndY = kDefDoub;
  fLepEndZ = kDefDoub;
  fLepEndT = kDefDoub;
  fLepNuAngle = kDefDoub;
  //Selection stuff
  //Track stuff
  //true bits
  fSelTrackTruePDG = kDefInt;
  fSelTrackTruePrimary = kDefInt;
  fSelTrackTrueMomX = kDefDoub;
  fSelTrackTrueMomY = kDefDoub;
  fSelTrackTrueMomZ = kDefDoub;
  fSelTrackTrueMomT = kDefDoub;
  fSelTrackTrueStartX = kDefDoub;
  fSelTrackTrueStartY = kDefDoub;
  fSelTrackTrueStartZ = kDefDoub;
  fSelTrackTrueStartT = kDefDoub;
  fSelTrackTrueEndX = kDefDoub;
  fSelTrackTrueEndY = kDefDoub;
  fSelTrackTrueEndZ = kDefDoub;
  fSelTrackTrueEndT = kDefDoub;
  //reco bits
  fSelTrackRecoNHits = kDefInt;
  fSelTrackRecoCompleteness = kDefDoub;
  fSelTrackRecoHitPurity = kDefDoub;
  fSelTrackRecoMomX = kDefDoub;
  fSelTrackRecoMomY = kDefDoub;
  fSelTrackRecoMomZ = kDefDoub;
  fSelTrackRecoMomT = kDefDoub;
  fSelTrackRecoStartX = kDefDoub;
  fSelTrackRecoStartY = kDefDoub;
  fSelTrackRecoStartZ = kDefDoub;
  fSelTrackRecoStartT = kDefDoub;
  fSelTrackRecoEndX = kDefDoub;
  fSelTrackRecoEndY = kDefDoub;
  fSelTrackRecoEndZ = kDefDoub;
  fSelTrackRecoEndT = kDefDoub;
  fSelTrackRecoUpstreamX = kDefDoub;
  fSelTrackRecoUpstreamY = kDefDoub;
  fSelTrackRecoUpstreamZ = kDefDoub;
  fSelTrackRecoUpstreamT = kDefDoub;
  fSelTrackRecoDownstreamX = kDefDoub;
  fSelTrackRecoDownstreamY = kDefDoub;
  fSelTrackRecoDownstreamZ = kDefDoub;
  fSelTrackRecoDownstreamT = kDefDoub;
  fSelTrackRecoEndClosestToVertexX = kDefDoub;
  fSelTrackRecoEndClosestToVertexY = kDefDoub;
  fSelTrackRecoEndClosestToVertexZ = kDefDoub;
  fSelTrackRecoLength = kDefDoub;
  fSelTrackRecoContained = 0;
  fSelTrackRecoCharge = kDefDoub;
  fSelTrackRecoMomMCS = kDefDoub;
  fSelTrackRecoMomContained = kDefDoub;
  //fSelTrackRecoNMatchedVertices = kDefInt;
  fSelTrackRecoVertexX = kDefDoub;
  fSelTrackRecoVertexY = kDefDoub;
  fSelTrackRecoVertexZ = kDefDoub;
  fSelTrackRecoNChildPFP = kDefInt;
  fSelTrackRecoNChildTrackPFP = kDefInt;
  fSelTrackRecoNChildShowerPFP = kDefInt;
  //MVA bits
  fSelTrackMVAElectron = kDefDoub;
  fSelTrackMVAPion = kDefDoub;
  fSelTrackMVAMuon = kDefDoub;
  fSelTrackMVAProton = kDefDoub;
  fSelTrackMVAPhoton = kDefDoub;
  //Pandizzle
  fSelTrackPandizzleVar = kDefDoub;
  fTMVAPFPMichelNHits            = kDefDoub;
  fTMVAPFPMichelElectronMVA      = kDefDoub;
  fTMVAPFPMichelRecoEnergyPlane2 = kDefDoub;
  fTMVAPFPTrackDeflecAngleSD     = kDefDoub;
  fTMVAPFPTrackLength            = kDefDoub;
  fTMVAPFPTrackEvalRatio         = kDefDoub;
  fTMVAPFPTrackConcentration     = kDefDoub;
  fTMVAPFPTrackCoreHaloRatio     = kDefDoub;
  fTMVAPFPTrackConicalness       = kDefDoub;
  fTMVAPFPTrackdEdxStart         = kDefDoub;
  fTMVAPFPTrackdEdxEnd           = kDefDoub;
  fTMVAPFPTrackdEdxEndRatio      = kDefDoub;
  fTMVAPFPTrackPIDA              = kDefDoub;

  //Reco energy bits
  fNumuRecoEHad = kDefDoub;
  fNumuRecoMomLep = kDefDoub;
  fNumuRecoENu = kDefDoub; //Neutrino reco energy

  //Shower bits
  //true bits
  fSelShowerTruePDG = kDefInt;
  fSelShowerTruePrimary = kDefInt;
  fSelShowerTrueMomX = kDefDoub;
  fSelShowerTrueMomY = kDefDoub;
  fSelShowerTrueMomZ = kDefDoub;
  fSelShowerTrueMomT = kDefDoub;
  fSelShowerTrueStartX = kDefDoub;
  fSelShowerTrueStartY = kDefDoub;
  fSelShowerTrueStartZ = kDefDoub;
  fSelShowerTrueStartT = kDefDoub;
  fSelShowerTrueEndX = kDefDoub;
  fSelShowerTrueEndY = kDefDoub;
  fSelShowerTrueEndZ = kDefDoub;
  fSelShowerTrueEndT = kDefDoub;
  //reco bits
  fSelShowerRecoNHits = kDefInt;
  fSelShowerRecoCompleteness = kDefDoub;
  fSelShowerRecoHitPurity = kDefDoub;
  fSelShowerRecoMomX = kDefDoub;
  fSelShowerRecoMomY = kDefDoub;
  fSelShowerRecoMomZ = kDefDoub;
  fSelShowerRecoMomT = kDefDoub;
  fSelShowerRecoStartX = kDefDoub;
  fSelShowerRecoStartY = kDefDoub;
  fSelShowerRecoStartZ = kDefDoub;
  fSelShowerRecoStartT = kDefDoub;
  fSelShowerRecoEndX = kDefDoub;
  fSelShowerRecoEndY = kDefDoub;
  fSelShowerRecoEndZ = kDefDoub;
  fSelShowerRecoEndT = kDefDoub;
  fSelShowerRecoUpstreamX = kDefDoub;
  fSelShowerRecoUpstreamY = kDefDoub;
  fSelShowerRecoUpstreamZ = kDefDoub;
  fSelShowerRecoUpstreamT = kDefDoub;
  fSelShowerRecoDownstreamX = kDefDoub;
  fSelShowerRecoDownstreamY = kDefDoub;
  fSelShowerRecoDownstreamZ = kDefDoub;
  fSelShowerRecoDownstreamT = kDefDoub;
  fSelShowerRecoEndClosestToVertexX = kDefDoub;
  fSelShowerRecoEndClosestToVertexY = kDefDoub;
  fSelShowerRecoEndClosestToVertexZ = kDefDoub;
  fSelShowerRecoLength = kDefDoub;
  fSelShowerRecoContained = 0;
  fSelShowerRecoCharge = kDefDoub;
  fSelShowerRecoMomMCS = kDefDoub;
  fSelShowerRecoMomContained = kDefDoub;
  fSelShowerRecoNMatchedVertices = kDefInt;
  fSelShowerRecoVertexX = kDefDoub;
  fSelShowerRecoVertexY = kDefDoub;
  fSelShowerRecoVertexZ = kDefDoub;
  fSelShowerRecoNChildPFP = kDefInt;
  fSelShowerRecoNChildTrackPFP = kDefInt;
  fSelShowerRecoNChildShowerPFP = kDefInt;
  fSelShowerRecoDirX = kDefDoub;
  fSelShowerRecoDirY = kDefDoub;
  fSelShowerRecoDirZ = kDefDoub;
  fSelShowerRecodEdx = kDefDoub;
  fSelShowerRecoIsPrimaryPFPDaughter = 0;
  //MVA bits
  fSelShowerMVAElectron = kDefDoub;
  fSelShowerMVAPion = kDefDoub;
  fSelShowerMVAMuon = kDefDoub;
  fSelShowerMVAProton = kDefDoub;
  fSelShowerMVAPhoton = kDefDoub;
  //Reco nu bits
  fRecoNuVtxX = kDefDoub;
  fRecoNuVtxY = kDefDoub;
  fRecoNuVtxZ = kDefDoub;
  fRecoNuVtxNShowers = 0;
  fRecoNuVtxNTracks = 0;
  fRecoNuVtxNChildren = 0;



  //Reco energy bits
  fNueRecoEHad = kDefDoub;
  fNueRecoMomLep = kDefDoub;
  fNueRecoENu = kDefDoub; //Neutrino reco energy

  //Event level stuff
  fRecoEventCharge = kDefDoub;

}

void FDSelection::CCNuSelection::GetEventInfo(art::Event const & evt){
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fT0 = detprop->TriggerOffset();

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  //Loop over the hits to get the total charge
  fRecoEventCharge = 0.;
  for (unsigned int i_hit = 0; i_hit < hitList.size(); i_hit++){
    if (hitList[i_hit]->WireID().Plane == 2){
      fRecoEventCharge += hitList[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(hitList[i_hit]->PeakTime(), fT0);
    }
  }

  return;
}

void FDSelection::CCNuSelection::GetTruthInfo(art::Event const & evt){
  //Get the generator record
  art::Handle< std::vector<simb::MCTruth> > mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle)){
    art::fill_ptr_vector(mcList, mcTruthListHandle);
  }
  //Get the flux record
  art::Handle< std::vector<simb::MCFlux> > mcFluxListHandle;
  std::vector<art::Ptr<simb::MCFlux> > mcFlux;
  if (evt.getByLabel(fNuGenModuleLabel, mcFluxListHandle)){
    art::fill_ptr_vector(mcFlux, mcFluxListHandle);
  }

  //Chuck out a warning if there are multiple truths (i.e. multiple neutrinos)
  if (mcList.size() > 1){
    mf::LogWarning("CCNuSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one!!!!";
  }
  for (unsigned int i_mctruth = 0; i_mctruth < mcList.size(); i_mctruth++){
    fNuPdg    = mcList[i_mctruth]->GetNeutrino().Nu().PdgCode();
    fBeamPdg  = mcFlux[i_mctruth]->fntype;
    fNC       = mcList[i_mctruth]->GetNeutrino().CCNC();
    fMode     = mcList[i_mctruth]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    fTargetZ  = mcList[i_mctruth]->GetNeutrino().Target()%100000000/10000;
    fENu      = mcList[i_mctruth]->GetNeutrino().Nu().E();
    fQ2       = mcList[i_mctruth]->GetNeutrino().QSqr();
    fW        = mcList[i_mctruth]->GetNeutrino().W();
    fX        = mcList[i_mctruth]->GetNeutrino().X();
    fY        = mcList[i_mctruth]->GetNeutrino().Y();
    fNuMomX   = mcList[i_mctruth]->GetNeutrino().Nu().Momentum().X();
    fNuMomY   = mcList[i_mctruth]->GetNeutrino().Nu().Momentum().Y();
    fNuMomZ   = mcList[i_mctruth]->GetNeutrino().Nu().Momentum().Z();
    fNuMomT   = mcList[i_mctruth]->GetNeutrino().Nu().Momentum().T();
    //Lepton stuff
    fLepPDG     = mcList[i_mctruth]->GetNeutrino().Lepton().PdgCode();
    fMomLepX    = mcList[i_mctruth]->GetNeutrino().Lepton().Momentum().X();
    fMomLepY    = mcList[i_mctruth]->GetNeutrino().Lepton().Momentum().Y();
    fMomLepZ    = mcList[i_mctruth]->GetNeutrino().Lepton().Momentum().Z();
    fMomLepT    = mcList[i_mctruth]->GetNeutrino().Lepton().Momentum().T();
    fLepNuAngle = mcList[i_mctruth]->GetNeutrino().Nu().Momentum().Vect().Angle(mcList[i_mctruth]->GetNeutrino().Lepton().Momentum().Vect());
    fNuX = mcList[i_mctruth]->GetNeutrino().Nu().Vx();
    fNuY = mcList[i_mctruth]->GetNeutrino().Nu().Vy();
    fNuZ = mcList[i_mctruth]->GetNeutrino().Nu().Vz();
    fNuT = mcList[i_mctruth]->GetNeutrino().Nu().T();
    //Loop over the final state particles
    for (int i_part = 0; i_part < mcList[i_mctruth]->NParticles(); i_part++){
      const simb::MCParticle& final_state_particle = mcList[i_mctruth]->GetParticle(i_part);
      if (final_state_particle.StatusCode() != 1) continue;
      int pdg = final_state_particle.PdgCode();
      if (pdg >= 2000000000) continue;
      if (std::abs(pdg) >= 11 && std::abs(pdg) <= 16) continue;
      if (pdg==211) fNPiP++;
      else if (pdg==-211) fNPim++;
      else if (pdg==111) fNPi0++;
      else if (pdg==22) fNPhotons++;
      else if (pdg==2212) fNProtons++;
      else if (pdg==2112) fNNeutrons++;
      else fNOther++;
    }

    //20/11/18 DBrailsford
    //Get the associated g4 particles so that we can find the stop position of the lepton
    art::FindManyP<simb::MCParticle> fmpt(mcTruthListHandle, evt, fLargeantModuleLabel);
    const std::vector<art::Ptr<simb::MCParticle> > associated_particles = fmpt.at(mcList[i_mctruth].key());
    //Also get the MCParticle according to the MCTruth object so we can match
    for (unsigned int i_part = 0; i_part < associated_particles.size(); i_part++){
      art::Ptr<simb::MCParticle> particle = associated_particles[i_part];
      if (particle->PdgCode() == fLepPDG && particle->Mother()==0){
        //We should have found the primary lepton so let's get its end position
        fLepEndX = particle->EndPosition().X();
        fLepEndY = particle->EndPosition().Y();
        fLepEndZ = particle->EndPosition().Z();
        fLepEndT = particle->EndPosition().T();
      }
//      if (particle->PdgCode = fLepPDG)
    }

  }
}

void FDSelection::CCNuSelection::RunTrackSelection(art::Event const & evt){
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(evt.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return;
  }

  //Get the selected track
  art::Ptr<recob::Track> sel_track = fRecoTrackSelector->FindSelectedTrack(evt);

  //If we didn't find a selected track then what's the point?
  if (!(sel_track.isAvailable())) {
    std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no track returned from selection" << std::endl; 
    return;
  }

  //24/07/18 DBrailsford Get the reco energy data product for neutrinos
  art::Handle<dune::EnergyRecoOutput> energyRecoHandle;
  if (!evt.getByLabel(fNumuEnergyRecoModuleLabel, energyRecoHandle)) {
    std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no energy reconstruction found with label " << fNumuEnergyRecoModuleLabel << std::endl;
    return;
  }

  //Get the hits for said track
  art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > sel_track_hits = fmht.at(sel_track.key());
  fSelTrackRecoNHits = sel_track_hits.size();

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }


  //Start filling some variables
  recob::Track::Point_t trackStart, trackEnd;
  std::tie(trackStart, trackEnd) = sel_track->Extent(); 
  fSelTrackRecoMomX = kDefDoub; //temp
  fSelTrackRecoMomY = kDefDoub; //temp
  fSelTrackRecoMomZ = kDefDoub; //temp
  fSelTrackRecoMomT = kDefDoub; //temp
  fSelTrackRecoStartX = trackStart.X();
  fSelTrackRecoStartY = trackStart.Y();
  fSelTrackRecoStartZ = trackStart.Z();
  fSelTrackRecoEndX = trackEnd.X();
  fSelTrackRecoEndY = trackEnd.Y();
  fSelTrackRecoEndZ = trackEnd.Z();
  if (fSelTrackRecoEndZ > fSelTrackRecoStartZ){
    fSelTrackRecoUpstreamX = fSelTrackRecoStartX;
    fSelTrackRecoUpstreamY = fSelTrackRecoStartY;
    fSelTrackRecoUpstreamZ = fSelTrackRecoStartZ;
    fSelTrackRecoDownstreamX = fSelTrackRecoEndX;
    fSelTrackRecoDownstreamY = fSelTrackRecoEndY;
    fSelTrackRecoDownstreamZ = fSelTrackRecoEndZ;
  }
  else{
    fSelTrackRecoDownstreamX = fSelTrackRecoStartX;
    fSelTrackRecoDownstreamY = fSelTrackRecoStartY;
    fSelTrackRecoDownstreamZ = fSelTrackRecoStartZ;
    fSelTrackRecoUpstreamX = fSelTrackRecoEndX;
    fSelTrackRecoUpstreamY = fSelTrackRecoEndY;
    fSelTrackRecoUpstreamZ = fSelTrackRecoEndZ;
  }
  fSelTrackRecoLength = sel_track->Length();
  fSelTrackRecoCharge = CalculateTrackCharge(sel_track, sel_track_hits);
  //Now that the vertex information has been filled.  Calculate which end of th etrack is closest
  TVector3 upstream_end(fSelTrackRecoUpstreamX,fSelTrackRecoUpstreamY,fSelTrackRecoUpstreamZ);
  TVector3 downstream_end(fSelTrackRecoDownstreamX,fSelTrackRecoDownstreamY,fSelTrackRecoDownstreamZ);
  TVector3 vertex_pos(fSelTrackRecoVertexX,fSelTrackRecoVertexY,fSelTrackRecoVertexZ);
  if ((vertex_pos-upstream_end).Mag() < (vertex_pos-downstream_end).Mag()){
    fSelTrackRecoEndClosestToVertexX = fSelTrackRecoUpstreamX;
    fSelTrackRecoEndClosestToVertexY = fSelTrackRecoUpstreamY;
    fSelTrackRecoEndClosestToVertexZ = fSelTrackRecoUpstreamZ;
  }
  else{
    fSelTrackRecoEndClosestToVertexX = fSelTrackRecoDownstreamX;
    fSelTrackRecoEndClosestToVertexY = fSelTrackRecoDownstreamY;
    fSelTrackRecoEndClosestToVertexZ = fSelTrackRecoDownstreamZ;
  }
  //28/11/18 DBrailsford Fill child PFP info
  FillChildPFPInformation(sel_track, evt);
  //24/07/18 DBrailsford Use the data product to get the neutrino energy
  //fSelTrackRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);
  fNumuRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNumuRecoEHad = energyRecoHandle->fHadLorentzVector.E();

  fSelTrackRecoContained = energyRecoHandle->longestTrackContained; 
  if (energyRecoHandle->trackMomMethod==1){ //momentum by range was used to calculate ENu
    fSelTrackRecoMomContained = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fNumuRecoMomLep = fSelTrackRecoMomContained;
  }
  else if (energyRecoHandle->trackMomMethod==0){//momentum by MCS
    fSelTrackRecoMomMCS = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fNumuRecoMomLep = fSelTrackRecoMomMCS;
  }

  int g4id = FDSelectionUtils::TrueParticleIDFromTotalRecoHits(sel_track_hits);
  fSelTrackRecoCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(sel_track_hits, hitList, g4id);
  fSelTrackRecoHitPurity = FDSelectionUtils::HitPurityFromTrueParticleID(sel_track_hits, g4id);

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
  if (matched_mcparticle){
    //Fill variables
    fSelTrackTruePDG = matched_mcparticle->PdgCode();
    if (matched_mcparticle->Mother()==0) fSelTrackTruePrimary = 1;
    else fSelTrackTruePrimary = 0;
    fSelTrackTrueMomX = matched_mcparticle->Momentum().X();
    fSelTrackTrueMomY = matched_mcparticle->Momentum().Y();
    fSelTrackTrueMomZ = matched_mcparticle->Momentum().Z();
    fSelTrackTrueMomT = matched_mcparticle->Momentum().T();
    fSelTrackTrueStartX = matched_mcparticle->Position(0).X();
    fSelTrackTrueStartY = matched_mcparticle->Position(0).Y();
    fSelTrackTrueStartZ = matched_mcparticle->Position(0).Z();
    fSelTrackTrueStartT = matched_mcparticle->Position(0).T();
    fSelTrackTrueEndX = matched_mcparticle->EndPosition().X();
    fSelTrackTrueEndY = matched_mcparticle->EndPosition().Y();
    fSelTrackTrueEndZ = matched_mcparticle->EndPosition().Z();
    fSelTrackTrueEndT = matched_mcparticle->EndPosition().T();
  }
  //Now get the pid stuff
  art::FindManyP<anab::MVAPIDResult> fmpidt(trackListHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(sel_track.key());
  std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
  //Get the PIDs
  fSelTrackMVAElectron = mvaOutMap["electron"];
  fSelTrackMVAPion = mvaOutMap["pich"];
  fSelTrackMVAMuon = mvaOutMap["muon"];
  fSelTrackMVAProton = mvaOutMap["proton"];
  fSelTrackMVAPhoton = mvaOutMap["photon"];

  //11/04/19 DBrailsford
  //Get the pandizzle variables
  art::Ptr<recob::PFParticle> track_pfp = GetPFParticleMatchedToTrack(sel_track, evt);
  fPandizzleAlg.ProcessPFParticle(track_pfp, evt);
  fTMVAPFPMichelNHits = (float)(fPandizzleAlg.GetIntVar("PFPMichelNHits"));
  fTMVAPFPMichelElectronMVA = fPandizzleAlg.GetFloatVar("PFPMichelElectronMVA");
  fTMVAPFPMichelRecoEnergyPlane2 = fPandizzleAlg.GetFloatVar("PFPMichelRecoEnergyPlane2");
  fTMVAPFPTrackDeflecAngleSD = fPandizzleAlg.GetFloatVar("PFPTrackDeflecAngleSD");
  fTMVAPFPTrackLength = fPandizzleAlg.GetFloatVar("PFPTrackLength");
  fTMVAPFPTrackEvalRatio = fPandizzleAlg.GetFloatVar("PFPTrackEvalRatio");
  fTMVAPFPTrackConcentration = fPandizzleAlg.GetFloatVar("PFPTrackConcentration");
  fTMVAPFPTrackCoreHaloRatio = fPandizzleAlg.GetFloatVar("PFPTrackCoreHaloRatio");
  fTMVAPFPTrackConicalness = fPandizzleAlg.GetFloatVar("PFPTrackConicalness");
  fTMVAPFPTrackdEdxStart = fPandizzleAlg.GetFloatVar("PFPTrackdEdxStart");
  fTMVAPFPTrackdEdxEnd = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEnd");
  fTMVAPFPTrackdEdxEndRatio = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEndRatio");
  fTMVAPFPTrackPIDA = fPandizzleAlg.GetFloatVar("PFPTrackPIDA");

  fSelTrackPandizzleVar = fReader.EvaluateMVA("BDTG");
}

double FDSelection::CCNuSelection::CalculateTrackCharge(art::Ptr<recob::Track> const track, std::vector< art::Ptr< recob::Hit> > const track_hits){
  double charge = 0;
  for (unsigned int i_hit = 0; i_hit < track_hits.size(); i_hit++){
    if (track_hits[i_hit]->WireID().Plane != 2) continue; //Collection only
    charge += track_hits[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(track_hits[i_hit]->PeakTime(), fT0);
  }
  return charge;
}


bool FDSelection::CCNuSelection::IsTrackContained(art::Ptr<recob::Track> const track, std::vector< art::Ptr<recob::Hit > > const track_hits, art::Event const & evt){
  //Get the space points for each of the hits
  //Annoyingly we have to go from the start of the handle for the hits...
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }
  art::FindManyP<recob::SpacePoint> fmhs(hitListHandle, evt, fTrackModuleLabel);
  for (unsigned int i_hit = 0; i_hit < track_hits.size(); i_hit++){
    if (track_hits[i_hit]->WireID().Plane != 2) continue;
    std::vector<art::Ptr<recob::SpacePoint> > space_points = fmhs.at(track_hits[i_hit].key());
    if (space_points.size()){
      //Make a TVector3
      TVector3 position(space_points[0]->XYZ()[0],space_points[0]->XYZ()[1],space_points[0]->XYZ()[2]);
      bool is_in_tpc = FDSelectionUtils::IsInsideTPC(position,20); //20cm buffer from the wall
      if (!is_in_tpc){
        return false;
      }
    }
  }
  return true;
}

void FDSelection::CCNuSelection::FillVertexInformation(art::Event const & evt){
  //Have to get the PFP and then get the VTX
  //Need the PFP vector handle
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);



  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);
  //Loop over the neutrinos
  if (nu_pfps.size() != 1){
    std::cout<<"FDSelection::CCNuSelection: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return; //do nothing
  }
  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  //Count the number of children
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);
  for (int i_childpfp = 0; i_childpfp < nu_pfp->NumDaughters(); i_childpfp++){
    int child_id = nu_pfp->Daughter(i_childpfp);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[child_id];
    int pdg = child_pfp->PdgCode();
    if (pdg == 11) fRecoNuVtxNShowers++;
    else if (pdg == 13) fRecoNuVtxNTracks++;
    fRecoNuVtxNChildren++;
  }

  art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());
  /*
  fSelTrackRecoNMatchedVertices = sel_pfp_vertices.size();
  if (fSelTrackRecoNMatchedVertices == 0){ //Nothing to do
    return;
  }
  else if (fSelTrackRecoNMatchedVertices > 1){
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << fSelTrackRecoNMatchedVertices << std::endl;
  }
  */
  if (sel_pfp_vertices.size() == 0){ //Nothing to do
    return;
  }
  else if (sel_pfp_vertices.size() > 1){
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  return;
}

void FDSelection::CCNuSelection::FillChildPFPInformation(art::Ptr<recob::Track> const track, art::Event const & evt){
  //Get the PFP from the track
  art::Ptr<recob::PFParticle> matched_pfp = GetPFParticleMatchedToTrack(track, evt);
  if (!(matched_pfp.isAvailable())) return;

  //Now make the entire PFPPArticleMap
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  fSelTrackRecoNChildPFP = matched_pfp->NumDaughters();
  fSelTrackRecoNChildTrackPFP = 0;
  fSelTrackRecoNChildShowerPFP = 0;
  for (int i_child = 0; i_child < matched_pfp->NumDaughters(); i_child++){
    int child_id = matched_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[child_id];
    int pdg = child_pfp->PdgCode();
    if (pdg==13) fSelTrackRecoNChildTrackPFP++;
    else if (pdg==11) fSelTrackRecoNChildShowerPFP++;
    else std::cout<<"FillChildPFPInformation: found a child PFP with an unexpected pdg code: " << pdg << std::endl;
  }

  return;
}


art::Ptr<recob::PFParticle> FDSelection::CCNuSelection::GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt){
  art::Ptr<recob::PFParticle> matched_pfp;
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(evt.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToTrack Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return matched_pfp;
  }
  art::FindManyP<recob::PFParticle> fmpfpt(trackListHandle, evt, fTrackModuleLabel);
  const std::vector<art::Ptr<recob::PFParticle> > sel_track_pfps = fmpfpt.at(track.key());
  if (sel_track_pfps.size() != 1){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToTrack NUMBER OF PFP MATCHED TO A TRACK DOES NOT EQUAL 1: " << sel_track_pfps.size() << std::endl;
    return matched_pfp;
  }
  matched_pfp = sel_track_pfps[0];
  return matched_pfp;
}

art::Ptr<recob::Track> FDSelection::CCNuSelection::GetTrackMatchedPFParticle(art::Ptr<recob::PFParticle> const pfparticle, art::Event const & evt){
  art::Ptr<recob::Track> matched_track;
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return matched_track; //empty
  }
  art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Track> > sel_pfp_tracks = fmtpfp.at(pfparticle.key());
  if (sel_pfp_tracks.size() != 1){
    std::cout<<"CCNuSelection::GetTrackMatchedPFParticle number of tracks matched to PFP does not equal 1: " << sel_pfp_tracks.size() << std::endl;
    return matched_track; //empty
  }
  matched_track = sel_pfp_tracks[0];
  return matched_track;
}






void FDSelection::CCNuSelection::RunShowerSelection(art::Event const & evt){
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle))){
    std::cout<<"Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
    return;
  }
  //std::vector<art::Ptr<recob::Shower> > showers;
  //art::fill_ptr_vector(showers,showerListHandle);


  //Get the selected shower
  art::Ptr<recob::Shower> sel_shower = fRecoShowerSelector->FindSelectedShower(evt);

  //If we didn't find a selected track then what's the point?
  if (!(sel_shower.isAvailable())) {
    std::cout<<"FDSelection::CCNuSelection::RunShowerSelection - no shower selected by tool"<<std::endl;
    return;
  }

  //10/10/18 DBrailsford start assessing PFParticle stuff
  //Now get the PFParticles
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle > > pfparticleList;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  else art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  //Get the PFParticle associated with the selected shower
  art::FindOneP<recob::PFParticle> fospfp(showerListHandle, evt, fShowerModuleLabel);
  const art::Ptr<recob::PFParticle> sel_pfp = fospfp.at(sel_shower.key());
  if (!(sel_pfp.isAvailable())){
    std::cout<<"Did not find a matched PFPArticle for this shower.  I don't think this should ever happen"<<std::endl;
    return;
  }
  //Get the neutrino parent PFParticle.  To do this we have to build the LArPandoraHelper PFParticle map
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);
  art::Ptr<recob::PFParticle> nu_pfp = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfparticleMap, sel_pfp);
  if (!(nu_pfp.isAvailable())){
    std::cout<<"Was not able to find the primary neutrino PFP after selecting a shower.  I don't think this should ever happen"<<std::endl;
  }
  /*
  //Now get the associated vertex
  art::FindOneP<recob::Vertex> fopfpv(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const art::Ptr<recob::Vertex> nu_vertex = fopfpv.at(nu_pfp.key());
  if (!(nu_vertex.isAvailable())){
    std::cout<<"Was not able to find the reco vertex after finding the neutrino PFP.  I don't think this should ever happen"<<std::endl;
  }
  */



  //24/07/18 DBrailsford Get the reco energy data product for neutrinos
  art::Handle<dune::EnergyRecoOutput> energyRecoHandle;
  if (!evt.getByLabel(fNueEnergyRecoModuleLabel, energyRecoHandle)) {
    std::cout<<"FDSelection::CCNuSelection::RunShowerSelection - Not able to find energy reconstruction container with name " << fNueEnergyRecoModuleLabel << std::endl;
    return;
  }

  //Get the hits for said shower
  art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > sel_shower_hits = fmhs.at(sel_shower.key());
  fSelShowerRecoNHits = sel_shower_hits.size();

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  ////Start filling
  ////Get the reconstructed neutrino vertex position
  //fRecoNuVtxX = nu_vertex->position().X();
  //fRecoNuVtxY = nu_vertex->position().Y();
  //fRecoNuVtxZ = nu_vertex->position().Z();

  //Is this PFParticle a primary daughter of the neutrino
  fSelShowerRecoIsPrimaryPFPDaughter = IsPFParticlePrimaryDaughter(sel_pfp, pfparticleList);



  fSelShowerRecoDirX   = sel_shower->Direction().X();
  fSelShowerRecoDirY   = sel_shower->Direction().Y();
  fSelShowerRecoDirZ   = sel_shower->Direction().Z();
  fSelShowerRecoStartX = sel_shower->ShowerStart().X();
  fSelShowerRecoStartY = sel_shower->ShowerStart().Y();
  fSelShowerRecoStartZ = sel_shower->ShowerStart().Z();

  fSelShowerRecoCharge = CalculateShowerCharge(sel_shower, sel_shower_hits);
  //24/07/18 DBrailsford Use the data product to get the neutrino energy
  //fSelShowerRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);
  fNueRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNueRecoEHad = energyRecoHandle->fHadLorentzVector.E();
  fNueRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

  int g4id = FDSelectionUtils::TrueParticleIDFromTotalRecoHits(sel_shower_hits);
  fSelShowerRecoCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(sel_shower_hits, hitList, g4id);
  fSelShowerRecoHitPurity = FDSelectionUtils::HitPurityFromTrueParticleID(sel_shower_hits, g4id);

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  //const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
  const simb::MCParticle* matched_mcparticle = pi_serv->TrackIdToParticle_P(g4id);
  if (matched_mcparticle){
    //Fill variables
    fSelShowerTruePDG = matched_mcparticle->PdgCode();
    if (matched_mcparticle->Mother()==0) fSelShowerTruePrimary = 1;
    else fSelShowerTruePrimary = 0;
    fSelShowerTrueMomX = matched_mcparticle->Momentum().X();
    fSelShowerTrueMomY = matched_mcparticle->Momentum().Y();
    fSelShowerTrueMomZ = matched_mcparticle->Momentum().Z();
    fSelShowerTrueMomT = matched_mcparticle->Momentum().T();
    fSelShowerTrueStartX = matched_mcparticle->Position(0).X();
    fSelShowerTrueStartY = matched_mcparticle->Position(0).Y();
    fSelShowerTrueStartZ = matched_mcparticle->Position(0).Z();
    fSelShowerTrueStartT = matched_mcparticle->Position(0).T();
    fSelShowerTrueEndX = matched_mcparticle->EndPosition().X();
    fSelShowerTrueEndY = matched_mcparticle->EndPosition().Y();
    fSelShowerTrueEndZ = matched_mcparticle->EndPosition().Z();
    fSelShowerTrueEndT = matched_mcparticle->EndPosition().T();
  }
  //Now get the pid stuff
  art::FindManyP<anab::MVAPIDResult> fmpidt(showerListHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(sel_shower.key());
  std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
  //Get the PIDs
  fSelShowerMVAElectron = mvaOutMap["electron"];
  fSelShowerMVAPion = mvaOutMap["pich"];
  fSelShowerMVAMuon = mvaOutMap["muon"];
  fSelShowerMVAProton = mvaOutMap["proton"];
  fSelShowerMVAPhoton = mvaOutMap["photon"];

}

double FDSelection::CCNuSelection::CalculateShowerCharge(art::Ptr<recob::Shower> const shower, std::vector< art::Ptr< recob::Hit> > const shower_hits){
  double charge = 0;
  for (unsigned int i_hit = 0; i_hit < shower_hits.size(); i_hit++){
    if (shower_hits[i_hit]->WireID().Plane != 2) continue; //Collection only
    charge += shower_hits[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(shower_hits[i_hit]->PeakTime(), fT0);
  }
  return charge;
}

bool FDSelection::CCNuSelection::IsPFParticlePrimaryDaughter(art::Ptr<recob::PFParticle> const & pfparticle, std::vector<art::Ptr<recob::PFParticle > > const & pfparticles){
  size_t parent_id = pfparticle->Parent();
  for (unsigned int i_pfp = 0; i_pfp < pfparticles.size(); i_pfp++){
    art::Ptr<recob::PFParticle> parent_candidate = pfparticles[i_pfp];
    if (parent_candidate->Self() == parent_id && parent_candidate->IsPrimary()){
      return true;
    }
  }
  return false;
}







DEFINE_ART_MODULE(FDSelection::CCNuSelection)

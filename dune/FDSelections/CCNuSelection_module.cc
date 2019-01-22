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


//DUNE
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

//Custom
//#include "PIDAnaAlg.h"
#include "PandizzleAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"


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
  void RunSelection(art::Event const & evt);  //Run the selection and dump relevant info to the truee
  double CalculateTrackCharge(art::Ptr<recob::Track> const track, std::vector< art::Ptr< recob::Hit> > const track_hits); //Calculate hit charge on track as measured by the collection plane
  bool IsTrackContained(art::Ptr<recob::Track> const track, std::vector< art::Ptr<recob::Hit > > const track_hits, art::Event const & evt); // check if the track is contained in the detector
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt);
  art::Ptr<recob::Track> GetTrackMatchedPFParticle(art::Ptr<recob::PFParticle> const pfparticle, art::Event const & evt);

  void FillVertexInformation(art::Ptr<recob::Track> const track, art::Event const & evt); //Grab the vertex parameters associated with this track
  void FillChildPFPInformation(art::Ptr<recob::Track> const track, art::Event const & evt); //Grab the vertex parameters associated with this track


  // Declare member data here.

  //Algs
  //PIDAnaAlg fPIDAnaAlg;
  PandizzleAlg fPandizzleAlg;
  calo::CalorimetryAlg fCalorimetryAlg;

  //Tools
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;

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
  int fSelTrackRecoNMatchedVertices;
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
  //Event-level bits
  double fRecoEventCharge; //Total collected charge (as measured by the collection planes)
  //reco energy bits
  double fNumuRecoMomLep; //Reco lepton momentum
  double fNumuRecoEHad; //Reco hadronic energy
  double fNumuRecoENu; //Reco neutrino energy

  //POT tree stuff
  TTree* fPOTTree;
  double fPOT;

  //Fhicl pset labels
  std::string fNuGenModuleLabel;
  std::string fLargeantModuleLabel;
  std::string fTrackModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fVertexModuleLabel;
  std::string fPIDModuleLabel;
  std::string fHitsModuleLabel;
  std::string fPOTModuleLabel;
  std::string fEnergyRecoModuleLabel;
 
  //Processing flags
  bool fUsePandoraVertex;

};


FDSelection::CCNuSelection::CCNuSelection(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)   ,
  //fPIDAnaAlg(pset.get<fhicl::ParameterSet>("ModuleLabels"))   ,
  fPandizzleAlg(pset.get<fhicl::ParameterSet>("ModuleLabels")) ,
  fCalorimetryAlg          (pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fNuGenModuleLabel        (pset.get< std::string >("ModuleLabels.NuGenModuleLabel")),
  fLargeantModuleLabel     (pset.get< std::string >("ModuleLabels.LargeantModuleLabel")),
  fTrackModuleLabel        (pset.get< std::string >("ModuleLabels.TrackModuleLabel")),
  fPFParticleModuleLabel   (pset.get< std::string >("ModuleLabels.PFParticleModuleLabel")),
  fVertexModuleLabel   (pset.get< std::string >("ModuleLabels.VertexModuleLabel")),
  fPIDModuleLabel          (pset.get< std::string >("ModuleLabels.PIDModuleLabel")),
  fHitsModuleLabel         (pset.get< std::string >("ModuleLabels.HitsModuleLabel")),
  fPOTModuleLabel          (pset.get< std::string >("ModuleLabels.POTModuleLabel")),
  fEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.EnergyRecoModuleLabel")),
  fUsePandoraVertex        (pset.get< bool >("UsePandoraVertex")) {}

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
  //Get the generic stuff that can be pulled from the top of the record
  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();
  fIsMC = !(evt.isRealData());

  GetEventInfo(evt);
  if (fIsMC) GetTruthInfo(evt);
  RunSelection(evt);

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
    fTree->Branch("SelTrackRecoNMatchedVertices",&fSelTrackRecoNMatchedVertices);
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
    fTree->Branch("RecoEventCharge",&fRecoEventCharge);
    fTree->Branch("NumuRecoMomLep",&fNumuRecoMomLep);
    fTree->Branch("NumuRecoEHad",&fNumuRecoEHad);
    fTree->Branch("NumuRecoENu",&fNumuRecoENu);


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
  fSelTrackRecoNMatchedVertices = kDefInt;
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
  //Event level stuff
  fRecoEventCharge = kDefDoub;
  //Reco energy bits
  fNumuRecoEHad = kDefDoub;
  fNumuRecoMomLep = kDefDoub;
  fNumuRecoENu = kDefDoub; //Neutrino reco energy
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

void FDSelection::CCNuSelection::RunSelection(art::Event const & evt){
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(evt.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return;
  }

  //Get the selected track
  art::Ptr<recob::Track> sel_track = fRecoTrackSelector->FindSelectedTrack(evt);

  //If we didn't find a selected track then what's the point?
  if (!(sel_track.isAvailable())) return;

  //24/07/18 DBrailsford Get the reco energy data product for neutrinos
  art::Handle<dune::EnergyRecoOutput> energyRecoHandle;
  if (!evt.getByLabel(fEnergyRecoModuleLabel, energyRecoHandle)) return;

  //Get the hits for said track
  art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > sel_track_hits = fmht.at(sel_track.key());

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
  //27/11/18 DBrailsford Fill vertex information
  FillVertexInformation(sel_track, evt);
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

void FDSelection::CCNuSelection::FillVertexInformation(art::Ptr<recob::Track> const track, art::Event const & evt){
  if (fUsePandoraVertex){
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
    art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
    const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());
    fSelTrackRecoNMatchedVertices = sel_pfp_vertices.size();
    if (fSelTrackRecoNMatchedVertices == 0){ //Nothing to do
      return;
    }
    else if (fSelTrackRecoNMatchedVertices > 1){
      std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << fSelTrackRecoNMatchedVertices << std::endl;
    }
    //always take the first vertex, even if there's more than one
    art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
    fSelTrackRecoVertexX = matched_vertex->position().X();
    fSelTrackRecoVertexY = matched_vertex->position().Y();
    fSelTrackRecoVertexZ = matched_vertex->position().Z();

    /*
    art::Ptr<recob::PFParticle> matched_pfp =GetPFParticleMatchedToTrack(track, evt);

    art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
    const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(matched_pfp.key());
    fSelTrackRecoNMatchedVertices = sel_pfp_vertices.size();
    if (fSelTrackRecoNMatchedVertices == 0){ //Nothing to do
      return;
    }
    else if (fSelTrackRecoNMatchedVertices > 1){
      std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << fSelTrackRecoNMatchedVertices << std::endl;
    }
    //always take the first vertex, even if there's more than one
    art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
    fSelTrackRecoVertexX = matched_vertex->position().X();
    fSelTrackRecoVertexY = matched_vertex->position().Y();
    fSelTrackRecoVertexZ = matched_vertex->position().Z();
    //Count how many PFPs are matched to this vertex
    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    if (!(evt.getByLabel(fVertexModuleLabel, vertexListHandle))){
      std::cout<<"CCNuSelection::FillVertexInformation Unable to find std::vector<recob::Vertex> with module label: " << fVertexModuleLabel << std::endl;
      return;
    }

    std::vector<art::Ptr<recob::Vertex> > vertexList;
    art::fill_ptr_vector(vertexList, vertexListHandle);

    art::FindManyP<recob::PFParticle> fmpfpv(vertexListHandle, evt, fVertexModuleLabel);
    const std::vector<art::Ptr<recob::PFParticle> > sel_vertex_pfps = fmpfpv.at(matched_vertex.key());
    */

  }
  else{
    fSelTrackRecoVertexX = track->Start().X();
    fSelTrackRecoVertexY = track->Start().Y();
    fSelTrackRecoVertexZ = track->Start().Z();
  }
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


DEFINE_ART_MODULE(FDSelection::CCNuSelection)

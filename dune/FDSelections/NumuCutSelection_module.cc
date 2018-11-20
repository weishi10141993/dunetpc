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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
//DUNE
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

//Custom
#include "PIDAnaAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"


constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

namespace FDSelection {
  class NumuCutSelection;
}


class FDSelection::NumuCutSelection : public art::EDAnalyzer {
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

  // Declare member data here.

  //Algs
  PIDAnaAlg fPIDAnaAlg;
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
  int fSelTruePDG;
  int fSelTruePrimary;
  double fSelTrueMomX;
  double fSelTrueMomY;
  double fSelTrueMomZ;
  double fSelTrueMomT;
  double fSelTrueStartX;
  double fSelTrueStartY;
  double fSelTrueStartZ;
  double fSelTrueStartT;
  double fSelTrueEndX;
  double fSelTrueEndY;
  double fSelTrueEndZ;
  double fSelTrueEndT;
  //reco bits
  double fSelRecoMomX;
  double fSelRecoMomY;
  double fSelRecoMomZ;
  double fSelRecoMomT;
  double fSelRecoStartX;
  double fSelRecoStartY;
  double fSelRecoStartZ;
  double fSelRecoStartT;
  double fSelRecoEndX;
  double fSelRecoEndY;
  double fSelRecoEndZ;
  double fSelRecoEndT;
  double fSelRecoUpstreamX;
  double fSelRecoUpstreamY;
  double fSelRecoUpstreamZ;
  double fSelRecoUpstreamT;
  double fSelRecoDownstreamX;
  double fSelRecoDownstreamY;
  double fSelRecoDownstreamZ;
  double fSelRecoDownstreamT;
  double fSelRecoLength;
  bool   fSelRecoContained;
  double fSelRecoCharge;
  double fSelRecoMomMCS;
  double fSelRecoMomContained;
  //MVA bits
  double fSelMVAElectron;
  double fSelMVAPion;
  double fSelMVAMuon;
  double fSelMVAProton;
  double fSelMVAPhoton;
  //Event-level bits
  double fRecoEventCharge; //Total collected charge (as measured by the collection planes)
  //reco energy bits
  double fRecoMomLep; //Reco lepton momentum
  double fRecoEHad; //Reco hadronic energy
  double fRecoENu; //Reco neutrino energy

  //POT tree stuff
  TTree* fPOTTree;
  double fPOT;

  //Fhicl pset labels
  std::string fNuGenModuleLabel;
  std::string fLargeantModuleLabel;
  std::string fTrackModuleLabel;
  std::string fPIDModuleLabel;
  std::string fHitsModuleLabel;
  std::string fPOTModuleLabel;
  std::string fEnergyRecoModuleLabel;
 

};


FDSelection::NumuCutSelection::NumuCutSelection(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)   ,
  fPIDAnaAlg(pset.get<fhicl::ParameterSet>("ModuleLabels"))   ,
  fCalorimetryAlg          (pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fNuGenModuleLabel        (pset.get< std::string >("ModuleLabels.NuGenModuleLabel")),
  fLargeantModuleLabel     (pset.get< std::string >("ModuleLabels.LargeantModuleLabel")),
  fTrackModuleLabel        (pset.get< std::string >("ModuleLabels.TrackModuleLabel")),
  fPIDModuleLabel          (pset.get< std::string >("ModuleLabels.PIDModuleLabel")),
  fHitsModuleLabel         (pset.get< std::string >("ModuleLabels.HitsModuleLabel")),
  fPOTModuleLabel          (pset.get< std::string >("ModuleLabels.POTModuleLabel")),
  fEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.EnergyRecoModuleLabel"))
{}

void FDSelection::NumuCutSelection::analyze(art::Event const & evt)
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

  fTree->Fill();
  Reset(); //Reset at the end of the event
}

void FDSelection::NumuCutSelection::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("numucutsel","Numu cut selection");
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
    fTree->Branch("SelTruePDG",&fSelTruePDG);
    fTree->Branch("SelTruePrimary",&fSelTruePrimary);
    fTree->Branch("SelTrueMomX",&fSelTrueMomX);
    fTree->Branch("SelTrueMomY",&fSelTrueMomY);
    fTree->Branch("SelTrueMomZ",&fSelTrueMomZ);
    fTree->Branch("SelTrueMomT",&fSelTrueMomT);
    fTree->Branch("SelTrueStartX",&fSelTrueStartX);
    fTree->Branch("SelTrueStartY",&fSelTrueStartY);
    fTree->Branch("SelTrueStartZ",&fSelTrueStartZ);
    fTree->Branch("SelTrueStartT",&fSelTrueStartT);
    fTree->Branch("SelTrueEndX",&fSelTrueEndX);
    fTree->Branch("SelTrueEndY",&fSelTrueEndY);
    fTree->Branch("SelTrueEndZ",&fSelTrueEndZ);
    fTree->Branch("SelTrueEndT",&fSelTrueEndT);
    fTree->Branch("SelRecoMomX",&fSelRecoMomX);
    fTree->Branch("SelRecoMomY",&fSelRecoMomY);
    fTree->Branch("SelRecoMomZ",&fSelRecoMomZ);
    fTree->Branch("SelRecoMomT",&fSelRecoMomT);
    fTree->Branch("SelRecoStartX",&fSelRecoStartX);
    fTree->Branch("SelRecoStartY",&fSelRecoStartY);
    fTree->Branch("SelRecoStartZ",&fSelRecoStartZ);
    fTree->Branch("SelRecoStartT",&fSelRecoStartT);
    fTree->Branch("SelRecoEndX",&fSelRecoEndX);
    fTree->Branch("SelRecoEndY",&fSelRecoEndY);
    fTree->Branch("SelRecoEndZ",&fSelRecoEndZ);
    fTree->Branch("SelRecoEndT",&fSelRecoEndT);
    fTree->Branch("SelRecoUpstreamX",&fSelRecoUpstreamX);
    fTree->Branch("SelRecoUpstreamY",&fSelRecoUpstreamY);
    fTree->Branch("SelRecoUpstreamZ",&fSelRecoUpstreamZ);
    fTree->Branch("SelRecoUpstreamT",&fSelRecoUpstreamT);
    fTree->Branch("SelRecoDownstreamX",&fSelRecoDownstreamX);
    fTree->Branch("SelRecoDownstreamY",&fSelRecoDownstreamY);
    fTree->Branch("SelRecoDownstreamZ",&fSelRecoDownstreamZ);
    fTree->Branch("SelRecoDownstreamT",&fSelRecoDownstreamT);
    fTree->Branch("SelRecoLength",&fSelRecoLength);
    fTree->Branch("SelRecoContained",&fSelRecoContained);
    fTree->Branch("SelRecoCharge",&fSelRecoCharge);
    fTree->Branch("SelRecoMomMCS",&fSelRecoMomMCS);
    fTree->Branch("SelRecoMomContained",&fSelRecoMomContained);
    fTree->Branch("SelMVAElectron",&fSelMVAElectron);
    fTree->Branch("SelMVAPion",&fSelMVAPion);
    fTree->Branch("SelMVAMuon",&fSelMVAMuon);
    fTree->Branch("SelMVAProton",&fSelMVAProton);
    fTree->Branch("SelMVAPhoton",&fSelMVAPhoton);
    fTree->Branch("RecoEventCharge",&fRecoEventCharge);
    fTree->Branch("RecoMomLep",&fRecoMomLep);
    fTree->Branch("RecoEHad",&fRecoEHad);
    fTree->Branch("RecoENu",&fRecoENu);


    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("POT",&fPOT);
    fPOTTree->Branch("Run",&fRun);
    fPOTTree->Branch("SubRun",&fSubRun);


    Reset();  //Default value all variables now
}

void FDSelection::NumuCutSelection::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}
void FDSelection::NumuCutSelection::endSubRun(const art::SubRun& sr){
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

void FDSelection::NumuCutSelection::endJob()
{
  // Implementation of optional member function here.
}

void FDSelection::NumuCutSelection::Reset()
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
  fSelTruePDG = kDefInt;
  fSelTruePrimary = kDefInt;
  fSelTrueMomX = kDefDoub;
  fSelTrueMomY = kDefDoub;
  fSelTrueMomZ = kDefDoub;
  fSelTrueMomT = kDefDoub;
  fSelTrueStartX = kDefDoub;
  fSelTrueStartY = kDefDoub;
  fSelTrueStartZ = kDefDoub;
  fSelTrueStartT = kDefDoub;
  fSelTrueEndX = kDefDoub;
  fSelTrueEndY = kDefDoub;
  fSelTrueEndZ = kDefDoub;
  fSelTrueEndT = kDefDoub;
  //reco bits
  fSelRecoMomX = kDefDoub;
  fSelRecoMomY = kDefDoub;
  fSelRecoMomZ = kDefDoub;
  fSelRecoMomT = kDefDoub;
  fSelRecoStartX = kDefDoub;
  fSelRecoStartY = kDefDoub;
  fSelRecoStartZ = kDefDoub;
  fSelRecoStartT = kDefDoub;
  fSelRecoEndX = kDefDoub;
  fSelRecoEndY = kDefDoub;
  fSelRecoEndZ = kDefDoub;
  fSelRecoEndT = kDefDoub;
  fSelRecoUpstreamX = kDefDoub;
  fSelRecoUpstreamY = kDefDoub;
  fSelRecoUpstreamZ = kDefDoub;
  fSelRecoUpstreamT = kDefDoub;
  fSelRecoDownstreamX = kDefDoub;
  fSelRecoDownstreamY = kDefDoub;
  fSelRecoDownstreamZ = kDefDoub;
  fSelRecoDownstreamT = kDefDoub;
  fSelRecoLength = kDefDoub;
  fSelRecoContained = 0;
  fSelRecoCharge = kDefDoub;
  fSelRecoMomMCS = kDefDoub;
  fSelRecoMomContained = kDefDoub;

  //MVA bits
  fSelMVAElectron = kDefDoub;
  fSelMVAPion = kDefDoub;
  fSelMVAMuon = kDefDoub;
  fSelMVAProton = kDefDoub;
  fSelMVAPhoton = kDefDoub;
  //Event level stuff
  fRecoEventCharge = kDefDoub;
  //Reco energy bits
  fRecoEHad = kDefDoub;
  fRecoMomLep = kDefDoub;
  fRecoENu = kDefDoub; //Neutrino reco energy
}

void FDSelection::NumuCutSelection::GetEventInfo(art::Event const & evt){
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

void FDSelection::NumuCutSelection::GetTruthInfo(art::Event const & evt){
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
    mf::LogWarning("NumuCutSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one!!!!";
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

void FDSelection::NumuCutSelection::RunSelection(art::Event const & evt){
  //art::Ptr<recob::Track> myseltrack = fRecoTrackSelector->FindSelectedTrack(evt);
  //if (myseltrack.isAvailable()) std::cout<<"Alright found a track"<<std::endl;
  //else std::cout<<"oh :("<<std::endl;

  /*
  int i_longest_track = -1;
  fSelRecoLength = -999;
  //Loop over the tracks to get the longest one NICE
  for (unsigned int i_track = 0; i_track < trackList.size(); i_track++){
    double current_track_length = trackList[i_track]->Length();
    if (current_track_length > fSelRecoLength){
      fSelRecoLength = current_track_length;
      i_longest_track = i_track;
    }
  }
  */
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
  fSelRecoMomX = kDefDoub; //temp
  fSelRecoMomY = kDefDoub; //temp
  fSelRecoMomZ = kDefDoub; //temp
  fSelRecoMomT = kDefDoub; //temp
  fSelRecoStartX = trackStart.X();
  fSelRecoStartY = trackStart.Y();
  fSelRecoStartZ = trackStart.Z();
  fSelRecoEndX = trackEnd.X();
  fSelRecoEndY = trackEnd.Y();
  fSelRecoEndZ = trackEnd.Z();
  if (fSelRecoEndZ > fSelRecoStartZ){
    fSelRecoUpstreamX = fSelRecoStartX;
    fSelRecoUpstreamY = fSelRecoStartY;
    fSelRecoUpstreamZ = fSelRecoStartZ;
    fSelRecoDownstreamX = fSelRecoEndX;
    fSelRecoDownstreamY = fSelRecoEndY;
    fSelRecoDownstreamZ = fSelRecoEndZ;
  }
  else{
    fSelRecoDownstreamX = fSelRecoStartX;
    fSelRecoDownstreamY = fSelRecoStartY;
    fSelRecoDownstreamZ = fSelRecoStartZ;
    fSelRecoUpstreamX = fSelRecoEndX;
    fSelRecoUpstreamY = fSelRecoEndY;
    fSelRecoUpstreamZ = fSelRecoEndZ;
  }
  fSelRecoLength = sel_track->Length();
  fSelRecoCharge = CalculateTrackCharge(sel_track, sel_track_hits);
  //24/07/18 DBrailsford Use the data product to get the neutrino energy
  //fSelRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);
  fRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fRecoEHad = energyRecoHandle->fHadLorentzVector.E();

  fSelRecoContained = energyRecoHandle->longestTrackContained; 
  if (energyRecoHandle->trackMomMethod==1){ //momentum by range was used to calculate ENu
    fSelRecoMomContained = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fRecoMomLep = fSelRecoMomContained;
  }
  else if (energyRecoHandle->trackMomMethod==0){//momentum by MCS
    fSelRecoMomMCS = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fRecoMomLep = fSelRecoMomMCS;
  }

  int g4id = FDSelectionUtils::TrueParticleID(sel_track_hits);
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
  if (matched_mcparticle){
    //Fill variables
    fSelTruePDG = matched_mcparticle->PdgCode();
    if (matched_mcparticle->Mother()==0) fSelTruePrimary = 1;
    else fSelTruePrimary = 0;
    fSelTrueMomX = matched_mcparticle->Momentum().X();
    fSelTrueMomY = matched_mcparticle->Momentum().Y();
    fSelTrueMomZ = matched_mcparticle->Momentum().Z();
    fSelTrueMomT = matched_mcparticle->Momentum().T();
    fSelTrueStartX = matched_mcparticle->Position(0).X();
    fSelTrueStartY = matched_mcparticle->Position(0).Y();
    fSelTrueStartZ = matched_mcparticle->Position(0).Z();
    fSelTrueStartT = matched_mcparticle->Position(0).T();
    fSelTrueEndX = matched_mcparticle->EndPosition().X();
    fSelTrueEndY = matched_mcparticle->EndPosition().Y();
    fSelTrueEndZ = matched_mcparticle->EndPosition().Z();
    fSelTrueEndT = matched_mcparticle->EndPosition().T();
  }
  //Now get the pid stuff
  art::FindManyP<anab::MVAPIDResult> fmpidt(trackListHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(sel_track.key());
  std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
  //Get the PIDs
  fSelMVAElectron = mvaOutMap["electron"];
  fSelMVAPion = mvaOutMap["pich"];
  fSelMVAMuon = mvaOutMap["muon"];
  fSelMVAProton = mvaOutMap["proton"];
  fSelMVAPhoton = mvaOutMap["photon"];
}

double FDSelection::NumuCutSelection::CalculateTrackCharge(art::Ptr<recob::Track> const track, std::vector< art::Ptr< recob::Hit> > const track_hits){
  double charge = 0;
  for (unsigned int i_hit = 0; i_hit < track_hits.size(); i_hit++){
    if (track_hits[i_hit]->WireID().Plane != 2) continue; //Collection only
    charge += track_hits[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(track_hits[i_hit]->PeakTime(), fT0);
  }
  return charge;
}


bool FDSelection::NumuCutSelection::IsTrackContained(art::Ptr<recob::Track> const track, std::vector< art::Ptr<recob::Hit > > const track_hits, art::Event const & evt){
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

DEFINE_ART_MODULE(FDSelection::NumuCutSelection)

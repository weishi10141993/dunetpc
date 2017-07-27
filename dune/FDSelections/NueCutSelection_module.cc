////////////////////////////////////////////////////////////////////////
// Class:       NueCutSelection
// Plugin Type: analyzer (art v2_07_03)
// File:        NueCutSelection_module.cc
//
// Generated at Tue Jun 27 06:07:56 2017 by Mike Wallbank using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

// framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

// custom
#include "FDSelectionUtils.h"
#include "PIDAnaAlg.h"

// c++
#include <iostream>

//ROOT
#include "TTree.h"

constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

namespace FDSelection {
  class NueCutSelection;
}

class FDSelection::NueCutSelection : public art::EDAnalyzer {
public:

  explicit NueCutSelection(const fhicl::ParameterSet& p);

  // Plugins should not be copied or assigned.
  NueCutSelection(NueCutSelection const &) = delete;
  NueCutSelection(NueCutSelection &&) = delete;
  NueCutSelection & operator = (NueCutSelection const &) = delete;
  NueCutSelection & operator = (NueCutSelection &&) = delete;

  // Required functions.
  void analyze(const art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginSubRun(const art::SubRun& sr) override;
  void endJob() override;
  void endSubRun(const art::SubRun& sr) override;

private:

  /// Resets all tree variables
  void Reset();

  /// Get the truth info from the art event record
  void GetTruthInfo(const art::Event& evt);

  /// Run the selection and dump relevant info to the tree
  void RunSelection(const art::Event& evt);

  // The selection tree
  TTree *fTree;
  // Generic stuff
  int fRun;
  int fSubRun;
  int fEvent;
  int fIsMC;
  // Neutrino stuff
  int fNuPdg;   // Interaction PDG
  int fBeamPdg; // PDG at point of creation
  int fNC;      // 1=is NC, 0=otherwise
  int fMode;    // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  double fQ2; 
  double fEnu; 
  double fW; // X-Sec params
  double fX;
  double fY;
  double fNuMomX; // Neutrino momentums
  double fNuMomY;
  double fNuMomZ;
  double fNuMomT;
  double fNuX; // Interaction positions
  double fNuY;
  double fNuZ;
  double fNuT;
  // Outgoing lepton stuff
  int fLepPDG;
  double fLepMomX;
  double fLepMomY;
  double fLepMomZ;
  double fLepMomT;
  double fLepNuAngle;
  // Selection stuff
  // true
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
  // reco
  double fSelRecoDirX;
  double fSelRecoDirY;
  double fSelRecoDirZ;
  double fSelRecoEnergy;
  double fSelRecoStartX;
  double fSelRecoStartY;
  double fSelRecoStartZ;
  double fSelRecodEdx;
  // MVA stuff
  double fSelMVAElectron;
  double fSelMVAPion;
  double fSelMVAMuon;
  double fSelMVAProton;
  double fSelMVAPhoton;
  double fRecoENu; //Reco neutrino energy

  // POT tree stuff
  TTree* fPOTTree;
  double fPOT;

  // Module fhicl labels
  std::string fNuGenModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fPIDModuleLabel, fPOTModuleLabel;

  // Algs
  PIDAnaAlg fPIDAnaAlg;

};

DEFINE_ART_MODULE(FDSelection::NueCutSelection)

FDSelection::NueCutSelection::NueCutSelection(const fhicl::ParameterSet& pset) : EDAnalyzer(pset),
                                                                                 fPIDAnaAlg(pset) {
  fNuGenModuleLabel  = pset.get<std::string>("NuGenModuleLabel");
  fTrackModuleLabel  = pset.get<std::string>("TrackModuleLabel");
  fShowerModuleLabel = pset.get<std::string>("ShowerModuleLabel");
  fPIDModuleLabel    = pset.get<std::string>("PIDModuleLabel");
  fPOTModuleLabel    = pset.get<std::string>("POTModuleLabel");
}

void FDSelection::NueCutSelection::analyze(const art::Event& evt) {

  // Generic stuff that can be pulled from the top of the record
  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();
  fIsMC = !(evt.isRealData());

  if (fIsMC)
    GetTruthInfo(evt);

  //RunSelection(evt);

  fPIDAnaAlg.Run(evt);

  fTree->Fill();
  Reset(); //Reset at the end of the event

  return;

}

void FDSelection::NueCutSelection::beginJob() {

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("nuecutsel","Nue cut selection");
  fTree->Branch("Run",			&fRun);
  fTree->Branch("SubRun",		&fSubRun);
  fTree->Branch("Event",		&fEvent);
  fTree->Branch("IsMC",			&fIsMC);
  fTree->Branch("NuPdg",		&fNuPdg);
  fTree->Branch("BeamPdg",		&fBeamPdg);
  fTree->Branch("NC",			&fNC);
  fTree->Branch("Mode",			&fMode);
  fTree->Branch("Q2",			&fQ2);
  fTree->Branch("Enu",			&fEnu);
  fTree->Branch("W",			&fW);
  fTree->Branch("X",			&fX);
  fTree->Branch("Y",			&fY);
  fTree->Branch("NuMomX",		&fNuMomX);
  fTree->Branch("NuMomY",		&fNuMomY);
  fTree->Branch("NuMomZ",		&fNuMomZ);
  fTree->Branch("NuMomT",		&fNuMomT);
  fTree->Branch("NuX",			&fNuX);
  fTree->Branch("NuY",			&fNuY);
  fTree->Branch("NuZ",			&fNuZ);
  fTree->Branch("NuT",			&fNuT);
  fTree->Branch("LepPDG",		&fLepPDG);
  fTree->Branch("LepMomX",		&fLepMomX);
  fTree->Branch("LepMomY",		&fLepMomY);
  fTree->Branch("LepMomZ",		&fLepMomZ);
  fTree->Branch("LepMomT",		&fLepMomT);
  fTree->Branch("LepNuAngle",		&fLepNuAngle);
  fTree->Branch("SelTruePDG",		&fSelTruePDG);
  fTree->Branch("SelTruePrimary",	&fSelTruePrimary);
  fTree->Branch("SelTrueMomX",		&fSelTrueMomX);
  fTree->Branch("SelTrueMomY",		&fSelTrueMomY);
  fTree->Branch("SelTrueMomZ",		&fSelTrueMomZ);
  fTree->Branch("SelTrueMomT",		&fSelTrueMomT);
  fTree->Branch("SelTrueStartX",	&fSelTrueStartX);
  fTree->Branch("SelTrueStartY",	&fSelTrueStartY);
  fTree->Branch("SelTrueStartZ",	&fSelTrueStartZ);
  fTree->Branch("SelTrueStartT",	&fSelTrueStartT);
  fTree->Branch("SelRecoMomX",		&fSelRecoDirX);
  fTree->Branch("SelRecoMomY",		&fSelRecoDirY);
  fTree->Branch("SelRecoMomZ",		&fSelRecoDirZ);
  fTree->Branch("SelRecoMomT",		&fSelRecoEnergy);
  fTree->Branch("SelRecoStartX",	&fSelRecoStartX);
  fTree->Branch("SelRecoStartY",	&fSelRecoStartY);
  fTree->Branch("SelRecoStartZ",	&fSelRecoStartZ);
  fTree->Branch("SelRecoStartT",	&fSelRecodEdx);
  fTree->Branch("RecoENu",		&fRecoENu);
  fTree->Branch("SelMVAElectron",	&fSelMVAElectron);
  fTree->Branch("SelMVAPion",		&fSelMVAPion);
  fTree->Branch("SelMVAMuon",		&fSelMVAMuon);
  fTree->Branch("SelMVAProton",		&fSelMVAProton);
  fTree->Branch("SelMVAPhoton",		&fSelMVAPhoton);

  fPOTTree = tfs->make<TTree>("pottree","pot tree");
  fPOTTree->Branch("POT",     &fPOT);
  fPOTTree->Branch("Run",     &fRun);
  fPOTTree->Branch("SubRun",  &fSubRun);

  Reset();

}

void FDSelection::NueCutSelection::beginSubRun(const art::SubRun& sr) {
  // Implementation of optional member function here.
}

void FDSelection::NueCutSelection::endSubRun(const art::SubRun& sr) {

  // Need the run and subrun
  fRun = sr.run();
  fSubRun = sr.subRun();
  // Need the POT (obvs) -- MW: lol 
  art::Handle<sumdata::POTSummary> potListHandle;

  if (sr.getByLabel(fPOTModuleLabel,potListHandle))
    fPOT = potListHandle->totpot;
  else
    fPOT = 0.;

  if (fPOTTree)
    fPOTTree->Fill();

}

void FDSelection::NueCutSelection::endJob() {
  // Implementation of optional member function here.
}

void FDSelection::NueCutSelection::Reset() {
  // Generic stuff
  fRun			= kDefInt;
  fSubRun		= kDefInt;
  fEvent		= kDefInt;
  fIsMC			= kDefInt;
  // Neutrino stuff
  fNuPdg		= kDefInt; 
  fBeamPdg		= kDefInt; 
  fNC			= kDefInt;    
  fMode			= kDefInt; 
  fQ2			= kDefDoub;
  fEnu			= kDefDoub;
  fW			= kDefDoub;
  fX			= kDefDoub;
  fY			= kDefDoub;
  fNuMomX		= kDefDoub;
  fNuMomY		= kDefDoub;
  fNuMomZ		= kDefDoub;
  fNuMomT		= kDefDoub;
  fNuX			= kDefDoub;
  fNuY			= kDefDoub;
  fNuZ			= kDefDoub;
  fNuT			= kDefDoub;
  // Outgoing lepton stuff
  fLepPDG		= kDefInt;
  fLepMomX		= kDefDoub;
  fLepMomY		= kDefDoub;
  fLepMomZ		= kDefDoub;
  fLepMomT		= kDefDoub;
  fLepNuAngle		= kDefDoub;
  // Selection stuff
  fRecoENu		= kDefDoub; //Neutrino reco energy
  // true
  fSelTruePDG		= kDefInt;
  fSelTruePrimary	= kDefInt;
  fSelTrueMomX		= kDefDoub;
  fSelTrueMomY		= kDefDoub;
  fSelTrueMomZ		= kDefDoub;
  fSelTrueMomT		= kDefDoub;
  fSelTrueStartX	= kDefDoub;
  fSelTrueStartY	= kDefDoub;
  fSelTrueStartZ	= kDefDoub;
  fSelTrueStartT	= kDefDoub;
  // reco
  fSelRecoDirX		= kDefDoub;
  fSelRecoDirY		= kDefDoub;
  fSelRecoDirZ		= kDefDoub;
  fSelRecoEnergy	= kDefDoub;
  fSelRecoStartX	= kDefDoub;
  fSelRecoStartY	= kDefDoub;
  fSelRecoStartZ	= kDefDoub;
  fSelRecodEdx		= kDefDoub;
  // MVA stuff
  fSelMVAElectron	= kDefDoub;
  fSelMVAPion		= kDefDoub;
  fSelMVAMuon		= kDefDoub;
  fSelMVAProton		= kDefDoub;
  fSelMVAPhoton		= kDefDoub;

}

void FDSelection::NueCutSelection::GetTruthInfo(const art::Event& evt){

  // Get the truth record
  art::Handle<std::vector<simb::MCTruth> > mcTruthHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcTruth;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthHandle)){
    art::fill_ptr_vector(mcTruth, mcTruthHandle);
  }
  else
    mf::LogWarning("NueCutSelection") << "No MCTruth";

  // Get the flux record
  art::Handle<std::vector<simb::MCFlux> > fluxHandle;
  std::vector<art::Ptr<simb::MCFlux> > flux;
  if (evt.getByLabel(fNuGenModuleLabel, fluxHandle))
    art::fill_ptr_vector(flux, fluxHandle);
  else
    mf::LogWarning("NueCutSelection") << "No MCFlux";

  // Chuck out a warning if there are multiple truths (i.e. multiple neutrinos)
  if (mcTruth.size() > 1) {
    mf::LogWarning("NueCutSelection") << "There are " << mcTruth.size() << " MCTruth in this event.  Only taking the first one!!!!";
    std::cout<<"nice"<<std::endl;
  }

  // Get the neutrino
  const simb::MCNeutrino& mcneutrino = mcTruth.at(0)->GetNeutrino();
  const simb::MCParticle& nu = mcneutrino.Nu();
  const simb::MCParticle& lepton = mcneutrino.Lepton();

  // Neutrino stuff
  fNuPdg	= nu.PdgCode();
  fBeamPdg	= flux.at(0)->fntype;
  fNC		= mcneutrino.CCNC();
  fMode		= mcneutrino.Mode();
  fQ2		= mcneutrino.QSqr();
  fEnu		= nu.E(); 
  fW		= mcneutrino.W(); 
  fX		= mcneutrino.X();
  fY		= mcneutrino.Y();
  fNuMomX	= nu.Momentum().X();
  fNuMomY	= nu.Momentum().Y();
  fNuMomZ	= nu.Momentum().Z();
  fNuMomT	= nu.Momentum().T();
  fNuX		= nu.Position().X();
  fNuY		= nu.Position().Y();
  fNuZ		= nu.Position().Z();
  fNuT		= nu.Position().Y();
  // Outgoing lepton stuff
  fLepPDG	= lepton.PdgCode();
  fLepMomX	= lepton.Momentum().X();
  fLepMomY	= lepton.Momentum().Y();
  fLepMomZ	= lepton.Momentum().Z();
  fLepMomT	= lepton.Momentum().T();
  fLepNuAngle	= nu.Momentum().Vect().Angle(lepton.Momentum().Vect());

  return;

}

void FDSelection::NueCutSelection::RunSelection(art::Event const & evt) {

  // Get showers from the event
  art::Handle<std::vector<recob::Shower> > showerHandle;
  std::vector<art::Ptr<recob::Shower> > showers;
  if (evt.getByLabel(fShowerModuleLabel,showerHandle))
    art::fill_ptr_vector(showers,showerHandle);

  // Get the highest energy shower
  int i_energyist_shower = -1;
  double energy_energyist_shower = -999;
  for (std::vector<art::Ptr<recob::Shower> >::const_iterator showerIt = showers.begin(); showerIt != showers.end(); ++showerIt) {
    if ((*showerIt)->Energy().size() == 0)
      continue;
    double current_energy = (*showerIt)->Energy().at((*showerIt)->best_plane());
    if (current_energy > energy_energyist_shower) {
      energy_energyist_shower = current_energy;
      i_energyist_shower = showerIt->key();
    }
  }
  if (i_energyist_shower < 0)
    return;

  const art::Ptr<recob::Shower> sel_shower = showers.at(i_energyist_shower);
  fSelRecoDirX   = sel_shower->Direction().X();
  fSelRecoDirY   = sel_shower->Direction().Y();
  fSelRecoDirZ   = sel_shower->Direction().Z();
  fSelRecoStartX = sel_shower->ShowerStart().X();
  fSelRecoStartY = sel_shower->ShowerStart().Y();
  fSelRecoStartZ = sel_shower->ShowerStart().Z();
  if (sel_shower->Energy().size() != 0) {
    fSelRecoEnergy = sel_shower->Energy().at(sel_shower->best_plane());
    fSelRecodEdx   = sel_shower->dEdx().at(sel_shower->best_plane());
  }

  fRecoENu = kDefDoub; //temp

  // shower hits
  art::FindManyP<recob::Hit> fmhs(showerHandle, evt, fShowerModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > hits = fmhs.at(i_energyist_shower);
  int g4id = FDSelectionUtils::TrueParticleID(hits);
  art::ServiceHandle<cheat::BackTracker> bt;
  const simb::MCParticle* matched_mcparticle = bt->ParticleList().at(g4id);
  if (matched_mcparticle) {
    fSelTruePDG = matched_mcparticle->PdgCode();
    if (matched_mcparticle->Mother() == 0)
      fSelTruePrimary = 1;
    else
      fSelTruePrimary = 0;
    fSelTrueMomX   = matched_mcparticle->Momentum().X();
    fSelTrueMomY   = matched_mcparticle->Momentum().Y();
    fSelTrueMomZ   = matched_mcparticle->Momentum().Z();
    fSelTrueMomT   = matched_mcparticle->Momentum().T();
    fSelTrueStartX = matched_mcparticle->Position(0).X();
    fSelTrueStartY = matched_mcparticle->Position(0).Y();
    fSelTrueStartZ = matched_mcparticle->Position(0).Z();
    fSelTrueStartT = matched_mcparticle->Position(0).T();
  }

  // pid stuff
  art::FindManyP<anab::MVAPIDResult> fmpids(showerHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpids.at(i_energyist_shower);
  std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
  fSelMVAElectron = mvaOutMap["electron"];
  fSelMVAPion     = mvaOutMap["pich"];
  fSelMVAMuon     = mvaOutMap["muon"];
  fSelMVAProton   = mvaOutMap["proton"];
  fSelMVAPhoton   = mvaOutMap["photon"];
  // reco stuff

}

////////////////////////////////////////////////////////////////////////
// Class:       NueCutSelection
// Plugin Type: analyzer (art v2_07_03)
// File:        NueCutSelection_module.cc
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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"

//DUNE
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

//Custom
#include "PIDAnaAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoShowerSelector.h"


constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

namespace FDSelection {
  class NueCutSelection;
}


class FDSelection::NueCutSelection : public art::EDAnalyzer {
public:
  explicit NueCutSelection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NueCutSelection(NueCutSelection const &) = delete;
  NueCutSelection(NueCutSelection &&) = delete;
  NueCutSelection & operator = (NueCutSelection const &) = delete;
  NueCutSelection & operator = (NueCutSelection &&) = delete;

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
  double CalculateShowerCharge(art::Ptr<recob::Shower> const track, std::vector< art::Ptr< recob::Hit> > const track_hits); //Calculate hit charge on track as measured by the collection plane

  // Declare member data here.

  //Algs
  PIDAnaAlg fPIDAnaAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  shower::ShowerEnergyAlg fShowerEnergyAlg;


  //Tools
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
  double fSelRecoDirX;
  double fSelRecoDirY;
  double fSelRecoDirZ;
  double fSelRecoStartX;
  double fSelRecoStartY;
  double fSelRecoStartZ;
  double fSelRecodEdx;
  double fSelRecoCharge;

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
  std::string fShowerModuleLabel;
  std::string fPIDModuleLabel;
  std::string fHitsModuleLabel;
  std::string fPOTModuleLabel;
  std::string fEnergyRecoModuleLabel;
 

};


FDSelection::NueCutSelection::NueCutSelection(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)   ,
  fPIDAnaAlg(pset.get<fhicl::ParameterSet>("ModuleLabels"))   ,
  fCalorimetryAlg          (pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(pset.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))},
  fNuGenModuleLabel        (pset.get< std::string >("ModuleLabels.NuGenModuleLabel")),
  fShowerModuleLabel        (pset.get< std::string >("ModuleLabels.ShowerModuleLabel")),
  fPIDModuleLabel          (pset.get< std::string >("ModuleLabels.PIDModuleLabel")),
  fHitsModuleLabel         (pset.get< std::string >("ModuleLabels.HitsModuleLabel")),
  fPOTModuleLabel          (pset.get< std::string >("ModuleLabels.POTModuleLabel")),
  fEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.EnergyRecoModuleLabel"))
{}

void FDSelection::NueCutSelection::analyze(art::Event const & evt)
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

void FDSelection::NueCutSelection::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("nuecutsel","Nue cut selection");
    fTree->Branch("Run",&fRun);
    fTree->Branch("SubRun",&fSubRun);
    fTree->Branch("Event",&fEvent);
    fTree->Branch("IsMC",&fIsMC);
    fTree->Branch("T0",&fT0);
    fTree->Branch("NuPdg",&fNuPdg);
    fTree->Branch("BeamPdg",&fBeamPdg);
    fTree->Branch("NC",&fNC);
    fTree->Branch("Mode",&fMode);
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
    fTree->Branch("SelRecoDirX",		&fSelRecoDirX);
    fTree->Branch("SelRecoDirY",		&fSelRecoDirY);
    fTree->Branch("SelRecoDirZ",		&fSelRecoDirZ);
    fTree->Branch("SelRecoStartX",	&fSelRecoStartX);
    fTree->Branch("SelRecoStartY",	&fSelRecoStartY);
    fTree->Branch("SelRecoStartZ",	&fSelRecoStartZ);
    fTree->Branch("SelRecodEdx",	&fSelRecodEdx);
    fTree->Branch("SelRecoCharge",&fSelRecoCharge);
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

void FDSelection::NueCutSelection::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}
void FDSelection::NueCutSelection::endSubRun(const art::SubRun& sr){
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

void FDSelection::NueCutSelection::endJob()
{
  // Implementation of optional member function here.
}

void FDSelection::NueCutSelection::Reset()
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
  fSelRecoDirX		= kDefDoub;
  fSelRecoDirY		= kDefDoub;
  fSelRecoDirZ		= kDefDoub;
  fSelRecoStartX	= kDefDoub;
  fSelRecoStartY	= kDefDoub;
  fSelRecoStartZ	= kDefDoub;
  fSelRecodEdx		= kDefDoub;
  fSelRecoCharge = kDefDoub;

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

void FDSelection::NueCutSelection::GetEventInfo(art::Event const & evt){
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

void FDSelection::NueCutSelection::GetTruthInfo(art::Event const & evt){
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
    mf::LogWarning("NueCutSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one!!!!";
  }
  for (unsigned int i_mctruth = 0; i_mctruth < mcList.size(); i_mctruth++){
    fNuPdg    = mcList[i_mctruth]->GetNeutrino().Nu().PdgCode();
    fBeamPdg  = mcFlux[i_mctruth]->fntype;
    fNC       = mcList[i_mctruth]->GetNeutrino().CCNC();
    fMode     = mcList[i_mctruth]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
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
  }
}

void FDSelection::NueCutSelection::RunSelection(art::Event const & evt){
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
  if (!(sel_shower.isAvailable())) return;


  //24/07/18 DBrailsford Get the reco energy data product for neutrinos
  art::Handle<dune::EnergyRecoOutput> energyRecoHandle;
  if (!evt.getByLabel(fEnergyRecoModuleLabel, energyRecoHandle)) return;

  //Get the hits for said track
  art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > sel_shower_hits = fmhs.at(sel_shower.key());

  fSelRecoDirX   = sel_shower->Direction().X();
  fSelRecoDirY   = sel_shower->Direction().Y();
  fSelRecoDirZ   = sel_shower->Direction().Z();
  fSelRecoStartX = sel_shower->ShowerStart().X();
  fSelRecoStartY = sel_shower->ShowerStart().Y();
  fSelRecoStartZ = sel_shower->ShowerStart().Z();

  fSelRecoCharge = CalculateShowerCharge(sel_shower, sel_shower_hits);
  //24/07/18 DBrailsford Use the data product to get the neutrino energy
  //fSelRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);
  fRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fRecoEHad = energyRecoHandle->fHadLorentzVector.E();
  fRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

  int g4id = FDSelectionUtils::TrueParticleID(sel_shower_hits);
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
  art::FindManyP<anab::MVAPIDResult> fmpidt(showerListHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(sel_shower.key());
  std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
  //Get the PIDs
  fSelMVAElectron = mvaOutMap["electron"];
  fSelMVAPion = mvaOutMap["pich"];
  fSelMVAMuon = mvaOutMap["muon"];
  fSelMVAProton = mvaOutMap["proton"];
  fSelMVAPhoton = mvaOutMap["photon"];
}

double FDSelection::NueCutSelection::CalculateShowerCharge(art::Ptr<recob::Shower> const shower, std::vector< art::Ptr< recob::Hit> > const shower_hits){
  double charge = 0;
  for (unsigned int i_hit = 0; i_hit < shower_hits.size(); i_hit++){
    if (shower_hits[i_hit]->WireID().Plane != 2) continue; //Collection only
    charge += shower_hits[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(shower_hits[i_hit]->PeakTime(), fT0);
  }
  return charge;
}




DEFINE_ART_MODULE(FDSelection::NueCutSelection)

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
//LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"




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
  void endJob() override;

private:

  //Delcare private functions
  void Reset();      //Resets all tree vars

  void GetTruthInfo(art::Event const & evt);  //Grab the truth info from the art record
  // Declare member data here.

  TTree *fTree; //The selection tree
  //Generic stuff
  int fRun;
  int fSubRun;
  int fEvent;
  int fIsMC;
  //Neutrino stuff
  int fNuPdg; //Interaction PDG
  int fBeamPdg; //PDG at point of creation
  int fNC;    // 1=is NC, 0=otherwise
  int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  double fQ2; 
  double fEnu; 
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
  double fLepMomX;
  double fLepMomY;
  double fLepMomZ;
  double fLepMomT;
  double fLepNuAngle;

  //Fhicl pset labels
  std::string fNuGenModuleLabel;

 
};


FDSelection::NumuCutSelection::NumuCutSelection(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)   ,
  fNuGenModuleLabel        (pset.get< std::string >("NuGenModuleLabel"))
{}

void FDSelection::NumuCutSelection::analyze(art::Event const & evt)
{
  //Get the generic stuff that can be pulled from the top of the record
  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();
  fIsMC = !(evt.isRealData());

  if (fIsMC) GetTruthInfo(evt);

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
    fTree->Branch("NuPdg",&fNuPdg);
    fTree->Branch("BeamPdg",&fBeamPdg);
    fTree->Branch("NC",&fNC);
    fTree->Branch("Mode",&fMode);
    fTree->Branch("Q2",&fQ2);
    fTree->Branch("Enu",&fEnu);
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
    fTree->Branch("LepMomX",&fLepMomX);
    fTree->Branch("LepMomY",&fLepMomY);
    fTree->Branch("LepMomZ",&fLepMomZ);
    fTree->Branch("LepMomT",&fLepMomT);
    fTree->Branch("LepNuAngle",&fLepNuAngle);

    Reset();  //Default value all variables now
}

void FDSelection::NumuCutSelection::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
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
  //Neutrino stuff
  fNuPdg = kDefInt; 
  fBeamPdg = kDefInt; 
  fNC = kDefInt;    
  fMode = kDefInt; 
  fQ2 = kDefDoub; 
  fEnu = kDefDoub; 
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
  fLepMomX = kDefDoub;
  fLepMomY = kDefDoub;
  fLepMomZ = kDefDoub;
  fLepMomT = kDefDoub;
  fLepNuAngle = kDefDoub;
}

void FDSelection::NumuCutSelection::GetTruthInfo(art::Event const & evt){
  //Get the truth record
  art::Handle< std::vector<simb::MCTruth> > mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle)){
    art::fill_ptr_vector(mcList, mcTruthListHandle);
  }
  //Chuck out a warning if there are multiple truths (i.e. multiple neutrinos)
  if (mcList.size() > 1){
    mf::LogWarning("NumuCutSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one!!!!";
  }
}

DEFINE_ART_MODULE(FDSelection::NumuCutSelection)

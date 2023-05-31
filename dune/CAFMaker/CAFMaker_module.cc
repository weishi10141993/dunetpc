
////////////////////////////////////////////////////////////////////////
//
// \file CAFMaker_module.cc
//
// Chris Marshall's version
// Largely based on historical FDSensOpt/CAFMaker_module.cc
//
///////////////////////////////////////////////////////////////////////

#ifndef CAFMaker_H
#define CAFMaker_H

// Generic C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

//#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "dune/FDSensOpt/FDSensOptData/MVASelectPID.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/Result.h"
#include "dune/RegCVN/func/RegCVNResult.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// dunerw stuff
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/exceptions.hh"
//#include "systematicstools/utility/md5.hh"

// root
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// pdg
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

// genie
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"

// custom
#include "dune/FDSelections/FDSelectionData/PandSelectParams.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"

constexpr int knShifts = 100; // number of shifts
constexpr int kmaxRwgts = 100; // Largest number of reweights in a shift

namespace dunemva {

  class CAFMaker : public art::EDAnalyzer {

    public:

      explicit CAFMaker(fhicl::ParameterSet const& pset);
      virtual ~CAFMaker();
      void beginJob() override;
      void beginSubRun(const art::SubRun& sr) override;
      void endSubRun(const art::SubRun& sr) override;
      void reconfigure(fhicl::ParameterSet const& pset) /*override*/;
      void analyze(art::Event const & evt) override;


    private:
      void resetCAFvars();

      std::string fMVASelectLabel;
      std::string fMVASelectNueLabel;
      std::string fMVASelectNumuLabel;
      std::string fPandSelectParamsLabel;

      std::string fCVNLabel;
      std::string fRegCVNLabel;

      std::string fEnergyRecoNueLabel;
      std::string fEnergyRecoNumuLabel;
      std::string fMVAMethod;

      std::string fHitsModuleLabel;

      float fOscPro;
      double fWeight;
      TTree* fTree;  
      TTree* fPOT;  

      int runPOTTreeVariable, subrunPOTTreeVariable;
      double potPOTTreeVariable;

      // Get reweight knobs from fhicl file -- no hard-coded shifts
      int fNwgt[knShifts];
      double fCvWgts[knShifts];
      double fWgts[knShifts][kmaxRwgts];
      bool fisWgtSyst[knShifts];

      // CAF variables
      // configuration variables
      int fIsFD, fIsFHC;
      // event accounting
      int fRun, fSubrun, fEvent;
      // Truth information
      int fIsCC, fNuPDG, fNuPDGunosc, fMode, fLepPDG; 
      double fEv, fQ2, fW, fX, fY, fNuMomX, fNuMomY, fNuMomZ, fLepMomX, fLepMomY, fLepMomZ, fLepE, fLepNuAngle;
      // True particle counts
      int nP, nN, nPip, nPim, nPi0, nKp, nKm, nK0, nEM, nOtherHad, nNucleus, nUNKNOWN;
      double eP, eN, ePip, ePim, ePi0, eOther;
      double eRecoP, eRecoN, eRecoPip, eRecoPim, eRecoPi0, eRecoOther;
      double eDepP, eDepN, eDepPip, eDepPim, eDepPi0, eDepOther;

      // Reco information
      double fErecoNue;
      double fRecoLepEnNue;
      double fRecoHadEnNue;
      int fRecoMethodNue; // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
      double fRecoLepAngNue;
      double fErecoNumu; 
      double fRecoLepEnNumu;
      double fRecoHadEnNumu;
      int fRecoMethodNumu; // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
      double fRecoLepAngNumu;
      int fLongestTrackContNumu; // 1 = contained, 0 = exiting, -1 = not set
      int fTrackMomMethodNumu; // 1 = range, 0 = MCS, -1 = not set

      double fMVAResult;
      double fMVAResultNue;
      double fMVAResultNumu;

      double fSelTrackPandizzleScore;
      double fSelShowerPandrizzleScore;
      double fSelShowerJamPandrizzleScore;

      double fNuePandrizzleCut;
      double fNueJamPandrizzleCut;
      double fNuePandizzleCut;
      double fNumuPandizzleCut;

      double fAnuePandrizzleCut;
      double fAnueJamPandrizzleCut;
      double fAnuePandizzleCut;
      double fAnumuPandizzleCut;

      bool fUseFHCCut;

      double fNuVtxX;
      double fNuVtxY;
      double fNuVtxZ;

      double fRecoNuVtxX;
      double fRecoNuVtxY;
      double fRecoNuVtxZ;

      // CVN outputs
      double fCVNResultIsAntineutrino;
      double fCVNResultNue, fCVNResultNumu, fCVNResultNutau, fCVNResultNC; // flavour
      double fCVNResult0Protons, fCVNResult1Protons, fCVNResult2Protons, fCVNResultNProtons; // #protons
      double fCVNResult0Pions, fCVNResult1Pions, fCVNResult2Pions, fCVNResultNPions; // #pions
      double fCVNResult0Pizeros, fCVNResult1Pizeros, fCVNResult2Pizeros, fCVNResultNPizeros; // #pizeros
      double fCVNResult0Neutrons, fCVNResult1Neutrons, fCVNResult2Neutrons, fCVNResultNNeutrons; // #neutrons

      double fRegCVNNueE;

    //int meta_run, meta_subrun;
      double meta_pot;

      systtools::provider_list_t fSystProviders;

      calo::CalorimetryAlg fCaloAlg;
      double fRecombFactor;

  }; // class CAFMaker


  //------------------------------------------------------------------------------
  CAFMaker::CAFMaker(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset), fCaloAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  {
    this->reconfigure(pset);
  }

  dunemva::CAFMaker::~CAFMaker(){}

  void CAFMaker::resetCAFvars()
  {
    fIsFD = -999;
    fIsFHC = -999; // ?
    fRun = -999;
    fSubrun = -999;
    fEvent = -999;
    fIsCC = -999;
    fNuPDG = -999;
    fNuPDGunosc = -999;
    fMode = -999;
    fLepPDG = -999;
    fEv = -999;
    fQ2 = -999;
    fW = -999;
    fX = -999;
    fY = -999;
    fNuMomX= -999;
    fNuMomY= -999;
    fNuMomZ= -999;
    fLepMomX= -999;
    fLepMomY= -999;
    fLepMomZ= -999;
    fLepE= -999;
    fLepNuAngle = -999;
    nP= -999;
    nN= -999;
    nPip= -999;
    nPim= -999;
    nPi0= -999;
    nKp= -999;
    nKm= -999;
    nK0= -999;
    nEM= -999;
    nOtherHad= -999;
    nNucleus= -999;
    nUNKNOWN = -999;
    eP= -999;
    eN= -999;
    ePip= -999;
    ePim= -999;
    ePi0= -999;
    eOther = -999;
    eRecoP= -999;
    eRecoN= -999;
    eRecoPip= -999;
    eRecoPim= -999;
    eRecoPi0= -999;
    eRecoOther = -999;
    eDepP= -999;
    eDepN= -999;
    eDepPip= -999;
    eDepPim= -999;
    eDepPi0= -999;
    eDepOther = -999;
    fNuVtxX = -999;
    fNuVtxY = -999;
    fNuVtxZ = -999;

    fErecoNue= -999;
    fRecoLepEnNue= -999;
    fRecoHadEnNue = -999;
    fRecoMethodNue = -999;
    fRecoLepAngNue = -999;
    fErecoNumu = -999; 
    fRecoLepEnNumu= -999;
    fRecoHadEnNumu = -999;
    fRecoMethodNumu = -999;
    fRecoLepAngNumu = -999;
    fLongestTrackContNumu = -999;
    fTrackMomMethodNumu= -999;
    fMVAResult= -999;
    fMVAResultNue= -999;
    fMVAResultNumu= -999;
    fCVNResultIsAntineutrino= -999;
    fCVNResultNue = -999;
    fCVNResultNumu = -999;
    fCVNResultNutau = -999;
    fCVNResultNC= -999;
    fCVNResult0Protons = -999;
    fCVNResult1Protons = -999;
    fCVNResult2Protons = -999;
    fCVNResultNProtons = -999;
    fCVNResult0Pions = -999;
    fCVNResult1Pions = -999;
    fCVNResult2Pions = -999;
    fCVNResultNPions= -999;
    fCVNResult0Pizeros = -999;
    fCVNResult1Pizeros = -999;
    fCVNResult2Pizeros = -999;
    fCVNResultNPizeros= -999;
    fCVNResult0Neutrons = -999;
    fCVNResult1Neutrons = -999;
    fCVNResult2Neutrons = -999;
    fCVNResultNNeutrons = -999;
    fRegCVNNueE = -999;

    fRecoNuVtxX = -999;
    fRecoNuVtxY = -999;
    fRecoNuVtxZ = -999;

    fSelTrackPandizzleScore = -999;
    fSelShowerPandrizzleScore = -999;
    fSelShowerJamPandrizzleScore = -999;

    // Reweight variables
    for( int k_it = 0; k_it < knShifts; ++k_it ) {
      fNwgt[k_it] = 7; // CAFAna assumes -3,-2,-1,0,1,2,3
      fCvWgts[k_it] = fisWgtSyst[k_it] ? 1. : 0; // Non weight systs default to 0.
      for( int r_it = 0; r_it < kmaxRwgts; ++r_it ) {
        fWgts[k_it][r_it] = fisWgtSyst[r_it] ? 1. : 0;
      }
    }
  }

  //------------------------------------------------------------------------------
  void CAFMaker::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fMVASelectLabel = pset.get<std::string>("MVASelectLabel");
    fMVASelectNueLabel = pset.get<std::string>("MVASelectNueLabel");
    fMVASelectNumuLabel = pset.get<std::string>("MVASelectNumuLabel");
    fCVNLabel = pset.get<std::string>("CVNLabel");
    fRegCVNLabel = pset.get<std::string>("RegCVNLabel");
    fPandSelectParamsLabel = pset.get<std::string>("PandSelectParamsLabel");
    fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
    fEnergyRecoNumuLabel = pset.get<std::string>("EnergyRecoNumuLabel");

    fHitsModuleLabel = pset.get<std::string>("HitsModuleLabel");

    // Get DUNErw stuff from its fhicl, which should be included on the CAFMaker config file
    //if( !pset.has_key("generated_systematic_provider_configuration") ) {
    //  std::cout << "[ERROR]: Could not find producer key: "
    //               "\"generated_systematic_provider_configuration\". This should "
    //               "contain a list of configured systematic providers generated by "
    //               "GenerateSystProviderConfig." << std::endl;
    //  return;
    //}

    fhicl::ParameterSet syst_provider_config = pset.get<fhicl::ParameterSet>("generated_systematic_provider_configuration");

    fSystProviders = systtools::ConfigureISystProvidersFromParameterHeaders(syst_provider_config);

    fRecombFactor = pset.get<double>("RecombFactor");

    fNuePandrizzleCut = pset.get<double>("NuePandrizzleCut");
    fNueJamPandrizzleCut = pset.get<double>("NueJamPandrizzleCut");
    fNuePandizzleCut = pset.get<double>("NuePandizzleCut");
    fNumuPandizzleCut = pset.get<double>("NumuPandizzleCut");
    fAnuePandrizzleCut = pset.get<double>("AnuePandrizzleCut");
    fAnueJamPandrizzleCut = pset.get<double>("AnueJamPandrizzleCut");
    fAnuePandizzleCut = pset.get<double>("AnuePandizzleCut");
    fAnumuPandizzleCut = pset.get<double>("AnumuPandizzleCut");
    fUseFHCCut = pset.get<bool>("UseFHCCut");
  }


  //------------------------------------------------------------------------------
  void CAFMaker::beginJob()
  {

    art::ServiceHandle<art::TFileService> tfs;
    fTree =tfs->make<TTree>("caf", "caf");

    fPOT = tfs->make<TTree>("pottree","pot tree");
    fPOT->Branch("pot", &potPOTTreeVariable,"pot/D");
    fPOT->Branch("run", &runPOTTreeVariable,"run/I");
    fPOT->Branch("subrun",&subrunPOTTreeVariable,"subrun/I");

    // book-keeping
    fTree->Branch("run",         &fRun,        "run/I");
    fTree->Branch("subrun",      &fSubrun,     "subrun/I");
    fTree->Branch("event",       &fEvent,      "event/I");
    fTree->Branch("isFD",        &fIsFD,       "isFD/I");
    fTree->Branch("isFHC",       &fIsFHC,      "isFHC/I");
    fTree->Branch("isCC",        &fIsCC,       "isCC/I");

    // true interaction quantities
    fTree->Branch("nuPDG",        &fNuPDG,        "nuPDG/I");
    fTree->Branch("nuPDGunosc",   &fNuPDGunosc,   "nuPDGunosc/I");
    fTree->Branch("NuMomX",       &fNuMomX,       "NuMomX/D");
    fTree->Branch("NuMomY",       &fNuMomY,       "NuMomY/D");
    fTree->Branch("NuMomZ",       &fNuMomZ,       "NuMomZ/D");
    fTree->Branch("Ev",           &fEv,           "Ev/D");
    fTree->Branch("mode",         &fMode,         "mode/I");
    fTree->Branch("LepPDG",       &fLepPDG,       "LepPDG/I");
    fTree->Branch("LepMomX",      &fLepMomX,      "LepMomX/D");
    fTree->Branch("LepMomY",      &fLepMomY,      "LepMomY/D");
    fTree->Branch("LepMomZ",      &fLepMomZ,      "LepMomZ/D");
    fTree->Branch("LepE",         &fLepE,         "LepE/D");
    fTree->Branch("LepNuAngle",   &fLepNuAngle,   "LepNuAngle/D");
    fTree->Branch("Q2",           &fQ2,           "Q2/D");
    fTree->Branch("W",            &fW,            "W/D");
    fTree->Branch("X",            &fX,            "X/D");
    fTree->Branch("Y",            &fY,            "Y/D");

    // FS particle counts
    fTree->Branch("nP",        &nP,         "nP/I");
    fTree->Branch("nN",        &nN,         "nN/I");
    fTree->Branch("nipip",     &nPip,       "nipip/I");
    fTree->Branch("nipim",     &nPim,       "nipim/I");
    fTree->Branch("nipi0",     &nPi0,       "nipi0/I");
    fTree->Branch("nikp",      &nKp,        "nikp/I");
    fTree->Branch("nikm",      &nKm,        "nikm/I");
    fTree->Branch("nik0",      &nK0,        "nik0/I");
    fTree->Branch("niem",      &nEM,        "niem/I");
    fTree->Branch("niother",   &nOtherHad,  "niother/I");
    fTree->Branch("nNucleus",  &nNucleus,   "nNucleus/I");
    fTree->Branch("nUNKNOWN",  &nUNKNOWN,   "nUNKNOWN/I");
    fTree->Branch("eP",        &eP,         "eP/D");
    fTree->Branch("eN",        &eN,         "eN/D");
    fTree->Branch("ePip",      &ePip,       "ePip/D");
    fTree->Branch("ePim",      &ePim,       "ePim/D");
    fTree->Branch("ePi0",      &ePi0,       "ePi0/D");
    fTree->Branch("eOther",    &eOther,     "eOther/D");
    fTree->Branch("eRecoP",        &eRecoP,         "eRecoP/D");
    fTree->Branch("eRecoN",        &eRecoN,         "eRecoN/D");
    fTree->Branch("eRecoPip",      &eRecoPip,       "eRecoPip/D");
    fTree->Branch("eRecoPim",      &eRecoPim,       "eRecoPim/D");
    fTree->Branch("eRecoPi0",      &eRecoPi0,       "eRecoPi0/D");
    fTree->Branch("eRecoOther",    &eRecoOther,     "eRecoOther/D");
    fTree->Branch("eDepP",        &eDepP,         "eDepP/D");
    fTree->Branch("eDepN",        &eDepN,         "eDepN/D");
    fTree->Branch("eDepPip",      &eDepPip,       "eDepPip/D");
    fTree->Branch("eDepPim",      &eDepPim,       "eDepPim/D");
    fTree->Branch("eDepPi0",      &eDepPi0,       "eDepPi0/D");
    fTree->Branch("eDepOther",    &eDepOther,     "eDepOther/D");

    // Reco variables

    fTree->Branch("cvnisantineutrino", &fCVNResultIsAntineutrino, "cvnisantineutrino/D");
    fTree->Branch("cvnnue",            &fCVNResultNue,            "cvnnue/D");
    fTree->Branch("cvnnumu",           &fCVNResultNumu,           "cvnnumu/D");
    fTree->Branch("cvnnutau",          &fCVNResultNutau,          "cvnnutau/D");
    fTree->Branch("cvnnc",             &fCVNResultNC,             "cvnnc/D");
    fTree->Branch("cvn0protons",       &fCVNResult0Protons,       "cvn0protons/D");
    fTree->Branch("cvn1protons",       &fCVNResult1Protons,       "cvn1protons/D");
    fTree->Branch("cvn2protons",       &fCVNResult2Protons,       "cvn2protons/D");
    fTree->Branch("cvnNprotons",       &fCVNResultNProtons,       "cvnNprotons/D");
    fTree->Branch("cvn0pions",         &fCVNResult0Pions,         "cvn0pions/D");
    fTree->Branch("cvn1pions",         &fCVNResult1Pions,         "cvn1pions/D");
    fTree->Branch("cvn2pions",         &fCVNResult2Pions,         "cvn2pions/D");
    fTree->Branch("cvnNpions",         &fCVNResultNPions,         "cvnNpions/D");
    fTree->Branch("cvn0pizeros",       &fCVNResult0Pizeros,       "cvn0pizeros/D");
    fTree->Branch("cvn1pizeros",       &fCVNResult1Pizeros,       "cvn1pizeros/D");
    fTree->Branch("cvn2pizeros",       &fCVNResult2Pizeros,       "cvn2pizeros/D");
    fTree->Branch("cvnNpizeros",       &fCVNResultNPizeros,       "cvnNpizeros/D");
    fTree->Branch("cvn0neutrons",      &fCVNResult0Neutrons,      "cvn0neutrons/D");
    fTree->Branch("cvn1neutrons",      &fCVNResult1Neutrons,      "cvn1neutrons/D");
    fTree->Branch("cvn2neutrons",      &fCVNResult2Neutrons,      "cvn2neutrons/D");
    fTree->Branch("cvnNneutrons",      &fCVNResultNNeutrons,      "cvnNneutrons/D");

    fTree->Branch("RegCVNNueE",  &fRegCVNNueE,   "RegCVNNueE/D");
    fTree->Branch("weight",      &fWeight,     "weight/D");
    fTree->Branch("oscpro",      &fOscPro,     "oscpro/F");

    fTree->Branch("Ev_reco_nue",      &fErecoNue,        "Ev_reco_nue/D");
    fTree->Branch("RecoLepEnNue",     &fRecoLepEnNue,    "RecoLepEnNue/D");
    fTree->Branch("RecoHadEnNue",     &fRecoHadEnNue,    "RecoHadEnNue/D");
    fTree->Branch("RecoMethodNue",    &fRecoMethodNue,   "RecoMethodNue/I");
    fTree->Branch("RecoLepAngNue",    &fRecoLepAngNue,   "RecoLepAngNue/D");
    fTree->Branch("Ev_reco_numu",     &fErecoNumu,       "Ev_reco_numu/D");
    fTree->Branch("RecoLepEnNumu",    &fRecoLepEnNumu,   "RecoLepEnNumu/D");
    fTree->Branch("RecoHadEnNumu",    &fRecoHadEnNumu,   "RecoHadEnNumu/D");
    fTree->Branch("RecoMethodNumu",   &fRecoMethodNumu,  "RecoMethodNumu/I");
    fTree->Branch("RecoLepAngNumu",    &fRecoLepAngNumu,   "RecoLepAngNumu/D");
    fTree->Branch("LongestTrackContNumu",  &fLongestTrackContNumu, "LongestTrackContNumu/I");
    fTree->Branch("TrackMomMethodNumu",    &fTrackMomMethodNumu,   "TrackMomMethodNumu/I");

    fTree->Branch("SelTrackPandizzleScore",    &fSelTrackPandizzleScore,   "SelTrackPandizzleScore/D");
    fTree->Branch("SelShowerPandrizzleScore",    &fSelShowerPandrizzleScore,   "SelShowerPandrizzleScore/D");
    fTree->Branch("SelShowerJamPandrizzleScore",    &fSelShowerJamPandrizzleScore,   "SelShowerJamPandrizzleScore/D");

    fTree->Branch("NuePandrizzleCut", &fNuePandrizzleCut, "NuePandrizzleCut/D");
    fTree->Branch("NueJamPandrizzleCut", &fNueJamPandrizzleCut, "NueJamPandrizzleCut/D");
    fTree->Branch("NuePandizzleCut", &fNuePandizzleCut, "NuePandizzleCut/D");
    fTree->Branch("NumuPandizzleCut", &fNumuPandizzleCut, "NumuPandizzleCut/D");

    fTree->Branch("AnuePandrizzleCut", &fAnuePandrizzleCut, "AnuePandrizzleCut/D");
    fTree->Branch("AnueJamPandrizzleCut", &fAnueJamPandrizzleCut, "AnueJamPandrizzleCut/D");
    fTree->Branch("AnuePandizzleCut", &fAnuePandizzleCut, "AnuePandizzleCut/D");
    fTree->Branch("AnumuPandizzleCut", &fAnumuPandizzleCut, "AnumuPandizzleCut/D");

    fTree->Branch("vtx_x", &fNuVtxX);
    fTree->Branch("vtx_y", &fNuVtxY);
    fTree->Branch("vtx_z", &fNuVtxZ);

    fTree->Branch("RecoVertex_x", &fRecoNuVtxX);
    fTree->Branch("RecoVertex_y", &fRecoNuVtxY);
    fTree->Branch("RecoVertex_z", &fRecoNuVtxZ);

    fTree->Branch("totpot",       &meta_pot,       "totpot/D");

    std::fill_n(fisWgtSyst, knShifts, true);
    // make DUNErw variables
    for( auto &sp : fSystProviders ) {
      systtools::SystMetaData metaData = sp->GetSystMetaData();
      for( systtools::SystMetaData::iterator itMeta = metaData.begin(); itMeta != metaData.end(); ++itMeta ) {      
        systtools::SystParamHeader head = *itMeta;
        std::string name = head.prettyName;
        unsigned int parId = head.systParamId;
        fisWgtSyst[parId] = head.isWeightSystematicVariation;
        std::string wgt_var = fisWgtSyst[parId] ? "wgt" : "var";
        std::cout << "Adding reweight branch " << parId << " for " << name << " with " << head.paramVariations.size() << " shifts" << std::endl;
        fTree->Branch( Form("%s_nshifts", name.c_str()), &fNwgt[parId], Form("%s_nshifts/I", name.c_str()) );
        fTree->Branch( Form("%s_cv%s", name.c_str(),wgt_var.c_str()), &fCvWgts[parId], Form("%s_cv%s/D", name.c_str(),wgt_var.c_str()) );
        fTree->Branch( Form("%s_%s", wgt_var.c_str(),name.c_str()), fWgts[parId], Form("%s_%s[%s_nshifts]/D", wgt_var.c_str(),name.c_str(), name.c_str()) );
      }
    }

    // initialize weight variables -- some won't ever be set
    for( int i = 0; i < knShifts; ++i ) {
      fNwgt[i] = 7; // CAFAna assumes -3,-2,-1,0,1,2,3
      fCvWgts[i] = fisWgtSyst[i] ? 1. : 0; // Non weight systs default to 0.
      for( int j = 0; j < kmaxRwgts; ++j ) {
        fWgts[i][j] = fisWgtSyst[i] ? 1. : 0;
      }
    }
  }

  //------------------------------------------------------------------------------
  void CAFMaker::beginSubRun(const art::SubRun& sr)
  {
    art::Handle< sumdata::POTSummary > pots;
    if( sr.getByLabel("generator",pots) ) meta_pot = pots->totpot;
    else meta_pot = -1.;
  }

  //------------------------------------------------------------------------------
  void CAFMaker::analyze(art::Event const & evt)
  {
    CAFMaker::resetCAFvars();

    art::Handle<std::vector<cvn::Result>> cvnin;
    evt.getByLabel(fCVNLabel, "cvnresultSel", cvnin);

    art::Handle<std::vector<cvn::RegCVNResult>> regcvnin;
    evt.getByLabel(fRegCVNLabel, "regcvnresult", regcvnin);

    art::Handle<pandselect::PandSelectParams> pandSelectParams;;
    evt.getByLabel(fPandSelectParamsLabel, pandSelectParams);

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    fRun = evt.id().run();
    fSubrun = evt.id().subRun();
    fEvent = evt.id().event();

    if ( !pandSelectParams.failedToGet() )
    {
      TVector3 pe( pandSelectParams->energyRecoNue.fLepLorentzVector.Px(), pandSelectParams->energyRecoNue.fLepLorentzVector.Py(), pandSelectParams->energyRecoNue.fLepLorentzVector.Pz() );
      fErecoNue = pandSelectParams->energyRecoNue.fNuLorentzVector.E();
      fRecoLepEnNue = pandSelectParams->energyRecoNue.fLepLorentzVector.E();
      fRecoHadEnNue = pandSelectParams->energyRecoNue.fHadLorentzVector.E();
      fRecoMethodNue = pandSelectParams->energyRecoNue.recoMethodUsed;
      fRecoLepAngNue     = acos( pe.z() / pe.Mag() );

      TVector3 pmu( pandSelectParams->energyRecoNumu.fLepLorentzVector.Px(), pandSelectParams->energyRecoNumu.fLepLorentzVector.Py(), pandSelectParams->energyRecoNumu.fLepLorentzVector.Pz() );
      fErecoNumu = pandSelectParams->energyRecoNumu.fNuLorentzVector.E();
      fRecoLepEnNumu = pandSelectParams->energyRecoNumu.fLepLorentzVector.E();
      fRecoHadEnNumu = pandSelectParams->energyRecoNumu.fHadLorentzVector.E();
      fRecoMethodNumu = pandSelectParams->energyRecoNumu.recoMethodUsed;
      fLongestTrackContNumu = pandSelectParams->energyRecoNumu.longestTrackContained;
      fTrackMomMethodNumu = pandSelectParams->energyRecoNumu.trackMomMethod;
      fRecoLepAngNumu    = acos( pmu.z() / pmu.Mag() );

      std::cout << "pandSelectParams->selShowerJamPandrizzleScore: " << pandSelectParams->selShowerJamPandrizzleScore << std::endl;


      fSelShowerPandrizzleScore = pandSelectParams->selShowerPandrizzleScore;
      fSelShowerJamPandrizzleScore = ((std::fabs(fNueJamPandrizzleCut - 1.0) < std::numeric_limits<float>::epsilon()) && (std::fabs(fAnueJamPandrizzleCut - 1.0) < std::numeric_limits<float>::epsilon())) ? -999.0 : pandSelectParams->selShowerJamPandrizzleScore;
      fSelTrackPandizzleScore = pandSelectParams->selTrackPandizzleScore;

      std::cout << "fErecoNue: " << fErecoNue << std::endl;
      std::cout << "fSelShowerPandrizzleScore: " << fSelShowerPandrizzleScore << std::endl;
      std::cout << "fSelShowerJamPandrizzleScore: " << fSelShowerJamPandrizzleScore << std::endl;
      std::cout << "fErecoNumu: " << fErecoNumu << std::endl;
      std::cout << "fSelTrackPandizzleScore: " << fSelTrackPandizzleScore << std::endl;
    }

    if( !cvnin.failedToGet() ) {
      //using i = cvn::Interaction;
      //if(cvnin->empty() || (*cvnin)[0].fOutput.size() <= i::kNutauOther){
      if(cvnin->empty()){
        fCVNResultIsAntineutrino = fCVNResultNue = fCVNResultNumu = fCVNResultNutau = fCVNResultNC = \
        fCVNResult0Protons = fCVNResult1Protons = fCVNResult2Protons = fCVNResultNProtons = \
        fCVNResult0Pions = fCVNResult1Pions = fCVNResult2Pions = fCVNResultNPions = \
        fCVNResult0Pizeros = fCVNResult1Pizeros = fCVNResult2Pizeros = fCVNResultNPizeros = \
        fCVNResult0Neutrons = fCVNResult1Neutrons = fCVNResult2Neutrons = fCVNResultNNeutrons = -3;
      }
      else{
        //const std::vector<float>& v = (*cvnin)[0].fOutput;
        //fCVNResultNue = v[i::kNueQE] + v[i::kNueRes] + v[i::kNueDIS] + v[i::kNueOther];
        //fCVNResultNumu = v[i::kNumuQE] + v[i::kNumuRes] + v[i::kNumuDIS] + v[i::kNumuOther];
        //fCVNResultNutau = v[i::kNutauQE] + v[i::kNutauRes] + v[i::kNutauDIS] + v[i::kNutauOther]

        fCVNResultIsAntineutrino = (*cvnin)[0].GetIsAntineutrinoProbability();

        fCVNResultNue = (*cvnin)[0].GetNueProbability();
        fCVNResultNumu = (*cvnin)[0].GetNumuProbability();
        fCVNResultNutau = (*cvnin)[0].GetNutauProbability();
        fCVNResultNC = (*cvnin)[0].GetNCProbability();

        fCVNResult0Protons = (*cvnin)[0].Get0protonsProbability();
        fCVNResult1Protons = (*cvnin)[0].Get1protonsProbability();
        fCVNResult2Protons = (*cvnin)[0].Get2protonsProbability();
        fCVNResultNProtons = (*cvnin)[0].GetNprotonsProbability();

        fCVNResult0Pions = (*cvnin)[0].Get0pionsProbability();
        fCVNResult1Pions = (*cvnin)[0].Get1pionsProbability();
        fCVNResult2Pions = (*cvnin)[0].Get2pionsProbability();
        fCVNResultNPions = (*cvnin)[0].GetNpionsProbability();

        fCVNResult0Pizeros = (*cvnin)[0].Get0pizerosProbability();
        fCVNResult1Pizeros = (*cvnin)[0].Get1pizerosProbability();
        fCVNResult2Pizeros = (*cvnin)[0].Get2pizerosProbability();
        fCVNResultNPizeros = (*cvnin)[0].GetNpizerosProbability();

        fCVNResult0Neutrons = (*cvnin)[0].Get0neutronsProbability();
        fCVNResult1Neutrons = (*cvnin)[0].Get1neutronsProbability();
        fCVNResult2Neutrons = (*cvnin)[0].Get2neutronsProbability();
        fCVNResultNNeutrons = (*cvnin)[0].GetNneutronsProbability();

        //std::cout << "fCVNResultNue: " << fCVNResultNue << std::endl;
        //std::cout << "fCVNResultNumu: " << fCVNResultNumu << std::endl;
        //std::cout << "fCVNResultNutau: " << fCVNResultNutau << std::endl;
        //std::cout << "fCVNResultNC: " << fCVNResultNC << std::endl;
      }
    }

    fRegCVNNueE = -1.;  // initializing
    if(!regcvnin.failedToGet()){
      if (!regcvnin->empty()){
        const std::vector<float>& v = (*regcvnin)[0].fOutput;
        fRegCVNNueE = v[0];
      }
    }
  
    art::Handle< std::vector<simb::MCTruth> > mct;
    std::vector< art::Ptr<simb::MCTruth> > truth;
    if( evt.getByLabel("generator", mct) )
      art::fill_ptr_vector(truth, mct);
    else
      mf::LogWarning("CAFMaker") << "No MCTruth.";


    art::Handle< std::vector<simb::MCFlux> > mcf;
    std::vector< art::Ptr<simb::MCFlux> > flux;
    if( evt.getByLabel("generator", mcf) )
      art::fill_ptr_vector(flux, mcf);
    else
      mf::LogWarning("CAFMaker") << "No MCFlux.";

/*
    art::Handle< std::vector<simb::GTruth> > gt;
    std::vector< art::Ptr<simb::GTruth> > gtru;
    if( evt.getByLabel("generator", gt) )
      art::fill_ptr_vector(gtru, gt);
    else
      mf::LogWarning("CAFMaker") << "No GTruth.";
*/

    std::map<int,int> primary_pdg; // track ID to PDG code of primary particle
    std::map<int,int> tid_to_mother; // track ID to mother track ID
    for(size_t i=0; i<truth.size(); i++){

      if(i>1){
        mf::LogWarning("CAFMaker") << "Skipping MC truth index " << i;
        continue;
      }

      fIsFD     = 1; // always FD
      fIsFHC    = 999; // don't know how to get this?
      fIsCC     = !(truth[i]->GetNeutrino().CCNC());  // ccnc is 0=CC 1=NC
      fNuPDG    = truth[i]->GetNeutrino().Nu().PdgCode();
      fNuPDGunosc = flux[i]->fntype;
      fMode     = truth[i]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production; this is different than mode in ND
      fEv       = truth[i]->GetNeutrino().Nu().E();
      fQ2       = truth[i]->GetNeutrino().QSqr();
      fW        = truth[i]->GetNeutrino().W();
      fX        = truth[i]->GetNeutrino().X();
      fY        = truth[i]->GetNeutrino().Y();
      fNuMomX   = truth[i]->GetNeutrino().Nu().Momentum().X();
      fNuMomY   = truth[i]->GetNeutrino().Nu().Momentum().Y();
      fNuMomZ   = truth[i]->GetNeutrino().Nu().Momentum().Z();

      fNuVtxX     = truth[i]->GetNeutrino().Nu().Vx();
      fNuVtxY     = truth[i]->GetNeutrino().Nu().Vy();
      fNuVtxZ     = truth[i]->GetNeutrino().Nu().Vz();

      //Lepton stuff
      fLepPDG     = truth[i]->GetNeutrino().Lepton().PdgCode();
      fLepMomX    = truth[i]->GetNeutrino().Lepton().Momentum().X();
      fLepMomY    = truth[i]->GetNeutrino().Lepton().Momentum().Y();
      fLepMomZ    = truth[i]->GetNeutrino().Lepton().Momentum().Z();
      fLepE       = truth[i]->GetNeutrino().Lepton().Momentum().T();
      fLepNuAngle = truth[i]->GetNeutrino().Nu().Momentum().Vect().Angle(truth[i]->GetNeutrino().Lepton().Momentum().Vect());

      // TODO
      fWeight = 0.;
      fOscPro = 0.;
      //      fOscPro   = fMVAAlg.OscPro(fCCNC,fBeamPdg,fNuPDG,fEtrue);

      nP        = 0;
      nN        = 0;
      nPip      = 0;
      nPim      = 0;
      nPi0      = 0;
      nKp       = 0;
      nKm       = 0;
      nK0       = 0;
      nEM       = 0;
      nOtherHad = 0;
      nNucleus  = 0;
      nUNKNOWN  = 0;

      eP = 0.;
      eN = 0.;
      ePip = 0.;
      ePim = 0.;
      ePi0 = 0.;
      eOther = 0.;

      for( int p = 0; p < truth[i]->NParticles(); p++ ) {
        if( truth[i]->GetParticle(p).StatusCode() == genie::kIStStableFinalState ) {

          int pdg = truth[i]->GetParticle(p).PdgCode();
          double ke = truth[i]->GetParticle(p).E() - truth[i]->GetParticle(p).Mass();

          if     ( pdg == genie::kPdgProton ) {
            nP++;
            eP += ke;
          } else if( pdg == genie::kPdgNeutron ) {
            nN++;
            eN += ke;
          } else if( pdg == genie::kPdgPiP ) {
            nPip++;
            ePip += ke;
          } else if( pdg == genie::kPdgPiM ) {
            nPim++;
            ePim += ke;
          } else if( pdg == genie::kPdgPi0 ) {
            nPi0++;
            ePi0 += ke;
          } else if( pdg == genie::kPdgKP ) {
            nKp++;
            eOther += ke;
          } else if( pdg == genie::kPdgKM ) {
            nKm++;
            eOther += ke;
          } else if( pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0 || pdg == genie::kPdgK0L || pdg == genie::kPdgK0S ) {
            nK0++;
            eOther += ke;
          } else if( pdg == genie::kPdgGamma ) {
            nEM++;
            eOther += ke;
          } else if( genie::pdg::IsHadron(pdg) ) {
            nOtherHad++; // charm mesons, strange and charm baryons, antibaryons, etc.
            eOther += ke;
          } else if( genie::pdg::IsIon(pdg) ) {
            nNucleus++;
          } else {
            nUNKNOWN++;
          }
        }
      }

      /* OLD CAFMAKER
      for( int p = 0; p < truth[i]->NParticles(); p++ ) {
        if( truth[i]->GetParticle(p).StatusCode() == genie::kIStHadronInTheNucleus ) {

          int pdg = truth[i]->GetParticle(p).PdgCode();
          if     ( pdg == genie::kPdgProton ) nP++;
          else if( pdg == genie::kPdgNeutron ) nN++;
          else if( pdg == genie::kPdgPiP ) nPip++;
          else if( pdg == genie::kPdgPiM ) nPim++;
          else if( pdg == genie::kPdgPi0 ) nPi0++;
          else if( pdg == genie::kPdgKP ) nKp++;
          else if( pdg == genie::kPdgKM ) nKm++;
          else if( pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0 || pdg == genie::kPdgK0L || pdg == genie::kPdgK0S ) nK0++;
          else if( pdg == genie::kPdgGamma ) nEM++;
          else if( genie::pdg::IsHadron(pdg) ) nOtherHad++; // charm mesons, strange and charm baryons, antibaryons, etc.
          else if( genie::pdg::IsIon(pdg) ) nNucleus++;
          else nUNKNOWN++;

        }
      }
      */

      // Reweighting variables
      //systtools::ScrubUnityEventResponses(er);

      // struct ParamResponses { 
      //   paramId_t pid;
      //   std::vector<double> responses;
      // }
      // typedef std::vector<ParamResponses> event_unit_response_t;
      // typedef std::vector<event_unit_response_t> EventResponse;

      for( auto &sp : fSystProviders ) {
        std::unique_ptr<systtools::EventAndCVResponse> syst_resp = sp->GetEventVariationAndCVResponse(evt);
        if( !syst_resp ) {
          std::cout << "[ERROR]: Got nullptr systtools::EventResponse from provider "
                    << sp->GetFullyQualifiedName();
          continue;
        }

        for( systtools::EventAndCVResponse::const_iterator itResp = syst_resp->begin(); itResp != syst_resp->end(); ++itResp ) {
          systtools::event_unit_response_w_cv_t const &resp = *itResp;
          for( systtools::event_unit_response_w_cv_t::const_iterator it = resp.begin(); it != resp.end(); ++it ) {
            fNwgt[it->pid] = it->responses.size();
            fCvWgts[it->pid] = it->CV_response;
            std::copy_n(it->responses.begin(), it->responses.size(), fWgts[it->pid]);
          }
        }
      }

      const std::vector< const simb::MCParticle* > parts = pi_serv->MCTruthToParticles_Ps(truth[i]);
      for( size_t pp = 0; pp < parts.size(); ++pp ) {
        int tid = parts[pp]->TrackId();
        int mother = parts[pp]->Mother();
        tid_to_mother.emplace(tid, mother);
        if( mother == 0 ) primary_pdg.emplace(tid, parts[pp]->PdgCode());
        int primaryTid = tid;
        while( mother != 0 ) {
          primaryTid = mother;
          mother = tid_to_mother[primaryTid]; // find primary
        }
        if( primary_pdg.find(primaryTid) == primary_pdg.end() ) std::cout << "Something is wrong" << std::endl;
        else {
          primary_pdg.emplace(tid, primary_pdg[primaryTid]);
        }
      }
    } // loop through MC truth i

    // truth-matching hadronic energy to particles
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);

    eRecoP = 0.;
    eRecoN = 0.;
    eRecoPip = 0.;
    eRecoPim = 0.;
    eRecoPi0 = 0.;
    eRecoOther = 0.;

    eDepP = 0.;
    eDepN = 0.;
    eDepPip = 0.;
    eDepPim = 0.;
    eDepPi0 = 0.;
    eDepOther = 0.;
    
    // Need t0 for electron lifetime correction
    auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double t0 = detprop->TriggerOffset();

    std::map<int,double> tid_charge;
    std::map<int,double> tid_eDep;
    double total_charge = 0.;
    for( size_t i = 0; i < hitlist.size(); ++i ) 
    {
    	art::Ptr<recob::Hit> hit = hitlist[i];
        double charge = hit->Integral();

        double charge_eLifetimeCorrected = hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);

      double deposited_energy = 0;
      if (hit->WireID().Plane == 2) deposited_energy = fCaloAlg.ElectronsFromADCArea( charge_eLifetimeCorrected , 2) * (1.0 / fRecombFactor) / util::kGeVToElectrons;

      std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
      double tote = 0.;
      for(size_t e = 0; e < TrackIDs.size(); ++e) tote += TrackIDs[e].energy;
      for(size_t e = 0; e < TrackIDs.size(); ++e) {
        int primpdg = primary_pdg[TrackIDs[e].trackID];

        if( abs(primpdg) != 11 && abs(primpdg) != 13 && abs(primpdg) != 15 ) 
        {
          tid_charge[TrackIDs[e].trackID] += charge*TrackIDs[e].energy/tote;
          tid_eDep[TrackIDs[e].trackID] += deposited_energy*TrackIDs[e].energy/tote;
          total_charge += charge*TrackIDs[e].energy/tote;
        }
      }
    }

    // choose the hadronic energy for the best reco

    double nuePandrizzleCut = (fUseFHCCut ? fNuePandrizzleCut : fAnuePandrizzleCut);
    double nueJamPandrizzleCut = (fUseFHCCut ? fNueJamPandrizzleCut : fAnueJamPandrizzleCut);
    double nuePandizzleCut = (fUseFHCCut ? fNuePandizzleCut : fAnuePandizzleCut);
    //double numuPandizzleCut = (fUseFHCCut ? fNumuPandizzleCut : fAnumuPandizzleCut);

    const bool passNueSelection(((fSelShowerJamPandrizzleScore > nueJamPandrizzleCut) || (fSelShowerPandrizzleScore > nuePandrizzleCut)) && (fSelTrackPandizzleScore < nuePandizzleCut));
    //const bool passNumuSelection(passNueSelection > false : (fSelTrackPandizzleScore > numuPandizzleCut));

    double ehad = ( passNueSelection ? fRecoHadEnNue : fRecoHadEnNumu );
    for( std::map<int,double>::iterator itTid = tid_charge.begin(); itTid != tid_charge.end(); ++itTid ) {
      int tid = (*itTid).first;
      double energy = (*itTid).second * (ehad/total_charge);

      if( tid < 0 ) tid *= -1;

      int primpdg = primary_pdg[tid];
      if( primpdg == genie::kPdgProton ) eRecoP += energy;
      else if( primpdg == genie::kPdgNeutron ) eRecoN += energy;
      else if( primpdg == genie::kPdgPiP ) eRecoPip += energy;
      else if( primpdg == genie::kPdgPiM ) eRecoPim += energy;
      else if( primpdg == genie::kPdgPi0 ) eRecoPi0 += energy;
      else eRecoOther += energy;
    }

    // Deposited energy
    for( std::map<int,double>::iterator itTid = tid_eDep.begin(); itTid != tid_eDep.end(); ++itTid ) {
      int tid = (*itTid).first;
      double energy_deposited = (*itTid).second;

      if( tid < 0 ) tid *= -1;

      int primpdg = primary_pdg[tid];
      if( primpdg == genie::kPdgProton )       eDepP += energy_deposited;
      else if( primpdg == genie::kPdgNeutron ) eDepN += energy_deposited;
      else if( primpdg == genie::kPdgPiP )     eDepPip += energy_deposited;
      else if( primpdg == genie::kPdgPiM )     eDepPim += energy_deposited;
      else if( primpdg == genie::kPdgPi0 )     eDepPi0 += energy_deposited;
      else eRecoOther += energy_deposited;
    }

    // Reco stuff
    art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
    if (evt.getByLabel("pandoraSel", pfparticleListHandle)){

        std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
        art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

        std::vector<art::Ptr<recob::PFParticle> > nu_pfps;

        for (const art::Ptr<recob::PFParticle> particle : pfparticleList){
            if ((std::fabs(particle->PdgCode()) == 12) || (std::fabs(particle->PdgCode()) == 14) || (std::fabs(particle->PdgCode()) == 16))
                nu_pfps.push_back(particle);
        }

        //Loop over the neutrinos - demand only one
         if (nu_pfps.size() == 1){

             art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

             art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, "pandoraSel");
             const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());

             if (sel_pfp_vertices.size() > 0){ //Nothing to do

                 //always take the first vertex, even if there's more than one
                 art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
                 fRecoNuVtxX = matched_vertex->position().X();
                 fRecoNuVtxY = matched_vertex->position().Y();
                 fRecoNuVtxZ = matched_vertex->position().Z();
             }
         }
    }

    fTree->Fill();
    return;
  }

  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  void CAFMaker::endSubRun(const art::SubRun& sr)
  {
    runPOTTreeVariable = fRun;
    subrunPOTTreeVariable = fSubrun;

    art::Handle< sumdata::POTSummary > potListHandle;

    if (sr.getByLabel("generator", potListHandle))
    {
        std::cout << "YO YO YO" << std::endl;
        potPOTTreeVariable = potListHandle->totpot;
    }
    else
    {
        potPOTTreeVariable = 0.;
    }

    if (fPOT) fPOT->Fill();
  }

  DEFINE_ART_MODULE(CAFMaker)

} // namespace dunemva

#endif // CAFMaker_H

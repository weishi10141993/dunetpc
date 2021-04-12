////////////////////////////////////////////////////////////////////////
// Class:       TestAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        TestAnalysis_module.cc
//
// Generated at Mon Feb  8 15:19:16 2021 by Isobel Mawby,,,Isobel.Mawby@warwick.ac.uk, using cetskelgen
// from cetlib version v3_10_00.
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

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>

#include <vector>
#include <string>

namespace test {
  class TestAnalysis;
}


class test::TestAnalysis : public art::EDAnalyzer {
public:
  explicit TestAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestAnalysis(TestAnalysis const&) = delete;
  TestAnalysis(TestAnalysis&&) = delete;
  TestAnalysis& operator=(TestAnalysis const&) = delete;
  TestAnalysis& operator=(TestAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::string fMCParticleLabel;
  std::string fPFParticleLabel;
  std::string fPFParticleAssociationLabel;

  TTree * fTree;

  std::vector<int> *fTrackID_MC = new std::vector<int>;
  std::vector<int> *fPDG_MC = new std::vector<int>;
  std::vector<float> *fEnergy_MC = new std::vector<float>;

  std::vector<int> *fPfoID = new std::vector<int>;
  std::vector<int> *fHitsU_PFO = new std::vector<int>;
  std::vector<int> *fHitsV_PFO = new std::vector<int>;
  std::vector<int> *fHitsW_PFO = new std::vector<int>;
  std::vector<int> *fHits3D_PFO = new std::vector<int>;
  std::vector<float> *fT0_PFO = new std::vector<float>;
};


test::TestAnalysis::TestAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fMCParticleLabel = p.get<std::string>("MCParticleLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fPFParticleAssociationLabel = p.get<std::string>("PFParticleAssociationLabel");
}

void test::TestAnalysis::analyze(art::Event const& e)
{
  // Get MCParticles 
  const std::vector< art::Ptr<simb::MCParticle> > mcParticleVector = dune_ana::DUNEAnaEventUtils::GetMCParticles(e,fMCParticleLabel);

  if (mcParticleVector.empty())
  {
    std::cout << "ISOBEL: NO MC PARTICLES" << std::endl;
  }

  std::cout << "mcParticleVector.size(): " << mcParticleVector.size() << std::endl;

  // Get PFParticles
  const std::vector< art::Ptr<recob::PFParticle> > pfParticleVector = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);


  std::cout << "pfParticleVector.size(): " << pfParticleVector.size() << std::endl;

  // Get associations between PFParticles and best matched MCParticles
  art::FindManyP<recob::PFParticle> pfParticleAssociations(mcParticleVector, e, fPFParticleAssociationLabel);

  for (const art::Ptr<simb::MCParticle> pMCParticle : mcParticleVector)
  {
    std::vector< art::Ptr<recob::PFParticle> > pPFParticleAssociationsList = pfParticleAssociations.at(pMCParticle.key());

    if (pPFParticleAssociationsList.empty())
      continue;

    if (pPFParticleAssociationsList.size() > 1)
      std::cout << "ISOBEL THIS IS NOT MEANT TO HAPPEN BUDDY" << std::endl;

    art::Ptr<recob::PFParticle> pPFParticle(pPFParticleAssociationsList.front());

    // MC INFO
    ///////////////////////////////////
    int trackID = pMCParticle->TrackId();
    float energy = pMCParticle->E();
    int pdg = pMCParticle->PdgCode();

    fTrackID_MC->push_back(trackID);
    fEnergy_MC->push_back(energy);
    fPDG_MC->push_back(pdg);
    ///////////////////////////////////

    // RECO INFO
    ///////////////////////////////////
    int pfoID = pPFParticle->Self();

    std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 0);
    std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 1);
    std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 2);
    //std::vector<art::Ptr<recob::Hit> > pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pPFParticle, e, fPFParticleLabel); 
    std::vector<art::Ptr<recob::SpacePoint> > pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pPFParticle, e, fPFParticleLabel); 

    int hitsU = pfpHitsView0.size();
    int hitsV = pfpHitsView1.size();
    int hitsW = pfpHitsView2.size();
    int hits3D = pfpHits.size();

    std::vector<art::Ptr<anab::T0>> t0Vector(dune_ana::DUNEAnaPFParticleUtils::GetT0(pPFParticle, e, "pandoraJam"));

    if (t0Vector.size() > 1)
      std::cout << "ISOBEL T0 BIGGER THAN ONE" << std::endl;

    float t0(0.f);
    if (!t0Vector.empty())
      t0 = t0Vector.front()->Time();

    fPfoID->push_back(pfoID);
    fHitsU_PFO->push_back(hitsU);
    fHitsV_PFO->push_back(hitsV);
    fHitsW_PFO->push_back(hitsW);
    fHits3D_PFO->push_back(hits3D);
    fT0_PFO->push_back(t0);
    ///////////////////////////////////
  }

    // Fill Tree
    fTree->Fill();
}

void test::TestAnalysis::beginJob()
{
  fTrackID_MC->clear();
  fPDG_MC->clear();
  fEnergy_MC->clear();

  fPfoID->clear();
  fHitsU_PFO->clear();
  fHitsV_PFO->clear();
  fHitsW_PFO->clear();
  fHits3D_PFO->clear();
  fT0_PFO->clear();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "TestAnalysisTree");
  fTree->Branch("fTrackID_MC", &fTrackID_MC);
  fTree->Branch("fPDG_MC", &fPDG_MC);
  fTree->Branch("fEnergy_MC", &fEnergy_MC);
  fTree->Branch("fPfoID", &fPfoID);
  fTree->Branch("fHitsU_PFO", &fHitsU_PFO);
  fTree->Branch("fHitsV_PFO", &fHitsV_PFO);
  fTree->Branch("fHitsW_PFO", &fHitsW_PFO);
  fTree->Branch("fHits3D_PFO", &fHits3D_PFO);
  fTree->Branch("fT0_PFO", &fT0_PFO);
}

void test::TestAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::TestAnalysis)


////////////////////////////////////////////////////////////////////////
// Class:       DeltaRayVertexAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        DeltaRayVertexAnalysis_module.cc
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

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>
#include <TH1F.h>

#include <vector>
#include <string>

namespace test {
  class DeltaRayVertexAnalysis;
}


class test::DeltaRayVertexAnalysis : public art::EDAnalyzer {
public:
  explicit DeltaRayVertexAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DeltaRayVertexAnalysis(DeltaRayVertexAnalysis const&) = delete;
  DeltaRayVertexAnalysis(DeltaRayVertexAnalysis&&) = delete;
  DeltaRayVertexAnalysis& operator=(DeltaRayVertexAnalysis const&) = delete;
  DeltaRayVertexAnalysis& operator=(DeltaRayVertexAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::string fMCParticleLabel;
  std::string fPFParticleLabel;
  std::string fPFParticleAssociationLabel;

  bool fApplySpaceChargeCorrection;

  TTree * fTree;

  std::vector<int> *fID_DR = new std::vector<int>;

  std::vector<float> *fMCVertexX_DR = new std::vector<float>;
  std::vector<float> *fMCVertexY_DR = new std::vector<float>;
  std::vector<float> *fMCVertexZ_DR = new std::vector<float>;

  std::vector<float> *fMCShiftedVertexX_DR = new std::vector<float>;
  std::vector<float> *fMCShiftedVertexY_DR = new std::vector<float>;
  std::vector<float> *fMCShiftedVertexZ_DR = new std::vector<float>;

  std::vector<float> *fRecoVertexX_DR = new std::vector<float>;
  std::vector<float> *fRecoVertexY_DR = new std::vector<float>;
  std::vector<float> *fRecoVertexZ_DR = new std::vector<float>;

  std::vector<float> *fRecoVertexR_DR = new std::vector<float>;

  std::vector<int> *fIsOneView_DR = new std::vector<int>;
  std::vector<int> *fN3DHits_DR = new std::vector<int>;

  TH1F *fDeltaVertexX_DR;
  TH1F *fDeltaVertexY_DR;
  TH1F *fDeltaVertexZ_DR;
  TH1F *fDeltaVertexR_DR;
};

test::DeltaRayVertexAnalysis::DeltaRayVertexAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fApplySpaceChargeCorrection = p.get<bool>("ApplySpaceChargeCorrection");
  fMCParticleLabel = p.get<std::string>("MCParticleLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fPFParticleAssociationLabel = p.get<std::string>("PFParticleAssociationLabel");
}

void test::DeltaRayVertexAnalysis::analyze(art::Event const& e)
{
  std::cout << "DeltaRayVertexAnalysis - ApplySpaceChargeCorrection: " << fApplySpaceChargeCorrection << std::endl; 

  // Get MCParticles 
  const std::vector< art::Ptr<simb::MCParticle> > mcParticleVector = dune_ana::DUNEAnaEventUtils::GetMCParticles(e,fMCParticleLabel);

  if (mcParticleVector.empty())
  {
    std::cout << "ISOBEL: NO MC PARTICLES" << std::endl;
  }

  // Get associations between PFParticles and best matched MCParticles
  art::FindManyP<recob::PFParticle> pfParticleAssociations(mcParticleVector, e, fPFParticleAssociationLabel);

  for (const art::Ptr<simb::MCParticle> pMCParticle : mcParticleVector)
  {
    // Move on if not a DR
    if (pMCParticle->Process() != "muIoni")
      continue;

    // MC INFO
    ///////////////////////////////////
    // Find start and end trajectory points
    int firstT(-1), lastT(-1);
    lar_pandora::LArPandoraInput::GetTrueStartAndEndPoints(pMCParticle, firstT, lastT);

    if (firstT < 0 && lastT < 0) 
    {
      firstT = 0;
      lastT = 0;
    }

    // Lookup position and kinematics at start and end points
    const float trueVertexX(pMCParticle->Vx(firstT));
    const float trueVertexY(pMCParticle->Vy(firstT));
    const float trueVertexZ(pMCParticle->Vz(firstT));

    float trueShiftedVertexX(trueVertexX);
    float trueShiftedVertexY(trueVertexY);
    float trueShiftedVertexZ(trueVertexZ);

    // Apply space chare correction to the true vertex
    if (fApplySpaceChargeCorrection)
    {
      art::ServiceHandle<spacecharge::SpaceChargeService> SCEService;
      auto const* SCE = SCEService->provider();
      auto gptVertex = geo::Point_t(trueVertexX, trueVertexY, trueVertexZ);
      auto sceVertexOffset = geo::Point_t(0.f, 0.f, 0.f);

      sceVertexOffset = SCE->GetPosOffsets(gptVertex);

      trueShiftedVertexX -= sceVertexOffset.X();
      trueShiftedVertexY += sceVertexOffset.Y();
      trueShiftedVertexZ += sceVertexOffset.Z();
    }

    // Apply time offset corrections to the true vertex - ignore if vertex is not in detector
    try 
    {
      trueShiftedVertexX +=  lar_pandora::LArPandoraInput::GetTrueX0(e, pMCParticle, firstT);
    }
    catch (...)
    {
      continue;
    }
    ////////////////////////////////////

    // RECO INFO
    ///////////////////////////////////
    std::vector< art::Ptr<recob::PFParticle> > pPFParticleAssociationsList = pfParticleAssociations.at(pMCParticle.key());

    if (pPFParticleAssociationsList.empty())
      continue;

    if (pPFParticleAssociationsList.size() > 1)
      std::cout << "ISOBEL THIS IS NOT MEANT TO HAPPEN BUDDY" << std::endl;

    art::Ptr<recob::PFParticle> pPFParticle(pPFParticleAssociationsList.front());
    art::Ptr<recob::Vertex> recoVertex;

    unsigned int hitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 0).size();
    unsigned int hitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 1).size();
    unsigned int hitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 2).size();
    unsigned int hits3D = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pPFParticle, e, fPFParticleLabel).size();

    int isOneView(0);

    if ((hitsView0 == 0) && (hitsView1 * hitsView2 == 0))
      isOneView = 1;

    if ((hitsView1 == 0) && (hitsView0 * hitsView2 == 0))
      isOneView = 1;

    if ((hitsView2 == 0) && (hitsView0 * hitsView1 == 0))
      isOneView = 1;

    try
    { 
      recoVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pPFParticle, e, fPFParticleLabel);
    }
    catch (...)
    {
      std::cout << "no vertex" << std::endl;
      continue;
    }
    ////////////////////////////////////

    // FILL TREES AND HISTOGRAMS
    ///////////////////////////////////
    fID_DR->push_back(pMCParticle->TrackId());

    fMCVertexX_DR->push_back(trueVertexX);
    fMCVertexY_DR->push_back(trueVertexY);
    fMCVertexZ_DR->push_back(trueVertexZ);

    fMCShiftedVertexX_DR->push_back(trueShiftedVertexX);
    fMCShiftedVertexY_DR->push_back(trueShiftedVertexY);
    fMCShiftedVertexZ_DR->push_back(trueShiftedVertexZ);

    fRecoVertexX_DR->push_back(recoVertex->position().X());
    fRecoVertexY_DR->push_back(recoVertex->position().Y());
    fRecoVertexZ_DR->push_back(recoVertex->position().Z());

    const float deltaVertexX = recoVertex->position().X() - trueShiftedVertexX;
    const float deltaVertexY = recoVertex->position().Y() - trueShiftedVertexY;
    const float deltaVertexZ = recoVertex->position().Z() - trueShiftedVertexZ;
    const float deltaR = std::sqrt((deltaVertexX * deltaVertexX) + (deltaVertexY * deltaVertexY) + (deltaVertexZ * deltaVertexZ));

    fRecoVertexR_DR->push_back(deltaR);
    fIsOneView_DR->push_back(isOneView);
    fN3DHits_DR->push_back(hits3D);

    fDeltaVertexX_DR->Fill(deltaVertexX);
    fDeltaVertexY_DR->Fill(deltaVertexY);
    fDeltaVertexZ_DR->Fill(deltaVertexZ);
    fDeltaVertexR_DR->Fill(deltaR);
    ///////////////////////////////////
  }

    // Fill Tree
    fTree->Fill();
}

void test::DeltaRayVertexAnalysis::beginJob()
{
  fID_DR->clear();

  fMCVertexX_DR->clear();
  fMCVertexY_DR->clear();
  fMCVertexZ_DR->clear();

  fMCShiftedVertexX_DR->clear();
  fMCShiftedVertexY_DR->clear();
  fMCShiftedVertexZ_DR->clear();

  fRecoVertexX_DR->clear();
  fRecoVertexY_DR->clear();
  fRecoVertexZ_DR->clear();

  fRecoVertexR_DR->clear();

  fIsOneView_DR->clear();
  fN3DHits_DR->clear();

  art::ServiceHandle<art::TFileService> tfs;

  fDeltaVertexX_DR = tfs->make<TH1F>("deltaVertexXHist", "DeltaVertexXHist", 50, -100, 100);
  fDeltaVertexY_DR = tfs->make<TH1F>("deltaVertexYHist", "DeltaVertexYHist", 50, -100, 100);
  fDeltaVertexZ_DR = tfs->make<TH1F>("deltaVertexZHist", "DeltaVertexZHist", 50, -100, 100);
  fDeltaVertexR_DR = tfs->make<TH1F>("deltaVertexRHist", "DeltaVertexRHist", 50, -100, 100);

  fTree = tfs->make<TTree>("tree", "DeltaRayVertexAnalysisTree");

  fTree->Branch("ID_DR", &fID_DR);
  fTree->Branch("mcVertexX_DR", &fMCVertexX_DR);
  fTree->Branch("mcVertexY_DR", &fMCVertexY_DR);
  fTree->Branch("mcVertexZ_DR", &fMCVertexZ_DR);
  fTree->Branch("mcShiftedVertexX_DR", &fMCShiftedVertexX_DR);
  fTree->Branch("mcShiftedVertexY_DR", &fMCShiftedVertexY_DR);
  fTree->Branch("mcShiftedVertexZ_DR", &fMCShiftedVertexZ_DR);
  fTree->Branch("recoVertexX_DR", &fRecoVertexX_DR);
  fTree->Branch("recoVertexY_DR", &fRecoVertexY_DR);
  fTree->Branch("recoVertexZ_DR", &fRecoVertexZ_DR);
  fTree->Branch("recoVertexR_DR", &fRecoVertexR_DR);
  fTree->Branch("isOneView", &fIsOneView_DR);
  fTree->Branch("n3DHits", &fN3DHits_DR);
}

void test::DeltaRayVertexAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::DeltaRayVertexAnalysis)

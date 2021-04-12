////////////////////////////////////////////////////////////////////////
// Class:       DeltaRayDrawing
// Plugin Type: analyzer (art v3_05_01)
// File:        DeltaRayDrawing_module.cc
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

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>

#include <vector>
#include <string>

namespace test {
  class DeltaRayDrawing;
}

class test::DeltaRayDrawing : public art::EDAnalyzer {
public:
  explicit DeltaRayDrawing(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DeltaRayDrawing(DeltaRayDrawing const&) = delete;
  DeltaRayDrawing(DeltaRayDrawing&&) = delete;
  DeltaRayDrawing& operator=(DeltaRayDrawing const&) = delete;
  DeltaRayDrawing& operator=(DeltaRayDrawing&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::string fMCParticleLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleAssociationLabel;
  std::string fTreeName;
  bool        fApplySpaceChargeCorrection;
  int         fTrackID;

  TTree * fTree;

  std::vector<float> *fRecoPositionX_CR = new std::vector<float>;
  std::vector<float> *fRecoPositionY_CR = new std::vector<float>;
  std::vector<float> *fRecoPositionZ_CR = new std::vector<float>;

  float fMCDirectionX_CR;
  float fMCDirectionY_CR;
  float fMCDirectionZ_CR;
  float fRecoDirectionX_CR;
  float fRecoDirectionY_CR;
  float fRecoDirectionZ_CR;

  std::vector<float> *fRecoPositionX_DR = new std::vector<float>;
  std::vector<float> *fRecoPositionY_DR = new std::vector<float>;
  std::vector<float> *fRecoPositionZ_DR = new std::vector<float>;

  std::vector<float> *fTruePositionX_DR = new std::vector<float>;
  std::vector<float> *fTruePositionY_DR = new std::vector<float>;
  std::vector<float> *fTruePositionZ_DR = new std::vector<float>;

  float fMCDirectionX_DR;
  float fMCDirectionY_DR;
  float fMCDirectionZ_DR;
  float fRecoDirectionX_DR;
  float fRecoDirectionY_DR;
  float fRecoDirectionZ_DR;

  float fMCVertexX_DR;
  float fMCVertexY_DR;
  float fMCVertexZ_DR;
  float fRecoVertexX_DR;
  float fRecoVertexY_DR;
  float fRecoVertexZ_DR;
  float fRecoShowerStartX_DR;
  float fRecoShowerStartY_DR;
  float fRecoShowerStartZ_DR;
};


test::DeltaRayDrawing::DeltaRayDrawing(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fMCParticleLabel = p.get<std::string>("MCParticleLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fShowerLabel = p.get<std::string>("ShowerLabel");
  fPFParticleAssociationLabel = p.get<std::string>("PFParticleAssociationLabel");
  fTreeName = p.get<std::string>("AnalysisTreeName");
  fApplySpaceChargeCorrection = p.get<bool>("ApplySpaceChargeCorrection");
  fTrackID = p.get<int>("TrackID");
}

void test::DeltaRayDrawing::analyze(art::Event const& e)
{
  /////////////////////////////////
  // Get MCParticles 
  const std::vector< art::Ptr<simb::MCParticle> > mcParticleVector = dune_ana::DUNEAnaEventUtils::GetMCParticles(e,fMCParticleLabel);

  if (mcParticleVector.empty())
  {
    std::cout << "ISOBEL: NO MC PARTICLES" << std::endl;
    throw;
  }

  typedef std::map<int, art::Ptr<simb::MCParticle>> IDToMCParticlesMap;
  IDToMCParticlesMap idToMCParticlesMap;

  for (const art::Ptr<simb::MCParticle> pMCParticle : mcParticleVector)
    idToMCParticlesMap[pMCParticle->TrackId()] = pMCParticle;
  /////////////////////////////////

  /////////////////////////////////
  // Get PFParticles
  const std::vector< art::Ptr<recob::PFParticle> > pfParticleVector = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);

  typedef std::map<int, art::Ptr<recob::PFParticle>> IDToPFParticlesMap;
  IDToPFParticlesMap idToPFParticlesMap;

  for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector)
    idToPFParticlesMap[pPFParticle->Self()] = pPFParticle;
  /////////////////////////////////

  /////////////////////////////////
  // Get associations between PFParticles and best matched MCParticles
  art::FindManyP<recob::PFParticle> pfParticleAssociations(mcParticleVector, e, fPFParticleAssociationLabel);
  /////////////////////////////////

  /////////////////////////////////
  // Find MC Particle to draw
  for (const art::Ptr<simb::MCParticle> pMCParticle : mcParticleVector)
  {
    if (pMCParticle->TrackId() != fTrackID)
      continue;

    ///////////////////////////////////
    // Write in MC info

    // Vertex
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
    const float trueEndpointX(pMCParticle->Vx(lastT));
    const float trueEndpointY(pMCParticle->Vy(lastT));
    const float trueEndpointZ(pMCParticle->Vz(lastT));

    float trueShiftedVertexX(trueVertexX);
    float trueShiftedVertexY(trueVertexY);
    float trueShiftedVertexZ(trueVertexZ);
    float trueShiftedEndpointX(trueEndpointX);
    float trueShiftedEndpointY(trueEndpointY);
    float trueShiftedEndpointZ(trueEndpointZ);

    // Apply space chare correction to the true vertex
    if (fApplySpaceChargeCorrection)
    {
      art::ServiceHandle<spacecharge::SpaceChargeService> SCEService;
      auto const* SCE = SCEService->provider();
      auto gptVertex = geo::Point_t(trueVertexX, trueVertexY, trueVertexZ);
      auto gptEndpoint = geo::Point_t(trueEndpointX, trueEndpointY, trueEndpointZ);
      auto sceVertexOffset = geo::Point_t(0.f, 0.f, 0.f);
      auto sceEndpointOffset = geo::Point_t(0.f, 0.f, 0.f);

      sceVertexOffset = SCE->GetPosOffsets(gptVertex);
      sceEndpointOffset = SCE->GetPosOffsets(gptEndpoint);

      trueShiftedVertexX -= sceVertexOffset.X();
      trueShiftedVertexY += sceVertexOffset.Y();
      trueShiftedVertexZ += sceVertexOffset.Z();

      trueShiftedEndpointX -= sceEndpointOffset.X();
      trueShiftedEndpointY += sceEndpointOffset.Y();
      trueShiftedEndpointZ += sceEndpointOffset.Z();
    }

    // Apply time offset corrections to the true vertex - ignore if vertex is not in detector
    try 
    {
      trueShiftedVertexX +=  lar_pandora::LArPandoraInput::GetTrueX0(e, pMCParticle, firstT);
      trueShiftedEndpointX +=  lar_pandora::LArPandoraInput::GetTrueX0(e, pMCParticle, lastT);
    }
    catch (...)
    {
      continue;
    }

    fMCVertexX_DR = trueShiftedVertexX;
    fMCVertexY_DR = trueShiftedVertexY;
    fMCVertexZ_DR = trueShiftedVertexZ;

    // Direction
    int parentTrackID = pMCParticle->Mother();
    const art::Ptr<simb::MCParticle> pParentMCParticle = idToMCParticlesMap.at(parentTrackID);

    unsigned int closestParentTrajectoryIndex(0);

    const unsigned int nTrajectoryPoints = pMCParticle->NumberTrajectoryPoints();

    // ATTN: This MC data has two trajectory points and the second point is unphysical so always use first point
    for (unsigned int i = 0; i < nTrajectoryPoints; ++i)
    {
      const TLorentzVector &position(pMCParticle->Position(i));

      float x0 =  lar_pandora::LArPandoraInput::GetTrueX0(e, pMCParticle, i);

      fTruePositionX_DR->push_back(position[0] + x0);
      fTruePositionY_DR->push_back(position[1]);
      fTruePositionZ_DR->push_back(position[2]);
    }

    // ATTN: Not unit vectors
    //const TLorentzVector &mcDirection(pMCParticle->Momentum(closestTrajectoryIndex));
    const TVector3 mcDirection(trueShiftedEndpointX - trueShiftedVertexX, trueShiftedEndpointY - trueShiftedVertexY, trueShiftedEndpointZ - trueShiftedVertexZ);
    const float mcDirectionMag = std::sqrt(std::pow(mcDirection.X(), 2) + std::pow(mcDirection.Y(), 2) + std::pow(mcDirection.Z(), 2));

    const TLorentzVector &parentMCDirection(pParentMCParticle->Momentum(closestParentTrajectoryIndex));

    fMCDirectionX_CR = parentMCDirection[0];
    fMCDirectionY_CR = parentMCDirection[1];
    fMCDirectionZ_CR = parentMCDirection[2];

    fMCDirectionX_DR = mcDirection.X() / mcDirectionMag;
    fMCDirectionY_DR = mcDirection.Y() / mcDirectionMag;
    fMCDirectionZ_DR = mcDirection.Z() / mcDirectionMag;
    ////////////////////////////////////

    ///////////////////////////////////
    // RECO INFO
    std::vector< art::Ptr<recob::PFParticle> > pPFParticleAssociationsList = pfParticleAssociations.at(pMCParticle.key());

    if (pPFParticleAssociationsList.empty())
      continue;

    if (pPFParticleAssociationsList.size() > 1)
    {
      std::cout << "ISOBEL THIS IS NOT MEANT TO HAPPEN BUDDY" << std::endl;
      throw;
    }

    art::Ptr<recob::PFParticle> pPFParticle(pPFParticleAssociationsList.front());

    const std::vector<art::Ptr<recob::SpacePoint>> hits3D = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pPFParticle, e, fPFParticleLabel);

    for (const art::Ptr<recob::SpacePoint> spacePoint : hits3D)
    {
      fRecoPositionX_DR->push_back(spacePoint->position().X());
      fRecoPositionY_DR->push_back(spacePoint->position().Y());
      fRecoPositionZ_DR->push_back(spacePoint->position().Z());
    }

    art::Ptr<recob::PFParticle> pParentPFParticle = lar_pandora::LArPandoraHelper::GetParentPFParticle(idToPFParticlesMap, pPFParticle);

    const std::vector<art::Ptr<recob::SpacePoint>> parentHits3D = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pParentPFParticle, e, fPFParticleLabel);

    for (const art::Ptr<recob::SpacePoint> spacePoint : parentHits3D)
    {
      fRecoPositionX_CR->push_back(spacePoint->position().X());
      fRecoPositionY_CR->push_back(spacePoint->position().Y());
      fRecoPositionZ_CR->push_back(spacePoint->position().Z());
    }

    art::Ptr<recob::Vertex> recoVertex;
    try
    { 
      recoVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pPFParticle, e, fPFParticleLabel);
    }
    catch (...)
    {
      std::cout << "no vertex" << std::endl;
      continue;
    }

    fRecoVertexX_DR = recoVertex->position().X();
    fRecoVertexY_DR = recoVertex->position().Y();
    fRecoVertexZ_DR = recoVertex->position().Z();

    // ATTN: If not reconstructed as shower, move on
    TVector3 recoDirection;
    TVector3 recoShowerStart;
    if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pPFParticle, e, fPFParticleLabel, fShowerLabel))
    {
      art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pPFParticle, e, fPFParticleLabel, fShowerLabel);
      recoDirection = shower->Direction();
      recoShowerStart = shower->ShowerStart();
    }
    else
    {
      continue;
    }

    recob::Track::Vector_t recoParentDirection;
    unsigned int nRecoTrajectoryPoints(0);
    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pParentPFParticle, e, fPFParticleLabel, fTrackLabel))
    {
      art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pParentPFParticle, e, fPFParticleLabel, fTrackLabel);

      nRecoTrajectoryPoints = track->NumberTrajectoryPoints();

      float recoClosestDistance(std::numeric_limits<float>::max()); unsigned int closestRecoTrajectoryIndex(0);
      for (unsigned int i = 0; i < nRecoTrajectoryPoints; ++i)
      {
	auto position(track->LocationAtPoint(i));
	const TVector3 displacementVector(recoShowerStart.X() - position.x(), recoShowerStart.Y() - position.y(), recoShowerStart.Z() - position.z());
	const float separation(std::sqrt(std::pow(displacementVector[0], 2) + std::pow(displacementVector[1], 2) + std::pow(displacementVector[2], 2)));

        if (separation < recoClosestDistance)
	{
	  recoClosestDistance = separation;
	  closestRecoTrajectoryIndex = i;
	}
      }

      recoParentDirection = track->DirectionAtPoint(closestRecoTrajectoryIndex);
    }
    else
    {
      continue;
    }

    fRecoDirectionX_CR = recoParentDirection.x();
    fRecoDirectionY_CR = recoParentDirection.y();
    fRecoDirectionZ_CR = recoParentDirection.z();

    fRecoDirectionX_DR = recoDirection[0];
    fRecoDirectionY_DR = recoDirection[1];
    fRecoDirectionZ_DR = recoDirection[2];

    fRecoShowerStartX_DR = recoShowerStart[0]; 
    fRecoShowerStartY_DR = recoShowerStart[1];
    fRecoShowerStartZ_DR = recoShowerStart[2];
    ///////////////////////////////////
  }

    // Fill Tree
    fTree->Fill();
}

void test::DeltaRayDrawing::beginJob()
{
  fTruePositionX_DR->clear();
  fTruePositionY_DR->clear();
  fTruePositionZ_DR->clear();

  fRecoPositionX_CR->clear();
  fRecoPositionY_CR->clear();
  fRecoPositionZ_CR->clear();

  fRecoPositionX_DR->clear();
  fRecoPositionY_DR->clear();
  fRecoPositionZ_DR->clear();


  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(), "DeltaRayDrawingTree");
  fTree->Branch("trackID", &fTrackID);

  fTree->Branch("truePositionX_DR", &fTruePositionX_DR);
  fTree->Branch("truePositionY_DR", &fTruePositionY_DR);
  fTree->Branch("truePositionZ_DR", &fTruePositionZ_DR);

  fTree->Branch("recoPositionX_CR", &fRecoPositionX_CR);
  fTree->Branch("recoPositionY_CR", &fRecoPositionY_CR);
  fTree->Branch("recoPositionZ_CR", &fRecoPositionZ_CR);

  fTree->Branch("mcDirectionX_CR", &fMCDirectionX_CR);
  fTree->Branch("mcDirectionY_CR", &fMCDirectionY_CR);
  fTree->Branch("mcDirectionZ_CR", &fMCDirectionZ_CR);

  fTree->Branch("recoDirectionX_CR", &fRecoDirectionX_CR);
  fTree->Branch("recoDirectionY_CR", &fRecoDirectionY_CR);
  fTree->Branch("recoDirectionZ_CR", &fRecoDirectionZ_CR);

  fTree->Branch("recoPositionX_DR", &fRecoPositionX_DR);
  fTree->Branch("recoPositionY_DR", &fRecoPositionY_DR);
  fTree->Branch("recoPositionZ_DR", &fRecoPositionZ_DR);

  fTree->Branch("mcDirectionX_DR", &fMCDirectionX_DR);
  fTree->Branch("mcDirectionY_DR", &fMCDirectionY_DR);
  fTree->Branch("mcDirectionZ_DR", &fMCDirectionZ_DR);

  fTree->Branch("recoDirectionX_DR", &fRecoDirectionX_DR);
  fTree->Branch("recoDirectionY_DR", &fRecoDirectionY_DR);
  fTree->Branch("recoDirectionZ_DR", &fRecoDirectionZ_DR);

  fTree->Branch("mcVertexX_DR", &fMCVertexX_DR);
  fTree->Branch("mcVertexY_DR", &fMCVertexY_DR);
  fTree->Branch("mcVertexZ_DR", &fMCVertexZ_DR);

  fTree->Branch("recoVertexX_DR", &fRecoVertexX_DR);
  fTree->Branch("recoVertexY_DR", &fRecoVertexY_DR);
  fTree->Branch("recoVertexZ_DR", &fRecoVertexZ_DR);

  fTree->Branch("recoShowerStartX_DR", &fRecoShowerStartX_DR);
  fTree->Branch("recoShowerStartY_DR", &fRecoShowerStartY_DR);
  fTree->Branch("recoShowerStartZ_DR", &fRecoShowerStartZ_DR);
}

void test::DeltaRayDrawing::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::DeltaRayDrawing)


    // This would tell you if you got the parent links correct...
    /*
    std::vector< art::Ptr<recob::PFParticle> > pParentPFParticleAssociationsList = pfParticleAssociations.at(pParentMCParticle.key());
    art::Ptr<recob::PFParticle> pParentPFParticle(pParentPFParticleAssociationsList.front());

    if (pParentPFParticle.empty())
      continue;

    if (pParentPFParticle.size() > 1)
      std::cout << "ISOBEL THIS IS NOT MEANT TO HAPPEN BUDDY" << std::endl;
    */

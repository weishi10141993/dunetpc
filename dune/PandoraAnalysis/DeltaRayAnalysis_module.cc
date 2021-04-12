////////////////////////////////////////////////////////////////////////
// Class:       DeltaRayAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        DeltaRayAnalysis_module.cc
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

#include <vector>
#include <string>

namespace test {
  class DeltaRayAnalysis;
}

class test::DeltaRayAnalysis : public art::EDAnalyzer {
public:
  explicit DeltaRayAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DeltaRayAnalysis(DeltaRayAnalysis const&) = delete;
  DeltaRayAnalysis(DeltaRayAnalysis&&) = delete;
  DeltaRayAnalysis& operator=(DeltaRayAnalysis const&) = delete;
  DeltaRayAnalysis& operator=(DeltaRayAnalysis&&) = delete;

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

  TTree * fTree;

  bool fApplySpaceChargeCorrection;

  std::vector<int> *fID_DR = new std::vector<int>;
  std::vector<int> *fEventNumber = new std::vector<int>;

  std::vector<float> *fMCDirectionX_DR = new std::vector<float>;
  std::vector<float> *fMCDirectionY_DR = new std::vector<float>;
  std::vector<float> *fMCDirectionZ_DR = new std::vector<float>;
  std::vector<float> *fRecoDirectionX_DR = new std::vector<float>;
  std::vector<float> *fRecoDirectionY_DR = new std::vector<float>;
  std::vector<float> *fRecoDirectionZ_DR = new std::vector<float>;

  std::vector<float> *fMCDirectionX_CR = new std::vector<float>;
  std::vector<float> *fMCDirectionY_CR = new std::vector<float>;
  std::vector<float> *fMCDirectionZ_CR = new std::vector<float>;
  std::vector<float> *fRecoDirectionX_CR = new std::vector<float>;
  std::vector<float> *fRecoDirectionY_CR = new std::vector<float>;
  std::vector<float> *fRecoDirectionZ_CR = new std::vector<float>;

  std::vector<float> *fMCOpeningAngle = new std::vector<float>;
  std::vector<float> *fRecoOpeningAngle = new std::vector<float>;

  std::vector<float> *fCosmicRayOpeningAngle = new std::vector<float>;
  std::vector<float> *fDeltaRayOpeningAngle = new std::vector<float>;

  std::vector<int> *fNViews_DR = new std::vector<int>;
  std::vector<int> *fN3DHits_DR = new std::vector<int>;
};


test::DeltaRayAnalysis::DeltaRayAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fApplySpaceChargeCorrection = p.get<bool>("ApplySpaceChargeCorrection");
  fMCParticleLabel = p.get<std::string>("MCParticleLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fShowerLabel = p.get<std::string>("ShowerLabel");
  fPFParticleAssociationLabel = p.get<std::string>("PFParticleAssociationLabel");
  fTreeName = p.get<std::string>("AnalysisTreeName");
}

void test::DeltaRayAnalysis::analyze(art::Event const& e)
{
  std::cout << "DeltaRayDirectionAnalysis - ApplySpaceChargeCorrection: " << fApplySpaceChargeCorrection << std::endl; 

  // Get MCParticles 
  const std::vector< art::Ptr<simb::MCParticle> > mcParticleVector = dune_ana::DUNEAnaEventUtils::GetMCParticles(e,fMCParticleLabel);

  if (mcParticleVector.empty())
  {
    std::cout << "ISOBEL: NO MC PARTICLES" << std::endl;
  }

  typedef std::map<int, art::Ptr<simb::MCParticle>> IDToMCParticlesMap;
  IDToMCParticlesMap idToMCParticlesMap;

  for (const art::Ptr<simb::MCParticle> pMCParticle : mcParticleVector)
    idToMCParticlesMap[pMCParticle->TrackId()] = pMCParticle;

  // Get PFParticles
  const std::vector< art::Ptr<recob::PFParticle> > pfParticleVector = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);

  typedef std::map<int, art::Ptr<recob::PFParticle>> IDToPFParticlesMap;
  IDToPFParticlesMap idToPFParticlesMap;

  for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector)
    idToPFParticlesMap[pPFParticle->Self()] = pPFParticle;

  // Get associations between PFParticles and best matched MCParticles
  art::FindManyP<recob::PFParticle> pfParticleAssociations(mcParticleVector, e, fPFParticleAssociationLabel);

  for (const art::Ptr<simb::MCParticle> pMCParticle : mcParticleVector)
  {
    // Move on if not a DR
    if (pMCParticle->Process() != "muIoni")
      continue;

    // MC INFO
    ///////////////////////////////////
    int parentTrackID = pMCParticle->Mother();
    const art::Ptr<simb::MCParticle> pParentMCParticle = idToMCParticlesMap.at(parentTrackID);

    unsigned int closestParentTrajectoryIndex(0);
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

    // ATTN: Not unit vectors
    //const TLorentzVector &mcDirection(pMCParticle->Momentum(closestTrajectoryIndex));
    const TVector3 mcDirection(trueShiftedEndpointX - trueShiftedVertexX, trueShiftedEndpointY - trueShiftedVertexY, trueShiftedEndpointZ - trueShiftedVertexZ);
    const TLorentzVector &parentMCDirection(pParentMCParticle->Momentum(closestParentTrajectoryIndex));

    // MC opening angle
    const float mcDirectionMag = std::sqrt(std::pow(mcDirection.X(), 2) + std::pow(mcDirection.Y(), 2) + std::pow(mcDirection.Z(), 2));
    const TVector3 mcDirectionUnit = std::fabs(mcDirectionMag) < std::numeric_limits<float>::epsilon() ? TVector3(0.f, 0.f, 0.f) : 
      TVector3(mcDirection.X() / mcDirectionMag, mcDirection.Y() / mcDirectionMag, mcDirection.Z() / mcDirectionMag);
    const float parentMCDirectionMag = std::sqrt(std::pow(parentMCDirection[0], 2) + std::pow(parentMCDirection[1], 2) + std::pow(parentMCDirection[2], 2));
    const TVector3 parentMCDirectionUnit = std::fabs(parentMCDirectionMag) < std::numeric_limits<float>::epsilon() ? TVector3(0.f, 0.f, 0.f) :
      TVector3(parentMCDirection[0] / parentMCDirectionMag, parentMCDirection[1] / parentMCDirectionMag, parentMCDirection[2] / parentMCDirectionMag);
    const float mcDot = (mcDirection.X() * parentMCDirection[0]) + (mcDirection.Y() * parentMCDirection[1]) + (mcDirection.Z() * parentMCDirection[2]);
    const float mcOpeningAngle = ((std::fabs(mcDirectionMag) < std::numeric_limits<float>::epsilon()) || (std::fabs(parentMCDirectionMag) < std::numeric_limits<float>::epsilon())) ? 
      -4.f : std::acos(mcDot / (mcDirectionMag * parentMCDirectionMag));
    ////////////////////////////////////

    // RECO INFO
    ///////////////////////////////////
    std::vector< art::Ptr<recob::PFParticle> > pPFParticleAssociationsList = pfParticleAssociations.at(pMCParticle.key());

    if (pPFParticleAssociationsList.empty())
      continue;

    if (pPFParticleAssociationsList.size() > 1)
    {
      std::cout << "ISOBEL THIS IS NOT MEANT TO HAPPEN BUDDY" << std::endl;
      throw;
    }

    art::Ptr<recob::PFParticle> pPFParticle(pPFParticleAssociationsList.front());

    unsigned int hitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 0).size();
    unsigned int hitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 1).size();
    unsigned int hitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pPFParticle, e, fPFParticleLabel, 2).size();
    unsigned int hits3D = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pPFParticle, e, fPFParticleLabel).size();

    int nViews(0);

    if (hitsView0 != 0)
      ++nViews;

    if (hitsView1 != 0)
      ++nViews;

    if (hitsView2 != 0)
      ++nViews;

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

    art::Ptr<recob::PFParticle> pParentPFParticle = lar_pandora::LArPandoraHelper::GetParentPFParticle(idToPFParticlesMap, pPFParticle);

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

    // Reco opening angle
    const float recoDirectionMag = std::sqrt(std::pow(recoDirection.X(), 2) + std::pow(recoDirection.Y(), 2) + std::pow(recoDirection.Z(), 2));
    const float recoParentDirectionMag = std::sqrt(std::pow(recoParentDirection.x(), 2) + std::pow(recoParentDirection.y(), 2) + std::pow(recoParentDirection.z(), 2));
    const float recoDot = (recoDirection.X() * recoParentDirection.x()) + (recoDirection.Y() * recoParentDirection.y()) + (recoDirection.Z() * recoParentDirection.z());
    const float recoOpeningAngle = ((std::fabs(recoDirectionMag) < std::numeric_limits<float>::epsilon()) || (std::fabs(recoParentDirectionMag) < std::numeric_limits<float>::epsilon())) ? 
      -4.f : std::acos(recoDot / (recoDirectionMag * recoParentDirectionMag));

    const TVector3 recoDirectionUnit = std::fabs(recoDirectionMag) < std::numeric_limits<float>::epsilon() ? TVector3(0.f, 0.f, 0.f) :
         TVector3(recoDirection.X() / recoDirectionMag, recoDirection.Y() / recoDirectionMag, recoDirection.Z() / recoDirectionMag);
    const float deltaRayDot = (mcDirectionUnit.X() * recoDirectionUnit.X()) + (mcDirectionUnit.Y() * recoDirectionUnit.Y()) + (mcDirectionUnit.Z() * recoDirectionUnit.Z());
    const float deltaRayOpeningAngle = ((std::fabs(mcDirectionMag) < std::numeric_limits<float>::epsilon()) || (std::fabs(recoDirectionMag) < std::numeric_limits<float>::epsilon())) ?
      -4.f : std::acos(deltaRayDot);

    const TVector3 recoParentDirectionUnit = std::fabs(recoParentDirectionMag) < std::numeric_limits<float>::epsilon() ? TVector3(0.f, 0.f, 0.f) : 
        TVector3(recoParentDirection.X() / recoParentDirectionMag, recoParentDirection.Y() / recoParentDirectionMag, recoParentDirection.Z() / recoParentDirectionMag);
    const float cosmicRayDot = (parentMCDirectionUnit.X() * recoParentDirectionUnit.X()) + (parentMCDirectionUnit.Y() * recoParentDirectionUnit.Y()) + 
      (parentMCDirectionUnit.Z() * recoParentDirectionUnit.Z());
    const float cosmicRayOpeningAngle = ((std::fabs(parentMCDirectionMag) < std::numeric_limits<float>::epsilon()) || (std::fabs(recoParentDirectionMag) < std::numeric_limits<float>::epsilon())) ?
      -4.f : std::acos(cosmicRayDot);

    /*
    std::cout << "//////////////////////////////////" << std::endl;
    std::cout << "nTrajectoryPoints: " << nTrajectoryPoints << ", nParentTrajectoryPoints: " << nParentTrajectoryPoints << std::endl;
    std::cout << "closestTrajectoryIndex: " << closestTrajectoryIndex << ", closestParentTrajectoryIndex: " << closestParentTrajectoryIndex << std::endl;
    std::cout << "mcDirection: " << mcDirection[0] << ", " << mcDirection[1] << ", " << mcDirection[2] << std::endl;
    std::cout << "parentMCDirection: " << parentMCDirection[0] << ", " << parentMCDirection[1] << ", " << parentMCDirection[2] << std::endl;
    std::cout << "mcDirectionMag: " << mcDirectionMag << std::endl;
    std::cout << "parentMCDirectionMag: " << parentMCDirectionMag << std::endl;
    std::cout << "mcDot: " << mcDot << std::endl;
    std::cout << "mcOpeningAngle: " << mcOpeningAngle << std::endl; 
    std::cout << "recoDirection: " << recoDirection.X() << ", " << recoDirection.Y() << ", " << recoDirection.Z() << std::endl;
    std::cout << "recoShowerStart: " << recoShowerStart.X() << ", " << recoShowerStart.Y() << ", " << recoShowerStart.Z() << std::endl;
    std::cout << "nRecoTrajectoryPoints: " << nRecoTrajectoryPoints << std::endl;
    std::cout << "recoParentDirection: " << recoParentDirection.x() << ", " << recoParentDirection.y() << ", " << recoParentDirection.z() << std::endl;
    std::cout << "recoDirectionMag: " << recoDirectionMag << std::endl;
    std::cout << "recoParentDirectionMag: " << recoParentDirectionMag << std::endl;
    std::cout << "recoDot" << recoDot << std::endl;
    std::cout << "recoOpeningAngle: " << recoOpeningAngle << std::endl; 
    std::cout << "mcDirectionUnit: " << mcDirectionUnit.X() << ", " << mcDirectionUnit.Y() << ", " << mcDirectionUnit.Z() << std::endl;
    std::cout << "parentMCDirectionUnit: " << parentMCDirectionUnit.X() << ", " << parentMCDirectionUnit.Y() << ", " << parentMCDirectionUnit.Z() << std::endl;
    std::cout << "recoDirectionUnit: " << recoDirectionUnit.X() << ", " << recoDirectionUnit.Y() << ", " << recoDirectionUnit.Z() << std::endl;
    std::cout << "recoParentDirectionUnit: " << recoParentDirectionUnit.X() << ", " << recoParentDirectionUnit.Y() << ", " << recoParentDirectionUnit.Z() << std::endl;
    std::cout << "deltaRayDot: " << deltaRayDot << std::endl;
    std::cout << "cosmicRayDot: " << cosmicRayDot << std::endl;
    std::cout << "deltaRayOpeningAngle: " << deltaRayOpeningAngle << std::endl;
    std::cout << "cosmicRayOpeningAngle: " << cosmicRayOpeningAngle << std::endl;
    */
    fEventNumber->push_back(e.id().event());
    fID_DR->push_back(pMCParticle->TrackId());

    fMCDirectionX_CR->push_back(parentMCDirection[0]);
    fMCDirectionY_CR->push_back(parentMCDirection[1]);
    fMCDirectionZ_CR->push_back(parentMCDirection[2]);

    fMCDirectionX_DR->push_back(mcDirectionUnit.X());
    fMCDirectionY_DR->push_back(mcDirectionUnit.Y());
    fMCDirectionZ_DR->push_back(mcDirectionUnit.Z());

    fMCOpeningAngle->push_back(mcOpeningAngle);

    fRecoDirectionX_CR->push_back(recoParentDirection.x());
    fRecoDirectionY_CR->push_back(recoParentDirection.y());
    fRecoDirectionZ_CR->push_back(recoParentDirection.z());

    fRecoDirectionX_DR->push_back(recoDirection.X());
    fRecoDirectionY_DR->push_back(recoDirection.Y());
    fRecoDirectionZ_DR->push_back(recoDirection.Z());

    fRecoOpeningAngle->push_back(recoOpeningAngle);

    fDeltaRayOpeningAngle->push_back(deltaRayOpeningAngle);
    fCosmicRayOpeningAngle->push_back(cosmicRayOpeningAngle);

    fNViews_DR->push_back(nViews);
    fN3DHits_DR->push_back(hits3D);
    ///////////////////////////////////
  }

    // Fill Tree
    fTree->Fill();
}

void test::DeltaRayAnalysis::beginJob()
{
  fEventNumber->clear();
  fID_DR->clear();
  fMCDirectionX_CR->clear();
  fMCDirectionY_CR->clear();
  fMCDirectionZ_CR->clear();
  fMCDirectionX_DR->clear();
  fMCDirectionY_DR->clear();
  fMCDirectionZ_DR->clear();
  fMCOpeningAngle->clear();
  fRecoDirectionX_CR->clear();
  fRecoDirectionY_CR->clear();
  fRecoDirectionZ_CR->clear();
  fRecoDirectionX_DR->clear();
  fRecoDirectionY_DR->clear();
  fRecoDirectionZ_DR->clear();
  fRecoOpeningAngle->clear();
  fDeltaRayOpeningAngle->clear();
  fCosmicRayOpeningAngle->clear();
  fNViews_DR->clear();
  fN3DHits_DR->clear();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(), "DeltaRayAnalysisTree");
  fTree->Branch("ID_DR", &fID_DR);
  fTree->Branch("EventNumber", &fEventNumber);
  fTree->Branch("mcDirectionX_CR", &fMCDirectionX_CR);
  fTree->Branch("mcDirectionY_CR", &fMCDirectionY_CR);
  fTree->Branch("mcDirectionZ_CR", &fMCDirectionZ_CR);
  fTree->Branch("mcDirectionX_DR", &fMCDirectionX_DR);
  fTree->Branch("mcDirectionY_DR", &fMCDirectionY_DR);
  fTree->Branch("mcDirectionZ_DR", &fMCDirectionZ_DR);
  fTree->Branch("mcOpeningAngle", &fMCOpeningAngle);
  fTree->Branch("recoDirectionX_CR", &fRecoDirectionX_CR);
  fTree->Branch("recoDirectionY_CR", &fRecoDirectionY_CR);
  fTree->Branch("recoDirectionZ_CR", &fRecoDirectionZ_CR);
  fTree->Branch("recoDirectionX_DR", &fRecoDirectionX_DR);
  fTree->Branch("recoDirectionY_DR", &fRecoDirectionY_DR);
  fTree->Branch("recoDirectionZ_DR", &fRecoDirectionZ_DR);
  fTree->Branch("recoOpeningAngle", &fRecoOpeningAngle);
  fTree->Branch("deltaRayOpeningAngle", &fDeltaRayOpeningAngle);
  fTree->Branch("cosmicRayOpeningAngle", &fCosmicRayOpeningAngle);
  fTree->Branch("nViews", &fNViews_DR);
  fTree->Branch("n3DHits", &fN3DHits_DR);
}

void test::DeltaRayAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::DeltaRayAnalysis)


    // This would tell you if you got the parent links correct...
    /*
    std::vector< art::Ptr<recob::PFParticle> > pParentPFParticleAssociationsList = pfParticleAssociations.at(pParentMCParticle.key());
    art::Ptr<recob::PFParticle> pParentPFParticle(pParentPFParticleAssociationsList.front());

    if (pParentPFParticle.empty())
      continue;

    if (pParentPFParticle.size() > 1)
      std::cout << "ISOBEL THIS IS NOT MEANT TO HAPPEN BUDDY" << std::endl;
    */

#include "CheatRecoShowerSelector.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

FDSelectionTools::CheatRecoShowerSelector::CheatRecoShowerSelector(fhicl::ParameterSet const& ps) :
    fNuGenModuleLabel(ps.get< std::string> ("ModuleLabels.NuGenModuleLabel")),
    fShowerModuleLabel(ps.get< std::string> ("ModuleLabels.ShowerModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fPandrizzleAlg(ps),
    fRecoNuVtxX(-9999),
    fRecoNuVtxY(-9999),
    fRecoNuVtxZ(-9999)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Shower> FDSelectionTools::CheatRecoShowerSelector::SelectShower(art::Event const & evt)
{
  bool foundSignal(false);
  art::Ptr<recob::Shower> selSignalShower, selShower;

  // Get the showers
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerList;
  if (evt.getByLabel(fShowerModuleLabel,showerListHandle))
    art::fill_ptr_vector(showerList, showerListHandle);

  // Get the PFParticles
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  // Get the generator record
  art::Handle< std::vector<simb::MCTruth> > mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle))
    art::fill_ptr_vector(mcList, mcTruthListHandle);

  if (mcList.size() != 1)
  {
      std::cout << "CheatRecoShowerSelector: There are " << mcList.size() << " MCTruth in this event. No shower selected" << std::endl;
      return selShower;
  }

  const bool isFHC = (mcList[0]->GetNeutrino().Nu().PdgCode() > 0);
  const bool isNC = mcList[0]->GetNeutrino().CCNC();
  const bool isNue = (std::abs(mcList[0]->GetNeutrino().Nu().PdgCode()) == 12);
  const bool isCCNue = (isNue && !isNC);

  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

  // Check we have one neutrino...
  if (nu_pfps.size() != 1)
  {
    std::cout << "CheatRecoShowerSelector: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return selShower;
  }

  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  // Get neutrino vertex for pandrizzle MVA
  art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());

  if (sel_pfp_vertices.size() == 0)
  {
    return selShower;
  }
  else if (sel_pfp_vertices.size() > 1)
  {
    std::cout<< "CheatRecoShowerSelector: FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  //Loop over each PFP, find the associated showers and then find which one is the highest energy
  double highestSignalPandrizzleScore = std::numeric_limits<double>::lowest(), highestPandrizzleScore = std::numeric_limits<double>::lowest();

  for (art::Ptr<recob::PFParticle> pfParticle : pfparticleList)
  {
    //Now get the associated shower 
    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, fShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(pfParticle.key());

    if (pfp_shower_vector.size() > 1)
    { 
      //Found a PFP with more than one shower matched.  Complain and exit
      std::cout<< "Number of associated showers to a PFP is greater than 1: " << pfp_shower_vector.size() << std::endl;
      continue;
    } 
    else if (pfp_shower_vector.size() == 0)
    { 
      //Don't bother complaining of no shower was found.  It's either missing or its a shower.  No biggie here
      continue;
    }

    const art::Ptr<recob::Shower> shower = pfp_shower_vector[0];

    // Get hits to perform truth matching
    art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);
    const std::vector<art::Ptr<recob::Hit> > showerHits = fmhs.at(shower.key());

    // Find true PDG of reco shower
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, 1);

    int matchedPDG = -999;
    bool primary = false;

    if (TruthMatchUtils::Valid(g4id)){
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
        matchedPDG = matched_mcparticle->PdgCode();
        //std::cout << "matchedPDG: " << matchedPDG << std::endl;
        //std::cout << "primary shower process: " << matched_mcparticle->Process() << std::endl;
        primary = (matched_mcparticle->Process() == "primary");
    }

    const bool isPrimaryElectron = ((isFHC && (matchedPDG == 11) && primary) || (!isFHC && (matchedPDG == -11) && primary));

    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(shower, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), evt));
    const double pandrizzleScore = pandrizzleRecord.GetMVAScore();

    if (isPrimaryElectron && (pandrizzleScore > highestSignalPandrizzleScore))
    {
        foundSignal = true;
        highestSignalPandrizzleScore = pandrizzleScore;
        selSignalShower = shower;
    }

    if (pandrizzleScore > highestPandrizzleScore)
    {
        highestPandrizzleScore = pandrizzleScore;
        selShower = shower;
    }
  }

  //std::cout << "found? " << (foundSignal ? "yep" : "no") << std::endl;
  //std::cout << "highestSignalPandrizzleScore: " << highestSignalPandrizzleScore << std::endl;
  //std::cout << "highestPandrizzleScore: " << highestPandrizzleScore << std::endl;

  return (isCCNue ? (foundSignal ? selSignalShower : selShower) : selShower);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "CheatRecoPfoShowerSelector.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

FDSelectionTools::CheatRecoPfoShowerSelector::CheatRecoPfoShowerSelector(fhicl::ParameterSet const& ps) :
    fNuGenModuleLabel(ps.get< std::string > ("ModuleLabels.NuGenModuleLabel")),
    fAllShowerModuleLabel(ps.get< std::string > ("AllShowerModuleLabel")),
    fFittedShowerModuleLabel(ps.get< std::string > ("FittedShowerModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string > ("ModuleLabels.PFParticleModuleLabel")),
    fPandrizzleAlg(ps),
    fRecoNuVtxX(-9999),
    fRecoNuVtxY(-9999),
    fRecoNuVtxZ(-9999)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Shower> FDSelectionTools::CheatRecoPfoShowerSelector::SelectShower(art::Event const & evt)
{
  bool foundSignal(false);
  art::Ptr<recob::Shower> selSignalShower, selShower;

  // For later - to obtain shower hits to use in truth matching
  art::Handle< std::vector<recob::Shower> > showerListHandle_forced;
  std::vector<art::Ptr<recob::Shower> > showerList_forced;
  if (evt.getByLabel(fAllShowerModuleLabel, showerListHandle_forced))
    art::fill_ptr_vector(showerList_forced, showerListHandle_forced);

  // To loop through all pandora pfos
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  // Get the generator record (to find whether event is FHC/RHC)
  art::Handle< std::vector<simb::MCTruth> > mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle))
    art::fill_ptr_vector(mcList, mcTruthListHandle);

  if (mcList.size() != 1)
  {
      std::cout << "CheatRecoPfoShowerSelector: There are " << mcList.size() << " MCTruth in this event. No shower selected" << std::endl;
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
    std::cout << "CheatRecoPfoShowerSelector: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return selShower;
  }

  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  // Get neutrino vertex for pandrizzle MVA
  art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());

  if (sel_pfp_vertices.size() == 0)
  {
    std::cout << "CheatRecoPfoShowerSelector: No reconstructed Nu vertex " << std::endl;
    return selShower;
  }
  else if (sel_pfp_vertices.size() > 1)
  {
    std::cout<< "CheatRecoPfoShowerSelector: FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  //Loop over each PFP, find the associated showers
  double highestSignalPandrizzleScore = std::numeric_limits<double>::lowest(), highestPandrizzleScore = std::numeric_limits<double>::lowest();

  for (art::Ptr<recob::PFParticle> pfParticle : pfparticleList)
  {
    // Get the forced shower fit 
    art::FindManyP<recob::Shower> fmspfp_forced(pfparticleListHandle, evt, fAllShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector_forced = fmspfp_forced.at(pfParticle.key());

    //Now get the associated shower - may not exist
    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, fFittedShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(pfParticle.key());

    if ((pfp_shower_vector_forced.size() > 1) || (pfp_shower_vector.size() > 1))
    { 
      //Found a PFP with more than one shower matched.  Complain and exit
      std::cout<< "Number of associated showers to a PFP is greater than 1" << std::endl;
      continue;
    } 

    if (pfp_shower_vector_forced.size() == 0)
    {
        // Maybe a shower couldn't be fitted? Who knows the ways of larsoft
        continue;
    }

    const bool isFitted(pfp_shower_vector.size() != 0);
    const art::Ptr<recob::Shower> shower_forced = pfp_shower_vector_forced[0];

    // Get hits to perform truth matching
    art::FindManyP<recob::Hit> fmhs_forced(showerListHandle_forced, evt, fAllShowerModuleLabel);
    const std::vector<art::Ptr<recob::Hit> > showerHits_forced = fmhs_forced.at(shower_forced.key());

    // Find true PDG of forced shower
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits_forced, 1);

    int matchedPDG = -999;
    bool primary = false;

    if (TruthMatchUtils::Valid(g4id)){
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
        matchedPDG = matched_mcparticle->PdgCode();
        primary = (matched_mcparticle->Process() == "primary");
        //std::cout << "matchedPDG: " << matchedPDG << std::endl;
        //std::cout << "proces: " << matched_mcparticle->Process() << std::endl;
        //std::cout << "is fitted? " << (isFitted ? "yes" : "no") << std::endl;
    }

    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(shower_forced, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), evt));
    const double pandrizzleScore = pandrizzleRecord.GetMVAScore();
    //std::cout << "pandrizzleScore: " << pandrizzleScore << std::endl;

    const bool isForcedPrimaryElectron = ((isFHC && (matchedPDG == 11) && primary) || (!isFHC && (matchedPDG == -11) && primary));

    if (isForcedPrimaryElectron)
    {
        //std::cout << "found a primary electron" << std::endl;
        //std::cout << "is electron fitted? " << (isFitted ? "yes" : "no") << std::endl;

        if (pandrizzleScore > highestSignalPandrizzleScore)
        {
            foundSignal = true;
            highestSignalPandrizzleScore = pandrizzleScore;
            selSignalShower = shower_forced;
        }
    }
    else if (isFitted)
    {
        if (pandrizzleScore > highestPandrizzleScore)
        {
            highestPandrizzleScore = pandrizzleScore;
            selShower = shower_forced;
        }
    }
  }  

  /*
  std::cout << "found? " << (foundSignal ? "yep" : "no") << std::endl;
  std::cout << "highestSignalPandrizzleScore: " << highestSignalPandrizzleScore << std::endl;
  std::cout << "highestPandrizzleScore: " << highestPandrizzleScore << std::endl;
  */

  return (isCCNue ? (foundSignal ? selSignalShower : selShower) : selShower);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

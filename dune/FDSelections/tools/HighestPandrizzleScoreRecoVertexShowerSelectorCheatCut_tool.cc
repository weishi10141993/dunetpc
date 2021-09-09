#include "HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut.h"

//#include "FDSelectionUtils.h"
#include "dune/FDSelections/FDSelectionUtils.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut::HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut(fhicl::ParameterSet const& ps) :
    fHitsModuleLabel(ps.get<std::string>("ModuleLabels.HitModuleLabel")),
    fShowerModuleLabel(ps.get<std::string>("ModuleLabels.ShowerModuleLabel")),
    fPFParticleModuleLabel(ps.get<std::string>("ModuleLabels.PFParticleModuleLabel")),
    fShowerCompletenessCut(ps.get<double>("ShowerCompletenessCut")),
    fShowerPurityCut(ps.get<double>("ShowerPurityCut")),
    fPandrizzleAlg(ps),
    fRecoNuVtxX(-9999),
    fRecoNuVtxY(-9999),
    fRecoNuVtxZ(-9999)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Shower> FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut::SelectShower(art::Event const & evt){

  art::Ptr<recob::Shower> selShower;

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerList;
  if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle)))
  {
    std::cout<<"Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
    return selShower;
  }

  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitList, hitListHandle);

  art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);

  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

  // Check we have one neutrino...
  if (nu_pfps.size() != 1)
  {
    std::cout << "HighestEnergyRecoVertexShowerSelectorCheatCut: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
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
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  //Loop over each PFP, find the associated showers and then find which one has the highest pandrizzle score

  double highestPandrizzleScore = std::numeric_limits<double>::lowest();

  for (int i_child = 0; i_child < nu_pfp->NumDaughters(); i_child++)
  {
    //Use the PFParticle map to get the child PFPs
    int id_child = nu_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[id_child];

    //Now get the associated shower 
    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, fShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(child_pfp.key());

    if (pfp_shower_vector.size() > 1)
    { 
      //Found a PFP with more than one shower matched.  Complain and exit
        std::cout<< "Number of associated showers to a PFP is greater than 1: " << pfp_shower_vector.size() << std::endl;
        continue;
    } 
    else if (pfp_shower_vector.size() == 0)
    { 
      //Don't bother complaining of no shower was found.  It's either missing or its a track.  No biggie here
      continue;
    }

    const art::Ptr<recob::Shower> shower = pfp_shower_vector[0];
    const std::vector<art::Ptr<recob::Hit> > showerHits = fmhs.at(shower.key());

    /////////////////////
    // Calculate completeness and purity
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, 1);

    double showerCompleteness = -9999.0, showerPurity = -9999.0;

    if (g4id > 0)
    {
        showerCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, showerHits, hitList, g4id);
        showerPurity = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, showerHits, g4id);
    }

    if ((showerCompleteness < fShowerCompletenessCut) || (showerPurity < fShowerPurityCut))
        continue;
    ////////////////////////////////

    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(shower, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), evt));
    double pandrizzleScore = pandrizzleRecord.GetMVAScore();

    if (pandrizzleScore > highestPandrizzleScore)
    {
        highestPandrizzleScore = pandrizzleScore;
        selShower = shower;
    }
  }

  return selShower;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



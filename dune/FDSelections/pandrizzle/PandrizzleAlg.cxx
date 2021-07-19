///////////////////////////////////////////////
// PandrizzleAlg.cxx
//
// D Brailsford
///////////////////////////////////////////////

//STL
#include <limits>

//ART
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "lardataobj/AnalysisBase/MVAPIDResult.h"

//Custom
#include "PandrizzleAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace
{
    constexpr Float_t kDefValue(std::numeric_limits<Float_t>::lowest());

    using namespace FDSelection;
    void Reset(PandrizzleAlg::InputVarsToReader &inputVarsToReader)
    {
        for (PandrizzleAlg::Vars var = PandrizzleAlg::kEvalRatio; var < PandrizzleAlg::kTerminatingValue; var=static_cast<PandrizzleAlg::Vars>(static_cast<int>(var)+1)) 
        {
            auto [itr, inserted] = inputVarsToReader.try_emplace(var, std::make_unique<Float_t>(kDefValue));
            if (!inserted)
                *(itr->second.get()) = kDefValue;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandrizzleAlg::Record::Record(const InputVarsToReader &inputVarsToReader, const Float_t mvaScore, const bool isFilled) :
    fMVAScore(mvaScore),
    fIsFilled(isFilled)
{
    for (PandrizzleAlg::Vars var = PandrizzleAlg::kEvalRatio; var < PandrizzleAlg::kTerminatingValue; var=static_cast<PandrizzleAlg::Vars>(static_cast<int>(var)+1)) 
        fInputs.try_emplace(var, *(inputVarsToReader.at(var)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandrizzleAlg::PandrizzleAlg(const fhicl::ParameterSet& pset) :
    fPFParticleModuleLabel(pset.get<std::string>("ModuleLabels.PFParticleModuleLabel")),
    fShowerModuleLabel(pset.get<std::string>("ModuleLabels.ShowerModuleLabel")),
    fClusterModuleLabel(pset.get<std::string>("ModuleLabels.ClusterModuleLabel")),
    fPIDModuleLabel(pset.get<std::string>("ModuleLabels.PIDModuleLabel")),
    fPandrizzleWeightFileName(pset.get< std::string > ("PandrizzleWeightFileName")),
    fCheatCharacterisation(pset.get<bool>("CheatCharacterisation", false)),
    fReader("", 0),
    fMakeTree(pset.get<bool>("MakeTree", false))
{
    if (fCheatCharacterisation)
    {
        fCheatShowerModuleLabel = pset.get< std::string >("ModuleLabels.CheatShowerModuleLabel");
        fCheatPIDModuleLabel = pset.get<std::string>("ModuleLabels.CheatPIDModuleLabel");
        fNuGenModuleLabel = pset.get<std::string>("ModuleLabels.NuGenModuleLabel");
        fShowerPDGToCheat = pset.get<std::vector<int>> ("ShowerPDGToCheat");
    }

    Reset(fInputsToReader);
    fReader.AddVariable("EvalRatio",GetVarPtr(kEvalRatio));
    fReader.AddVariable("Concentration",GetVarPtr(kConcentration));
    fReader.AddVariable("CoreHaloRatio",GetVarPtr(kCoreHaloRatio));
    fReader.AddVariable("Conicalness",GetVarPtr(kConicalness));
    fReader.AddVariable("dEdxBestPlane",GetVarPtr(kdEdxBestPlane));
    fReader.AddVariable("Displacement",GetVarPtr(kDisplacement));
    fReader.AddVariable("DCA",GetVarPtr(kDCA));
    fReader.AddVariable("Wideness",GetVarPtr(kWideness));
    fReader.AddVariable("EnergyDensity",GetVarPtr(kEnergyDensity));
    const std::string weightFileName(fPandrizzleWeightFileName);
    std::string weightFilePath;
    cet::search_path sP("FW_SEARCH_PATH");
    sP.find_file(weightFileName, weightFilePath);

    fReader.BookMVA("BDTG",weightFilePath);

    if (fMakeTree)
    {
        InitialiseTrees();
        ResetTreeVariables();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::InitialiseTrees() {
  fSignalShowerTree = tfs->make<TTree>("DrizzleSigShowerTree","Pandrizzle Signal Shower Tree");
  fBackgroundShowerTree = tfs->make<TTree>("DrizzleBgShowerTree","Pandrizzle Background Shower Tree");

  std::map<std::string, TTree*> treeMap;
  treeMap["signalShower"] = fSignalShowerTree;
  treeMap["backgroundShower"] = fBackgroundShowerTree;
  for (std::map<std::string, TTree*>::iterator mapIt = treeMap.begin(); mapIt != treeMap.end(); mapIt++){
    TTree *tree = mapIt->second;
    BookTreeInt(tree, "Event");
    BookTreeInt(tree, "Run");
    BookTreeInt(tree, "SubRun");
    BookTreeInt(tree, "TruePDG");
    BookTreeFloat(tree, "TrueEnergy");
    BookTreeInt(tree, "PFPPDG");
    BookTreeInt(tree, "PFPNHits");
    BookTreeFloat(tree, "EvalRatio");
    BookTreeFloat(tree, "Concentration");
    BookTreeFloat(tree, "CoreHaloRatio");
    BookTreeFloat(tree, "Conicalness");
    BookTreeFloat(tree, "dEdxBestPlane");
    BookTreeFloat(tree, "Displacement");
    BookTreeFloat(tree, "DCA");
    BookTreeFloat(tree, "Wideness");
    BookTreeFloat(tree, "EnergyDensity");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::BookTreeInt(TTree *tree, std::string branch_name)
{
    fVarHolder.IntVars[branch_name] = -9999;
    tree->Branch(branch_name.c_str(), &(fVarHolder.IntVars[branch_name]));
    return;
}

void FDSelection::PandrizzleAlg::BookTreeFloat(TTree *tree, std::string branch_name)
{
    fVarHolder.FloatVars[branch_name] = -9999.f;
    tree->Branch(branch_name.c_str(), &(fVarHolder.FloatVars[branch_name]));
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::Run(const art::Event& evt) {

  //Get the PFPs out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)))
  {
      mf::LogWarning("PandrizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
      return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  //Ceate the full PFP map
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  //Grab the primary PFPs (the neutrinos) from the 
  std::vector<art::Ptr<recob::PFParticle> > nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

  //Now grab the primary children of these PFP
  for (unsigned int i_nupfp = 0; i_nupfp < nu_pfps.size(); i_nupfp++)
  {
      art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[i_nupfp];
      std::vector<art::Ptr<recob::PFParticle> > child_pfps = SelectChildPFParticles(nu_pfp, pfparticleMap);

      //Assess each child pfp
      for (unsigned int i_childpfp = 0; i_childpfp < child_pfps.size(); i_childpfp++)
      {
          art::Ptr<recob::PFParticle> child_pfp = child_pfps[i_childpfp];

          //Process the child PFP
          ProcessPFParticle(child_pfp, evt);
          FillTree();
          ResetTreeVariables();
      }
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<art::Ptr<recob::PFParticle> > FDSelection::PandrizzleAlg::SelectChildPFParticles(const art::Ptr<recob::PFParticle> parent_pfp, const lar_pandora::PFParticleMap & pfp_map)
{
  std::vector<art::Ptr<recob::PFParticle> > child_pfps;

  for (int i_child = 0; i_child < parent_pfp->NumDaughters(); i_child++)
  {
      int child_id = parent_pfp->Daughter(i_child);
      art::Ptr<recob::PFParticle> child_pfp = pfp_map.at(child_id);
      child_pfps.push_back(child_pfp);
  }

  return child_pfps;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
    art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
    if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
        art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

    //Get the matched MCParticle
    simb::MCParticle* matched_mcparticle(nullptr);
    std::vector<art::Ptr<recob::Hit> > pfp_hits = GetPFPHits(pfp, evt);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfp_hits, 1);
    if (g4id > 0)
    {
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        matched_mcparticle = pi_serv->ParticleList().at(g4id);
    }

    std::string showerModuleLabel(fShowerModuleLabel);
    std::string pidModuleLabel(fPIDModuleLabel);

    if (matched_mcparticle && !fShowerPDGToCheat.empty())
    {
        const int absPdg(std::abs(matched_mcparticle->PdgCode()));

        std::cout << "absPdg: " << absPdg << std::endl;

        for (int cheatPDG : fShowerPDGToCheat)
        {
            const int absCheatPdg(std::abs(cheatPDG));

            if (absPdg != absCheatPdg)
                continue;

            if ((absCheatPdg == 11 || absCheatPdg == 13))
                continue;

            showerModuleLabel = fCheatShowerModuleLabel;
            pidModuleLabel = fCheatPIDModuleLabel;
            break;
        }

       // Get the generator record
        art::Handle<std::vector<simb::MCTruth> > mcTruthListHandle;
        std::vector<art::Ptr<simb::MCTruth> > mcList;

        if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle))
            art::fill_ptr_vector(mcList, mcTruthListHandle);

        if (mcList.size() == 1)
        {
            const bool isNC = mcList[0]->GetNeutrino().CCNC();
            const bool isNue = (std::abs(mcList[0]->GetNeutrino().Nu().PdgCode()) == 12);
            const bool isNumu = (std::abs(mcList[0]->GetNeutrino().Nu().PdgCode()) == 14);
            const bool isCCNue = (isNue && !isNC);
            const bool isCCNumu = (isNumu && !isNC);

            for (int cheatPDG : fShowerPDGToCheat)
            {
                const int absCheatPdg(std::abs(cheatPDG));

                if (absPdg != absCheatPdg)
                    continue;

                if (absPdg == 11 && isCCNue)
                {
                    showerModuleLabel = fCheatShowerModuleLabel;
                    pidModuleLabel = fCheatPIDModuleLabel;
                    break;
                }

                if (absPdg == 13 && isCCNumu)
                {
                    showerModuleLabel = fCheatShowerModuleLabel;
                    pidModuleLabel = fCheatPIDModuleLabel;
                    break;
                }
            }
        }
    }

    std::cout << "showerModuleLabel: " << showerModuleLabel << std::endl;
    std::cout << "pidModuleLabel: " << pidModuleLabel << std::endl;

    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, showerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(pfp.key());

    if (pfp_shower_vector.size() > 1)
    { 
        //Found a PFP with more than one shower matched.  Complain and exit
        std::cout<< "Number of associated showers to a PFP is greater than 1: " << pfp_shower_vector.size() << std::endl;
        return;
    } 
    else if (pfp_shower_vector.size() == 0)
    { 
        //Don't bother complaining of no shower was found.  It's either missing or its a track.  No biggie here
        return;
    }

    const art::Ptr<recob::Shower> pShower = pfp_shower_vector[0];

    lar_pandora::PFParticleVector nu_pfps;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

    // Check we have one neutrino...
    if (nu_pfps.size() != 1)
    {
        std::cout << "PandrizzleAlg: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
        return;
    }

    art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

    art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
    const std::vector<art::Ptr<recob::Vertex> > nu_vertices = fmvpfp.at(nu_pfp.key());

    if (nu_vertices.size() == 0)
    {
        return;
    }
    else if (nu_vertices.size() > 1)
    {
        std::cout << "PandrizzleAlg: Number of matched neutrino vertices bigger than 1: " << nu_vertices.size() << std::endl;
    }

    //always take the first vertex, even if there's more than one
    TVector3 nuVertex = TVector3(nu_vertices[0]->position().X(), nu_vertices[0]->position().Y(), nu_vertices[0]->position().Z());

    art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower> >{pShower}, evt, pidModuleLabel);
    art::Ptr<anab::MVAPIDResult> mvaPIDResult(findPIDResult.at(0));

    //Get event,subrun,run
    fVarHolder.IntVars["Run"] = evt.run();
    fVarHolder.IntVars["SubRun"] = evt.subRun();
    fVarHolder.IntVars["Event"] = evt.event();

    //Fill the PFP hits
    fVarHolder.IntVars["PFPPDG"] = pfp->PdgCode();
    fVarHolder.IntVars["PFPNHits"] = pfp_hits.size();

    if (matched_mcparticle)
    {
        fVarHolder.IntVars["TruePDG"] = matched_mcparticle->PdgCode();
        fVarHolder.FloatVars["TrueEnergy"] = matched_mcparticle->Momentum().T();
    }

    //MVAPID vars
    if (mvaPIDResult.isAvailable())
    {
        fVarHolder.FloatVars["EvalRatio"] = (isnan(mvaPIDResult->evalRatio) ? -0.5f : static_cast<Float_t>(mvaPIDResult->evalRatio));
        fVarHolder.FloatVars["Concentration"] = (isnan(mvaPIDResult->concentration) ? -2.f : std::min(static_cast<Float_t>(mvaPIDResult->concentration), 50.f));
        fVarHolder.FloatVars["CoreHaloRatio"] = static_cast<Float_t>(mvaPIDResult->coreHaloRatio);
        fVarHolder.FloatVars["Conicalness"] = std::min(static_cast<Float_t>(mvaPIDResult->conicalness), 100.f);
    }

    //dEdx
    if (pShower->dEdx().size() > 0)
        fVarHolder.FloatVars["dEdxBestPlane"] = std::max(std::min(static_cast<Float_t>(pShower->dEdx().at(pShower->best_plane())), 20.f), -2.f);

    //Displacement
    fVarHolder.FloatVars["Displacement"] = std::min(static_cast<Float_t>((pShower->ShowerStart() - nuVertex).Mag()), 100.f);

    //Distance of closest approach
    double alpha((pShower->ShowerStart() - nuVertex).Dot(pShower->Direction()));
    TVector3 r(pShower->ShowerStart() + alpha*pShower->Direction());

    fVarHolder.FloatVars["DCA"] = std::min(static_cast<Float_t>((r-nuVertex).Mag()), 50.f);

    //Wideness
    Float_t wideness(static_cast<Float_t>(pShower->OpenAngle()/pShower->Length()));

    fVarHolder.FloatVars["Wideness"] = isnan(wideness) ? -0.01f : std::min(wideness, 0.1f);

    //Energy density
    if (pShower->Energy().size() > 0)
    {
        Float_t volume(static_cast<Float_t>((M_PI * pShower->Length() * pShower->Length() * pShower->Length() * std::tan(pShower->OpenAngle()))/3.f));
        Float_t energyDensity(std::min(std::max(static_cast<Float_t>(pShower->Energy().at(2))/volume, -0.1f), 5.f));

        fVarHolder.FloatVars["EnergyDensity"] = isnan(energyDensity) ? -0.1f : energyDensity;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::FillTree(){
  if (std::abs(fVarHolder.IntVars["TruePDG"]) == 11)
  { 
      if (fVarHolder.IntVars["PFPPDG"] == 11)
          fSignalShowerTree->Fill();
  }
  else 
  { 
      if (fVarHolder.IntVars["PFPPDG"] == 11)
          fBackgroundShowerTree->Fill();
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::ResetTreeVariables()
{
    for (std::map<std::string, int>::iterator mapIt = fVarHolder.IntVars.begin(); mapIt != fVarHolder.IntVars.end(); mapIt++)
        mapIt->second = -9999;

    for (std::map<std::string, float>::iterator mapIt = fVarHolder.FloatVars.begin(); mapIt != fVarHolder.FloatVars.end(); mapIt++)
        mapIt->second = -9999.f;

    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<art::Ptr<recob::Hit> > FDSelection::PandrizzleAlg::GetPFPHits(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
  std::vector<art::Ptr<recob::Hit> > pfp_hits;

  //Get the PFP handle out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)))
  {
      mf::LogWarning("PandrizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
      return pfp_hits;
  }

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  if (!(evt.getByLabel(fClusterModuleLabel, clusterListHandle)))
  {
      mf::LogWarning("PandrizzleAlg") << "Unable to find std::vector<recob::Cluster> with module label: " << fClusterModuleLabel;
      return pfp_hits;
  }

  art::FindManyP<recob::Cluster> fmcpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Cluster> > sel_pfp_clusters = fmcpfp.at(pfp.key());
  art::FindManyP<recob::Hit> fmhc(clusterListHandle, evt, fClusterModuleLabel);

  //Loop over the clusters and retrieve the hits
  for (unsigned i_cluster = 0; i_cluster < sel_pfp_clusters.size(); i_cluster++)
  {
      art::Ptr<recob::Cluster> cluster = sel_pfp_clusters[i_cluster];
      const std::vector<art::Ptr<recob::Hit> > sel_pfp_hits = fmhc.at(cluster.key());

    for (unsigned i_hit = 0; i_hit < sel_pfp_hits.size(); i_hit++)
        pfp_hits.push_back(sel_pfp_hits[i_hit]);
  }

  return pfp_hits;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////


FDSelection::PandrizzleAlg::Record FDSelection::PandrizzleAlg::RunPID(const art::Ptr<recob::Shower> pShower, const TVector3& nuVertex, const art::Event& evt) 
{
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector <art::Ptr<recob::Shower>> showerList;
  if (evt.getByLabel(fShowerModuleLabel, showerListHandle))
      art::fill_ptr_vector(showerList, showerListHandle);

  bool cheat(false);
  art::Handle< std::vector<recob::Shower> > cheatShowerListHandle;
  if (fCheatCharacterisation && (std::find(showerList.begin(), showerList.end(), pShower) == showerList.end()))
    cheat = true;

  art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower> >{pShower}, evt, cheat ? fCheatPIDModuleLabel : fPIDModuleLabel);
    art::Ptr<anab::MVAPIDResult> mvaPIDResult(findPIDResult.at(0));
    //MVAPID vars
    if (mvaPIDResult.isAvailable())
    {
        if (isnan(mvaPIDResult->evalRatio))
            SetVar(kEvalRatio, -0.5f);
        else
            SetVar(kEvalRatio, static_cast<Float_t>(mvaPIDResult->evalRatio));

        if (isnan(mvaPIDResult->concentration))
            SetVar(kConcentration, -2.f);
        else
            SetVar(kConcentration, std::min(static_cast<Float_t>(mvaPIDResult->concentration), 50.f));
        SetVar(kCoreHaloRatio, static_cast<Float_t>(mvaPIDResult->coreHaloRatio));
        SetVar(kConicalness, std::min(static_cast<Float_t>(mvaPIDResult->conicalness), 100.f));
    }
    else 
        ReturnEmptyRecord();

    //dEdx
    if (pShower->dEdx().size() > 0)
    {
        SetVar(kdEdxBestPlane, std::max(std::min(static_cast<Float_t>(pShower->dEdx().at(pShower->best_plane())), 20.f), -2.f));
    }
    else
        ReturnEmptyRecord();

    //Displacement
    SetVar(kDisplacement, std::min(static_cast<Float_t>((pShower->ShowerStart() - nuVertex).Mag()), 100.f));

    //Distance of closest approach
    double alpha((pShower->ShowerStart() - nuVertex).Dot(pShower->Direction()));
    TVector3 r(pShower->ShowerStart() + alpha*pShower->Direction());
    SetVar(kDCA, std::min(static_cast<Float_t>((r-nuVertex).Mag()), 50.f));

    //Wideness
    Float_t wideness(static_cast<Float_t>(pShower->OpenAngle()/pShower->Length()));
    if (isnan(wideness))
        SetVar(kWideness, -0.01f);
    else
        SetVar(kWideness, std::min( wideness, 0.1f));

    //Energy density
    if (pShower->Energy().size() > 0)
    {
        Float_t volume(static_cast<Float_t>((M_PI * pShower->Length() * pShower->Length() * pShower->Length() * std::tan(pShower->OpenAngle()))/3.));
        Float_t energyDensity(std::min(std::max(static_cast<Float_t>(pShower->Energy().at(2))/volume, -0.1f), 5.f));
        if (isnan(energyDensity))
            energyDensity = -0.1f;
        SetVar(kEnergyDensity, energyDensity);
    }
    else
        ReturnEmptyRecord();

    return Record(fInputsToReader, fReader.EvaluateMVA("BDTG"), true);
}

FDSelection::PandrizzleAlg::Record FDSelection::PandrizzleAlg::ReturnEmptyRecord()
{
    Reset(fInputsToReader);
    return Record(fInputsToReader, kDefValue, false);
}

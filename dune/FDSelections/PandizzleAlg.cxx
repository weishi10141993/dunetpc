///////////////////////////////////////////////
// PandizzleAlg.cxx
//
// Reco and true PID stuff up
// D Brailsford & M Wallbank, June 2017
///////////////////////////////////////////////

#include "PandizzleAlg.h"

FDSelection::PandizzleAlg::PandizzleAlg(const fhicl::ParameterSet& pset) {
  fTrackModuleLabel  = pset.get<std::string>("TrackModuleLabel");
  fShowerModuleLabel = pset.get<std::string>("ShowerModuleLabel");
  fPIDModuleLabel    = pset.get<std::string>("PIDModuleLabel");
  fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
  fSpacePointModuleLabel = pset.get<std::string>("SpacePointModuleLabel");
  InitialiseTrees();
  ResetTreeVariables();
}

void FDSelection::PandizzleAlg::InitialiseTrees() {
  fSignalTrackTree = tfs->make<TTree>("DizzleSigTrackTree","Pandizzle Signal Track Tree");
  fBackgroundTrackTree = tfs->make<TTree>("DizzleBgTrackTree","Pandizzle Background Track Tree");
  fSignalShowerTree = tfs->make<TTree>("DizzleSigShowerTree","Pandizzle Signal Shower Tree");
  fBackgroundShowerTree = tfs->make<TTree>("DizzleBgShowerTree","Pandizzle Background Shower Tree");

  //I am lazy.  Sue me.
  std::map<std::string, TTree*> treeMap;
  treeMap["signalTrack"] = fSignalTrackTree;
  treeMap["backgroundTrack"] = fBackgroundTrackTree;
  treeMap["signalShower"] = fSignalShowerTree;
  treeMap["backgroundShower"] = fBackgroundShowerTree;
  for (std::map<std::string, TTree*>::iterator mapIt = treeMap.begin(); mapIt != treeMap.end(); mapIt++){
    TTree *tree = mapIt->second;
    BookTreeInt(tree, "Event");
    BookTreeInt(tree, "Run");
    BookTreeInt(tree, "SubRun");
    BookTreeInt(tree, "PFPPDG");
    BookTreeInt(tree, "PFPNHits");
    BookTreeInt(tree, "PFPTrueID");
    BookTreeInt(tree, "PFPTruePDG");
    BookTreeInt(tree, "PFPTrueMotherID");
    BookTreeFloat(tree, "PFPTrueMomT");
    BookTreeFloat(tree, "PFPTrueMomX");
    BookTreeFloat(tree, "PFPTrueMomY");
    BookTreeFloat(tree, "PFPTrueMomZ");
    BookTreeInt(tree, "PFPNChildren");
    BookTreeInt(tree, "PFPNShowerChildren");
    BookTreeInt(tree, "PFPNTrackChildren");
    BookTreeFloat(tree, "PFPMichelDist");
    BookTreeInt(tree, "PFPMichelNHits");
    BookTreeInt(tree, "PFPMichelTrueID");
    BookTreeInt(tree, "PFPMichelTrueMotherID");
    BookTreeInt(tree, "PFPMichelTruePDG");






  }
}

void FDSelection::PandizzleAlg::Run(const art::Event& evt) {

  //Get the PFPs out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
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
  for (unsigned int i_nupfp = 0; i_nupfp < nu_pfps.size(); i_nupfp++){
    art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[i_nupfp];
    std::vector<art::Ptr<recob::PFParticle> > child_pfps = SelectChildPFParticles(nu_pfp, pfparticleMap);
    //Assess each child pfp
    for (unsigned int i_childpfp = 0; i_childpfp < child_pfps.size(); i_childpfp++){
      art::Ptr<recob::PFParticle> child_pfp = child_pfps[i_childpfp];
      //Process the child PFP
      ProcessPFParticle(child_pfp, evt);
    }
  }
  //mf::LogWarning("PandizzleAlg") << "No tracks";

  return;

}

std::vector<art::Ptr<recob::PFParticle> > FDSelection::PandizzleAlg::SelectChildPFParticles(const art::Ptr<recob::PFParticle> parent_pfp, const lar_pandora::PFParticleMap & pfp_map){
  std::vector<art::Ptr<recob::PFParticle> > child_pfps;
  for (int i_child = 0; i_child < parent_pfp->NumDaughters(); i_child++){
    int child_id = parent_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfp_map.at(child_id);
    child_pfps.push_back(child_pfp);
  }
  return child_pfps;
}

void FDSelection::PandizzleAlg::ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt){
  //Get event,subrun,run
  fVarHolder.IntVars["Run"] = evt.run();
  fVarHolder.IntVars["SubRun"] = evt.subRun();
  fVarHolder.IntVars["Event"] = evt.event();

  //Get the PFPs out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
    return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  //Ceate the full PFP map
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  //PDG
  fVarHolder.IntVars["PFPPDG"] = pfp->PdgCode();


  //Get the PFP hits
  std::vector<art::Ptr<recob::Hit> > pfp_hits = GetPFPHits(pfp, evt);
  fVarHolder.IntVars["PFPNHits"] = pfp_hits.size();

  //Get the true PDG
  int g4id = FDSelectionUtils::TrueParticleIDFromTotalRecoHits(pfp_hits); 
  if (g4id > 0){
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);

    if (matched_mcparticle){
      fVarHolder.IntVars["PFPTrueID"] = g4id;
      fVarHolder.IntVars["PFPTruePDG"] = matched_mcparticle->PdgCode();
      fVarHolder.IntVars["PFPTrueMotherID"] = matched_mcparticle->Mother();
      fVarHolder.FloatVars["PFPTrueMomX"] = matched_mcparticle->Momentum(0).X();
      fVarHolder.FloatVars["PFPTrueMomY"] = matched_mcparticle->Momentum(0).Y();
      fVarHolder.FloatVars["PFPTrueMomZ"] = matched_mcparticle->Momentum(0).Z();
      fVarHolder.FloatVars["PFPTrueMomT"] = matched_mcparticle->Momentum(0).T();
    }
  }

  //Count all of the children
  std::vector<art::Ptr<recob::PFParticle> > child_pfps = SelectChildPFParticles(pfp, pfparticleMap);
  fVarHolder.IntVars["PFPNChildren"] = pfp->NumDaughters();
  fVarHolder.IntVars["PFPNShowerChildren"] = CountPFPWithPDG(child_pfps, 11);
  fVarHolder.IntVars["PFPNTrackChildren"] = CountPFPWithPDG(child_pfps, 13);
  
  FillMichelElectronVariables(pfp, child_pfps, evt);
  /*
  if (fVarHolder.IntVars["PFPTruePDG"] == 13 && fVarHolder.IntVars["PFPPDG"] ==13){
    for (unsigned int i_child = 0; i_child < child_pfps.size(); i_child++){
      art::Ptr<recob::PFParticle> child_pfp = child_pfps[i_child];
      std::vector<art::Ptr<recob::Hit> > child_pfp_hits = GetPFPHits(child_pfp, evt); 
      int child_g4id = FDSelectionUtils::TrueParticleIDFromTotalTrueEnergy(child_pfp_hits);
      if (child_g4id <= 0) {
        std::cout<<"child PFP has no true ID"<<std::endl;
        continue;
      }
      const simb::MCParticle* matched_childmcparticle = pi_serv->ParticleList().at(child_g4id);

      if (matched_childmcparticle){
        std::cout<<"child " << i_child << "  True PDG: " << matched_childmcparticle->PdgCode() << "   PFP PDG: " << child_pfp->PdgCode() << std::endl;
      }
    }
  } 
  */
  FillTree();
  ResetTreeVariables();
  return;
}

void FDSelection::PandizzleAlg::ResetTreeVariables(){
  for (std::map<std::string, int>::iterator mapIt = fVarHolder.IntVars.begin(); mapIt != fVarHolder.IntVars.end(); mapIt++){
    mapIt->second = -9999;
  }
  for (std::map<std::string, float>::iterator mapIt = fVarHolder.FloatVars.begin(); mapIt != fVarHolder.FloatVars.end(); mapIt++){
    mapIt->second = -9999.;
  }
  return;
}

int FDSelection::PandizzleAlg::CountPFPWithPDG(const std::vector<art::Ptr<recob::PFParticle> > & pfps, int pdg){
  int NPFP = 0;
  for (unsigned int i_pfp = 0; i_pfp < pfps.size(); i_pfp++){
    int pfp_pdg = pfps[i_pfp]->PdgCode();
    if (pfp_pdg == pdg) NPFP++;
  }
  return NPFP;
}

std::vector<art::Ptr<recob::Hit> > FDSelection::PandizzleAlg::GetPFPHits(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt){
  std::vector<art::Ptr<recob::Hit> > pfp_hits;
  //Get the PFP handle out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
    return pfp_hits; //empty
  }
  //Get the spacepoint handle out of the event
  art::Handle< std::vector<recob::SpacePoint> > spacePointListHandle;
  if (!(evt.getByLabel(fSpacePointModuleLabel, spacePointListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::SpacePoint> with module label: " << fSpacePointModuleLabel;
    return pfp_hits; //empty
  }

  art::FindManyP<recob::SpacePoint> fmsppfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::SpacePoint> > sel_pfp_spacepoints = fmsppfp.at(pfp.key());
  art::FindManyP<recob::Hit> fmhsp(spacePointListHandle, evt, fSpacePointModuleLabel);
  //Loop over the space points and retrieve the hits
  for (unsigned i_sp = 0; i_sp < sel_pfp_spacepoints.size(); i_sp++){
    art::Ptr<recob::SpacePoint> sp = sel_pfp_spacepoints[i_sp];
    const std::vector<art::Ptr<recob::Hit> > sel_pfp_hits = fmhsp.at(sp.key());
    for (unsigned i_hit = 0; i_hit < sel_pfp_hits.size(); i_hit++){
      pfp_hits.push_back(sel_pfp_hits[i_hit]);
    }
  }
  return pfp_hits;
}

void FDSelection::PandizzleAlg::BookTreeInt(TTree *tree, std::string branch_name){
  fVarHolder.IntVars[branch_name] = -9999;
  tree->Branch(branch_name.c_str(), &(fVarHolder.IntVars[branch_name]));
  return;
}

void FDSelection::PandizzleAlg::BookTreeFloat(TTree *tree, std::string branch_name){
  fVarHolder.FloatVars[branch_name] = -9999.;
  tree->Branch(branch_name.c_str(), &(fVarHolder.FloatVars[branch_name]));
  return;
}
void FDSelection::PandizzleAlg::FillTree(){
  if (fVarHolder.IntVars["PFPTruePDG"] == 13){ //signal
    if (fVarHolder.IntVars["PFPPDG"] == 13){ //track
      fSignalTrackTree->Fill();
    }
    else if (fVarHolder.IntVars["PFPPDG"] == 11){ //shower
      fSignalShowerTree->Fill();
    }
    else{ //Don't know
      mf::LogWarning("PandizzleAlg") << "Unknown PFP PDG when filling tree"<< fVarHolder.IntVars["PFPPDG"];
    }
  }
  else { //background
    if (fVarHolder.IntVars["PFPPDG"] == 13){ //track
      fBackgroundTrackTree->Fill();
    }
    else if (fVarHolder.IntVars["PFPPDG"] == 11){ //shower
      fBackgroundShowerTree->Fill();
    }
    else{ //Don't know
      mf::LogWarning("PandizzleAlg") << "Unknown PFP PDG when filling tree"<< fVarHolder.IntVars["PFPPDG"];

    }
  }
  return;
}


void FDSelection::PandizzleAlg::FillMichelElectronVariables(const art::Ptr<recob::PFParticle> mu_pfp, const std::vector<art::Ptr<recob::PFParticle> > & child_pfps, const art::Event& evt){
  //Find closest PFP to end of mu
  //If we don't have a track, sack everything off
  if (mu_pfp->PdgCode() != 13){
    return;
  }
  //Assign the branch value to a new default value to indicate that we have a track, but don't necessarily have a Michel candidate (rather than a global default value)
  fVarHolder.FloatVars["PFPMichelDist"] = -100.;
  fVarHolder.IntVars["PFPMichelNHits"] = -100;

  //If no child PFPs, sack everything off
  if (child_pfps.size() == 0){
    return;
  }
  //Get the PFP handle
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
    return;
  }

  art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, fShowerModuleLabel);

  const std::vector<art::Ptr<recob::Track> > sel_pfp_tracks = fmtpfp.at(mu_pfp.key());
  if (sel_pfp_tracks.size() != 1){
    return;
  }
  art::Ptr<recob::Track> mu_track = sel_pfp_tracks[0];
  TVector3 mu_end_position = mu_track->End();

  art::Ptr<recob::PFParticle> closest_child_pfp;
  double closest_distance = 99999999999;
  for (unsigned int i_child = 0; i_child < child_pfps.size(); i_child++){
    art::Ptr<recob::PFParticle> child_pfp = child_pfps[i_child];
    TVector3 child_start_pos;
    if (child_pfp->PdgCode()==13){
      const std::vector<art::Ptr<recob::Track> > child_pfp_tracks = fmtpfp.at(child_pfp.key());
      if (child_pfp_tracks.size() != 1) return;
      art::Ptr<recob::Track> child_track = child_pfp_tracks[0];
      child_start_pos.SetXYZ(child_track->Start().X(),child_track->Start().Y(),child_track->Start().Z());

    }
    else if (child_pfp->PdgCode()==11){
      const std::vector<art::Ptr<recob::Shower> > child_pfp_showers = fmspfp.at(child_pfp.key());
      if (child_pfp_showers.size() != 1) return;
      art::Ptr<recob::Shower> child_shower = child_pfp_showers[0];
      child_start_pos.SetXYZ(child_shower->ShowerStart().X(),child_shower->ShowerStart().Y(),child_shower->ShowerStart().Z());
    }
    else {
      mf::LogWarning("PandizzleAlg") << "Unknown PFP PDG when finding Michel candidate"<< child_pfp->PdgCode();

    }
    double curr_distance = (child_start_pos-mu_end_position).Mag();
    if (curr_distance < closest_distance){
      closest_distance = curr_distance;
      closest_child_pfp = child_pfp;
    }
  }

  fVarHolder.FloatVars["PFPMichelDist"] = closest_distance;

  std::vector<art::Ptr<recob::Hit> > michel_hits = GetPFPHits(closest_child_pfp, evt);
  fVarHolder.IntVars["PFPMichelNHits"] = michel_hits.size();
  int g4id = FDSelectionUtils::TrueParticleIDFromTotalRecoHits(michel_hits); 
  if (g4id > 0){
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
    if (matched_mcparticle){
      fVarHolder.IntVars["PFPMichelTrueID"] = matched_mcparticle->TrackId();
      fVarHolder.IntVars["PFPMichelTrueMotherID"] = matched_mcparticle->Mother();
      fVarHolder.IntVars["PFPMichelTruePDG"] = matched_mcparticle->PdgCode();
    }
  }

  return;
}


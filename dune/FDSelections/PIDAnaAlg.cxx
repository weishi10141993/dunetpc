///////////////////////////////////////////////
// PIDAnaAlg.cxx
//
// Reco and true PID stuff up
// D Brailsford & M Wallbank, June 2017
///////////////////////////////////////////////

#include "PIDAnaAlg.h"

FDSelection::PIDAnaAlg::PIDAnaAlg(const fhicl::ParameterSet& pset) {
  fTrackModuleLabel  = pset.get<std::string>("TrackModuleLabel");
  fShowerModuleLabel = pset.get<std::string>("ShowerModuleLabel");
  fPIDModuleLabel    = pset.get<std::string>("PIDModuleLabel");
  this->InitialiseTree();
}

void FDSelection::PIDAnaAlg::InitialiseTree() {

  fTree = tfs->make<TTree>("PIDAna","PIDAna");
  fTree->Branch("Run",         &fRun);
  fTree->Branch("SubRun",      &fSubRun);
  fTree->Branch("NObjects",    &fNObjects);
  fTree->Branch("Track",       fTrack,       "fTrack[NObjects]/O");
  fTree->Branch("Shower",      fShower,      "fShower[NObjects]/O");
  fTree->Branch("TruePDG",     fTruePDG,     "fTruePDG[NObjects]/I");
  fTree->Branch("RecoPDG",     fRecoPDG,     "fRecoPDG[NObjects]/I");
  fTree->Branch("ElectronMVA", fElectronMVA, "fElectronMVA[NObjects]/D");
  fTree->Branch("MuonMVA",     fMuonMVA,     "fMuonMVA[NObjects]/D");
  fTree->Branch("PhotonMVA",   fPhotonMVA,   "fPhotonMVA[NObjects]/D");
  fTree->Branch("ProtonMVA",   fProtonMVA,   "fProtonMVA[NObjects]/D");
  fTree->Branch("PionMVA",     fPionMVA,     "fPionMVA[NObjects]/D");
  fTree->Branch("ObjPurity",   fObjPurity,   "fObjPurity[NObjects]/D");
  fTree->Branch("Primary",     fPrimary,     "fPrimary[NObjects]/B");
  fTree->Branch("TrueEnergy",  fTrueEnergy,  "fTrueEnergy[NObjects]/D");
  fTree->Branch("RecoEnergy",  fRecoEnergy,  "fRecoEnergy[NObjects]/D");
  fTree->Branch("TrueX",       fTrueX,       "fTrueX[NObjects]/D");
  fTree->Branch("TrueY",       fTrueY,       "fTrueY[NObjects]/D");
  fTree->Branch("TrueZ",       fTrueZ,       "fTrueZ[NObjects]/D");
  fTree->Branch("TrueEndX",    fTrueEndX,    "fTrueEndX[NObjects]/D");
  fTree->Branch("TrueEndY",    fTrueEndY,    "fTrueEndY[NObjects]/D");
  fTree->Branch("TrueEndZ",    fTrueEndZ,    "fTrueEndZ[NObjects]/D");
  fTree->Branch("RecoX",       fRecoX,       "fRecoX[NObjects]/D");
  fTree->Branch("RecoY",       fRecoY,       "fRecoY[NObjects]/D");
  fTree->Branch("RecoZ",       fRecoZ,       "fRecoZ[NObjects]/D");
  fTree->Branch("RecoEndX",    fRecoEndX,    "fRecoEndX[NObjects]/D");
  fTree->Branch("RecoEndY",    fRecoEndY,    "fRecoEndY[NObjects]/D");
  fTree->Branch("RecoEndZ",    fRecoEndZ,    "fRecoEndZ[NObjects]/D");
  fTree->Branch("RecoLength",  fRecoLength,  "fRecoLength[NObjects]/D");
  fTree->Branch("RecoPoints",  fRecoPoints,  "fRecoPoints[NObjects]/I");

}

void FDSelection::PIDAnaAlg::Run(const art::Event& evt) {

  //Grab the run and subrun info before doing anything
  fRun = evt.run();
  fSubRun = evt.subRun();

  std::cout << std::endl << "Running PIDAna for event " << evt.event() << std::endl;

  // Get tracks
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel(fTrackModuleLabel, trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);
  else
    mf::LogWarning("PIDAnaAlg") << "No tracks";

  // Get showers
  art::Handle<std::vector<recob::Shower> > showerHandle;
  std::vector<art::Ptr<recob::Shower> > showers;
  if (evt.getByLabel(fShowerModuleLabel, showerHandle))
    art::fill_ptr_vector(showers, showerHandle);
  else
    mf::LogWarning("PIDAnaAlg") << "No showers";

  // Hit associations
  art::FindManyP<recob::Hit> fmht(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit> fmhs(showerHandle, evt, fShowerModuleLabel);

  // PID particle associations
  art::FindManyP<anab::MVAPIDResult> fmpidt(trackHandle, evt, fPIDModuleLabel);
  art::FindManyP<anab::MVAPIDResult> fmpids(showerHandle, evt, fPIDModuleLabel);

  // true particles
  const sim::ParticleList& trueParticleList = pi_serv->ParticleList();

  // reco
  //fNObjects = (int)tracks.size() + (int)showers.size(); // for the moment, don't assume all objects have associated pid
  int i_obj = 0;

  // tracks
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {
    const std::vector<art::Ptr<recob::Hit> > hits = fmht.at(trackIt->key());
    int trueParticle = TrueParticle(hits);
    if (!trueParticleList.HasParticle(trueParticle))
      continue;
    const std::map<int,double> trueParticles = TrueParticles(hits);
    double chargeTot = 0;
    for (std::map<int,double>::const_iterator trueParticleIt = trueParticles.begin(); trueParticleIt != trueParticles.end(); ++trueParticleIt)
      chargeTot += trueParticleIt->second;
    double purity = trueParticles.at(trueParticle)/chargeTot;
    std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(trackIt->key());
    // std::cout << std::endl << "Track " << trackIt->key() << std::endl;
    // std::cout << "  Associated with particles:" << std::endl;
    // for (std::map<int,double>::const_iterator trueParticlesIt = trueParticles.begin(); trueParticlesIt != trueParticles.end(); ++trueParticlesIt)
    //   std::cout << "    " << trueParticlesIt->first << " (PDG " << trueParticleList[trueParticlesIt->first]->PdgCode() << ") contributes charge " << trueParticlesIt->second << std::endl;
    // std::cout << "  True particle assumed to be " << trueParticle << " (PDG " << trueParticleList[trueParticle]->PdgCode() << ")" << std::endl;
    // std::cout << "  There are " << pids.size() << " mva pid result objects association with track " << trackIt->key() << std::endl;
    // for (std::vector<art::Ptr<anab::MVAPIDResult> >::const_iterator pidIt = pids.begin(); pidIt != pids.end(); ++pidIt) {
    //   const std::map<std::string,double> mvaOutput = (*pidIt)->mvaOutput;
    //   std::cout << "    PID " << pidIt->key() << " has " << mvaOutput.size() << " MVA outputs" << std::endl;
    //   for (std::map<std::string,double>::const_iterator mvaOutIt = mvaOutput.begin(); mvaOutIt != mvaOutput.end(); ++mvaOutIt)
    // 	std::cout << "      MVA out " << std::distance<std::map<std::string,double>::const_iterator>(mvaOutput.begin(),mvaOutIt) << " is " << mvaOutIt->first << " : " << mvaOutIt->second << std::endl;
    // }
    if (!pids.size())
      continue;
    const simb::MCParticle* mcparticle = trueParticleList.at(trueParticle);
    std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
    // fill tree variables
    fTrack[i_obj]       = true;
    fShower[i_obj]   	= false;
    fTruePDG[i_obj]    	= mcparticle->PdgCode();
    fRecoPDG[i_obj]    	= 13;//tmp
    fElectronMVA[i_obj] = mvaOutMap["electron"];
    fMuonMVA[i_obj]	= mvaOutMap["muon"]; 
    fPhotonMVA[i_obj]   = mvaOutMap["photon"];
    fProtonMVA[i_obj]   = mvaOutMap["proton"];
    fPionMVA[i_obj]     = mvaOutMap["pich"];
    fObjPurity[i_obj]   = purity;
    fPrimary[i_obj]     = mcparticle->Process() == "primary";
    fTrueEnergy[i_obj]  = mcparticle->E();
    fRecoEnergy[i_obj]  = (*trackIt)->StartMomentum();
    fTrueX[i_obj]       = mcparticle->Position().X();
    fTrueY[i_obj]       = mcparticle->Position().Y();
    fTrueZ[i_obj]       = mcparticle->Position().Z();
    fTrueEndX[i_obj]    = mcparticle->Position(mcparticle->NumberTrajectoryPoints()-1).X();
    fTrueEndY[i_obj]    = mcparticle->Position(mcparticle->NumberTrajectoryPoints()-1).Y();
    fTrueEndZ[i_obj]    = mcparticle->Position(mcparticle->NumberTrajectoryPoints()-1).Z();
    fRecoX[i_obj]       = (*trackIt)->Vertex().X();
    fRecoY[i_obj]       = (*trackIt)->Vertex().Y();
    fRecoZ[i_obj]       = (*trackIt)->Vertex().Z();
    fRecoEndX[i_obj]    = (*trackIt)->End().X();
    fRecoEndY[i_obj]    = (*trackIt)->End().Y();
    fRecoEndZ[i_obj]    = (*trackIt)->End().Z();
    fRecoLength[i_obj]  = (*trackIt)->Length();
    fRecoPoints[i_obj]  = (*trackIt)->NumberTrajectoryPoints();
    ++i_obj;
  }

  // showers
  for (std::vector<art::Ptr<recob::Shower> >::const_iterator showerIt = showers.begin(); showerIt != showers.end(); ++showerIt) {
    const std::vector<art::Ptr<recob::Hit> > hits = fmhs.at(showerIt->key());
    int trueParticle = TrueParticle(hits);
    if (!trueParticleList.HasParticle(trueParticle))
      continue;
    const std::map<int,double> trueParticles = TrueParticles(hits);
    double chargeTot = 0;
    for (std::map<int,double>::const_iterator trueParticleIt = trueParticles.begin(); trueParticleIt != trueParticles.end(); ++trueParticleIt)
      chargeTot += trueParticleIt->second;
    double purity = trueParticles.at(trueParticle)/chargeTot;
    std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpids.at(showerIt->key());
    if (!pids.size())
      continue;
    const simb::MCParticle* mcparticle = trueParticleList.at(trueParticle);
    std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
    // fill tree variables
    fTrack[i_obj]       = false;
    fShower[i_obj]   	= true;
    fTruePDG[i_obj]    	= mcparticle->PdgCode();
    fRecoPDG[i_obj]    	= 11;//tmp
    fElectronMVA[i_obj] = mvaOutMap["electron"];
    fMuonMVA[i_obj]	= mvaOutMap["muon"]; 
    fPhotonMVA[i_obj]   = mvaOutMap["photon"];
    fProtonMVA[i_obj]   = mvaOutMap["proton"];
    fPionMVA[i_obj]     = mvaOutMap["pich"];
    fObjPurity[i_obj]   = purity;
    fPrimary[i_obj]     = mcparticle->Process() == "primary";
    fTrueEnergy[i_obj]  = mcparticle->E();
    if ((*showerIt)->Energy().size() == 0)
      fRecoEnergy[i_obj] = 0;
    else
      fRecoEnergy[i_obj] = (*showerIt)->Energy().at((*showerIt)->best_plane());
    fTrueX[i_obj]       = mcparticle->Position().X();
    fTrueY[i_obj]       = mcparticle->Position().Y();
    fTrueZ[i_obj]       = mcparticle->Position().Z();
    fRecoX[i_obj]       = (*showerIt)->ShowerStart().X();
    fRecoY[i_obj]       = (*showerIt)->ShowerStart().Y();
    fRecoZ[i_obj]       = (*showerIt)->ShowerStart().Z();
    fRecoEndX[i_obj]    = -999.;  //Does not seem to exist for showers right now  
    fRecoEndY[i_obj]    = -999.;  //Does not seem to exist for showers right now
    fRecoEndZ[i_obj]    = -999.;  //Does not seem to exist for showers right now
    fRecoLength[i_obj]  = -999.;//tmp, this has only just been added to recob::Shower so older code won't have it
    fRecoPoints[i_obj]  = hits.size();
    ++i_obj;
  }

  fNObjects = i_obj;

  fTree->Fill();

  return;

}

std::map<int,double> FDSelection::PIDAnaAlg::TrueParticles(const std::vector<art::Ptr<recob::Hit> >& hits) {

  // Make a map of the tracks which are associated with this object and the charge each contributes
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = ParticleID(hit);
    trackMap[trackID] += hit->Integral();
  }

  return trackMap;

}

int FDSelection::PIDAnaAlg::TrueParticle(const std::vector<art::Ptr<recob::Hit> >& hits) {

  std::map<int,double> trackMap = TrueParticles(hits);

  // Pick the track with the highest charge as the 'true track'
  double highestCharge = 0;
  int objectTrack = 0;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      objectTrack  = trackIt->first;
    }
  }

  return objectTrack;

}

int FDSelection::PIDAnaAlg::ParticleID(const art::Ptr<recob::Hit>& hit) {

  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }

  return likelyTrackID;

}

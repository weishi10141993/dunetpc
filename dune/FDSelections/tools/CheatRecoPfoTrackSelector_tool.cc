#include "CheatRecoPfoTrackSelector.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Utils/TruthMatchUtils.h"

FDSelectionTools::CheatRecoPfoTrackSelector::CheatRecoPfoTrackSelector(fhicl::ParameterSet const& ps) :
    fNuGenModuleLabel(ps.get< std::string> ("ModuleLabels.NuGenModuleLabel")),
    fAllTrackModuleLabel(ps.get< std::string> ("AllTrackModuleLabel")),
    fFittedTrackModuleLabel(ps.get< std::string> ("FittedTrackModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fPandizzleWeightFileName(ps.get< std::string> ("PandizzleWeightFileName")),
    fPandizzleAlg(ps),
    fPandizzleReader("")
{
  fPandizzleReader.AddVariable("PFPMichelNHits",&fTMVAPFPMichelNHits);
  fPandizzleReader.AddVariable("PFPMichelElectronMVA",&fTMVAPFPMichelElectronMVA);
  fPandizzleReader.AddVariable("PFPMichelRecoEnergyPlane2",&fTMVAPFPMichelRecoEnergyPlane2);
  fPandizzleReader.AddVariable("PFPTrackDeflecAngleSD",&fTMVAPFPTrackDeflecAngleSD);
  fPandizzleReader.AddVariable("PFPTrackLength",&fTMVAPFPTrackLength);
  fPandizzleReader.AddVariable("PFPTrackEvalRatio",&fTMVAPFPTrackEvalRatio);
  fPandizzleReader.AddVariable("PFPTrackConcentration",&fTMVAPFPTrackConcentration);
  fPandizzleReader.AddVariable("PFPTrackCoreHaloRatio",&fTMVAPFPTrackCoreHaloRatio);
  fPandizzleReader.AddVariable("PFPTrackConicalness",&fTMVAPFPTrackConicalness);
  fPandizzleReader.AddVariable("PFPTrackdEdxStart",&fTMVAPFPTrackdEdxStart);
  fPandizzleReader.AddVariable("PFPTrackdEdxEnd",&fTMVAPFPTrackdEdxEnd);
  fPandizzleReader.AddVariable("PFPTrackdEdxEndRatio",&fTMVAPFPTrackdEdxEndRatio);
  std::string weight_file_name = fPandizzleWeightFileName;
  std::string weight_file_path;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(weight_file_name, weight_file_path);
  fPandizzleReader.BookMVA("BDTG",weight_file_path);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Track> FDSelectionTools::CheatRecoPfoTrackSelector::SelectTrack(art::Event const & evt){

  //std::cout << "fFittedTrackModuleLabel: " << fFittedTrackModuleLabel << std::endl;
  //std::cout << "fAllTrackModuleLabel: " << fAllTrackModuleLabel << std::endl;

  bool foundSignal(false);
  art::Ptr<recob::Track> selSignalTrack, selTrack;

  // For later - to obtain track hits
  art::Handle< std::vector<recob::Track> > trackListHandle_forced;
  std::vector<art::Ptr<recob::Track> > trackList_forced;
  if (evt.getByLabel(fAllTrackModuleLabel, trackListHandle_forced))
    art::fill_ptr_vector(trackList_forced, trackListHandle_forced);

  // Get PFParticles from event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  // Get the generator record (to find whether event is FHC/RHC)
  art::Handle< std::vector<simb::MCTruth> > mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle))
    art::fill_ptr_vector(mcList, mcTruthListHandle);

  if (mcList.size() != 1){
      std::cout << "CCNuSelection: There are " << mcList.size() << " MCTruth in this event. No track selected" << std::endl;
      return selTrack;
  }

  const bool isFHC = (mcList[0]->GetNeutrino().Nu().PdgCode() > 0);
  const bool isNC = mcList[0]->GetNeutrino().CCNC();
  const bool isNumu = (std::abs(mcList[0]->GetNeutrino().Nu().PdgCode()) == 14);
  const bool isCCNumu = (isNumu && !isNC);
  //std::cout << "isFHC? " << (isFHC ? "yep" : "no") << std::endl;

  // Search for the highest pandizzle signal and global tracks
  double highestPandizzleScore(std::numeric_limits<double>::lowest()), highestSignalPandizzleScore(std::numeric_limits<double>::lowest());

  for (art::Ptr<recob::PFParticle> pfParticle : pfparticleList)
  {
    // Get the forced track fit 
    art::FindManyP<recob::Track> fmtpfp_forced(pfparticleListHandle, evt, fAllTrackModuleLabel);
    const std::vector<art::Ptr<recob::Track> > pfp_track_vector_forced = fmtpfp_forced.at(pfParticle.key());

    // Get the associated track - may not exist
    art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fFittedTrackModuleLabel);
    const std::vector<art::Ptr<recob::Track> > pfp_track_vector = fmtpfp.at(pfParticle.key());

    if ((pfp_track_vector_forced.size() > 1) || (pfp_track_vector.size() > 1))
    { 
      //Found a PFP with more than one track matched.  Complain and exit
      std::cout<< "Number of associated tracks to a PFP is greater than 1" << std::endl;
      continue;
    } 

    if (pfp_track_vector_forced.size() == 0)
    {
        // Maybe a track couldn't be fitted? Who knows the ways of larsoft
        continue;
    }
  
    const bool isFitted(pfp_track_vector.size() != 0);
    const art::Ptr<recob::Track> track_forced = pfp_track_vector_forced[0];

    // Get hits to perform truth matching
    art::FindManyP<recob::Hit> fmht_forced(trackListHandle_forced, evt, fAllTrackModuleLabel);
    const std::vector<art::Ptr<recob::Hit> > trackHits_forced = fmht_forced.at(track_forced.key());

    // Find true PDG of reco track
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, trackHits_forced, 0);

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

    // Get Pandizzle score
    fPandizzleAlg.ProcessPFParticle(pfParticle, evt);

    fTMVAPFPMichelNHits = (float)(fPandizzleAlg.GetIntVar("PFPMichelNHits"));
    fTMVAPFPMichelElectronMVA = fPandizzleAlg.GetFloatVar("PFPMichelElectronMVA");
    fTMVAPFPMichelRecoEnergyPlane2 = fPandizzleAlg.GetFloatVar("PFPMichelRecoEnergyPlane2");
    fTMVAPFPTrackDeflecAngleSD = fPandizzleAlg.GetFloatVar("PFPTrackDeflecAngleSD");
    fTMVAPFPTrackLength = fPandizzleAlg.GetFloatVar("PFPTrackLength");
    fTMVAPFPTrackEvalRatio = fPandizzleAlg.GetFloatVar("PFPTrackEvalRatio");
    fTMVAPFPTrackConcentration = fPandizzleAlg.GetFloatVar("PFPTrackConcentration");
    fTMVAPFPTrackCoreHaloRatio = fPandizzleAlg.GetFloatVar("PFPTrackCoreHaloRatio");
    fTMVAPFPTrackConicalness = fPandizzleAlg.GetFloatVar("PFPTrackConicalness");
    fTMVAPFPTrackdEdxStart = fPandizzleAlg.GetFloatVar("PFPTrackdEdxStart");
    fTMVAPFPTrackdEdxEnd = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEnd");
    fTMVAPFPTrackdEdxEndRatio = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEndRatio");
    fTMVAPFPTrackPIDA = fPandizzleAlg.GetFloatVar("PFPTrackPIDA");

    /*
  std::cout << "PFPMichelNHits: " << fTMVAPFPMichelNHits << std::endl;
  std::cout << "PFPMichelElectronMVA: " << fTMVAPFPMichelElectronMVA << std::endl;
  std::cout << "PFPMichelRecoEnergyPlane2: " << fTMVAPFPMichelRecoEnergyPlane2 << std::endl;
  std::cout << "PFPTrackDeflecAngleSD: " << fTMVAPFPTrackDeflecAngleSD << std::endl;
  std::cout << "PFPTrackLength: " <<fTMVAPFPTrackLength << std::endl;
  std::cout << "PFPTrackEvalRatio: " <<fTMVAPFPTrackEvalRatio << std::endl;
  std::cout << "PFPTrackConcentration: " << fTMVAPFPTrackConcentration << std::endl;
  std::cout << "PFPTrackCoreHaloRatio: " << fTMVAPFPTrackCoreHaloRatio << std::endl;
  std::cout << "PFPTrackConicalness: " << fTMVAPFPTrackConicalness << std::endl;
  std::cout << "PFPTrackdEdxStart: " << fTMVAPFPTrackdEdxStart << std::endl;
  std::cout << "PFPTrackdEdxEnd: " << fTMVAPFPTrackdEdxEnd << std::endl;
  std::cout << "PFPTrackdEdxEndRatio: " << fTMVAPFPTrackdEdxEndRatio << std::endl;
    */

    const double pandizzleScore = fPandizzleReader.EvaluateMVA("BDTG");

    //std::cout << "pandizzleScore: " << pandizzleScore << std::endl;

    const bool isForcedPrimaryMuon = ((isFHC && (matchedPDG == 13) && primary) || (!isFHC && (matchedPDG == -13) && primary));
  
    if (isForcedPrimaryMuon)
    {
        //std::cout << "found a primary muon" << std::endl;
        //std::cout << "is moun fitted? " << (isFitted ? "yes" : "no") << std::endl;

        if (pandizzleScore > highestSignalPandizzleScore)
        {
            foundSignal = true;
            highestSignalPandizzleScore = pandizzleScore;
            selSignalTrack = track_forced;
        }
    }
    else if (isFitted)
    {
        if (pandizzleScore > highestPandizzleScore)
        {
            highestPandizzleScore = pandizzleScore;
            selTrack = track_forced;
        }
    }
  }

  /*
  std::cout << "found? " << (foundSignal ? "yep" : "no") << std::endl;
  std::cout << "highestSignalPandizzleScore: " << highestSignalPandizzleScore << std::endl;
  std::cout << "highestPandizzleScore: " << highestPandizzleScore << std::endl;
  */

  return (isCCNumu ? (foundSignal ? selSignalTrack : selTrack) : selTrack);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


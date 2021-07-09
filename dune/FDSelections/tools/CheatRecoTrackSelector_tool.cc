#include "CheatRecoTrackSelector.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Utils/TruthMatchUtils.h"

FDSelectionTools::CheatRecoTrackSelector::CheatRecoTrackSelector(fhicl::ParameterSet const& ps) :
    fNuGenModuleLabel(ps.get< std::string> ("ModuleLabels.NuGenModuleLabel")),
    fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
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

art::Ptr<recob::Track> FDSelectionTools::CheatRecoTrackSelector::SelectTrack(art::Event const & evt){

  bool foundSignal(false);
  art::Ptr<recob::Track> selSignalTrack, selTrack;

  // For later - to obtain track hits
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > trackList;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(trackList, trackListHandle);

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

  // Search for the highest pandizzle signal and global tracks
  double highestPandizzleScore(std::numeric_limits<double>::lowest()), highestSignalPandizzleScore(std::numeric_limits<double>::lowest());

  for (art::Ptr<recob::PFParticle> pfParticle : pfparticleList)
  {
    // Get the associated track 
    art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fTrackModuleLabel);
    const std::vector<art::Ptr<recob::Track> > pfp_track_vector = fmtpfp.at(pfParticle.key());

    if (pfp_track_vector.size() > 1)
    { 
      // Found a PFP with more than one track matched.  Complain and move on
      std::cout<< "Number of associated tracks to a PFP is greater than 1: " << pfp_track_vector.size() << std::endl;
      continue;
    } 
    else if (pfp_track_vector.size() == 0)
    { 
      // Don't bother complaining of no track was found.  It's either missing or its a shower.  No biggie here
      continue; 
    }

    const art::Ptr<recob::Track> track = pfp_track_vector[0];

    // Get hits to perform truth matching
    art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
    const std::vector<art::Ptr<recob::Hit> > current_track_hits = fmht.at(track.key());

    // Find true PDG of reco track
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, current_track_hits, 0);

    int matchedPDG = -999;
    bool primary = false;

    if (TruthMatchUtils::Valid(g4id)){
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
        matchedPDG = matched_mcparticle->PdgCode();
        primary = (matched_mcparticle->Process() == "primary");
    }

    bool isPrimaryMuon = ((isFHC && (matchedPDG == 13) && primary) || (!isFHC && (matchedPDG == -13) && primary));

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
    const double pandizzleScore = fPandizzleReader.EvaluateMVA("BDTG");
    /*
    std::cout << "fTMVAPFPMichelNHits: " << fTMVAPFPMichelNHits << std::endl;
    std::cout << "fTMVAPFPMichelElectronMVA: " << fTMVAPFPMichelElectronMVA << std::endl;
    std::cout << "fTMVAPFPMichelRecoEnergyPlane2: " << fTMVAPFPMichelRecoEnergyPlane2 << std::endl;
    std::cout << "fTMVAPFPTrackDeflecAngleSD: " << fTMVAPFPTrackDeflecAngleSD << std::endl;
    std::cout << "fTMVAPFPTrackLength" << fTMVAPFPTrackLength << std::endl;
    std::cout << "fTMVAPFPTrackEvalRatio: " << fTMVAPFPTrackEvalRatio << std::endl;
    std::cout << "fTMVAPFPTrackConcentration: " << fTMVAPFPTrackConcentration << std::endl;
    std::cout << "fTMVAPFPTrackCoreHaloRatio: " << fTMVAPFPTrackCoreHaloRatio << std::endl;
    std::cout << "fTMVAPFPTrackConicalness: " << fTMVAPFPTrackConicalness << std::endl;
    std::cout << "fTMVAPFPTrackdEdxStart: " << fTMVAPFPTrackdEdxStart << std::endl;
    std::cout << "fTMVAPFPTrackdEdxEnd: " << fTMVAPFPTrackdEdxEnd << std::endl;
    std::cout << "fTMVAPFPTrackdEdxEndRatio: " << fTMVAPFPTrackdEdxEndRatio << std::endl;
    std::cout << "fTMVAPFPTrackPIDA: " << fTMVAPFPTrackPIDA << std::endl;
    std::cout << "pandizzle score: " << pandizzleScore << std::endl;
    */

    if (isPrimaryMuon && (pandizzleScore > highestSignalPandizzleScore))
    {
        foundSignal = true;
        highestSignalPandizzleScore = pandizzleScore;
        selSignalTrack = track;
    }
        
    if (pandizzleScore > highestPandizzleScore)
    {
        highestPandizzleScore = pandizzleScore;
        selTrack = track;
    }
  }

  //std::cout << "found? " << (foundSignal ? "yep" : "no") << std::endl;
  //std::cout << "highestSignalPandizzleScore: " << highestSignalPandizzleScore << std::endl;
  //std::cout << "highestPandizzleScore: " << highestPandizzleScore << std::endl;

  return (isCCNumu ? (foundSignal ? selSignalTrack : selTrack) : selTrack);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


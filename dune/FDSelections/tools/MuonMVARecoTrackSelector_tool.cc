#include "MuonMVARecoTrackSelector.h"


art::Ptr<recob::Track> FDSelectionTools::MuonMVARecoTrackSelector::SelectTrack(art::Event const & evt){
  art::Ptr<recob::Track> mytrack;

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > trackList;
  if (evt.getByLabel(fTrackModuleLabel, trackListHandle)){
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  art::FindManyP<anab::MVAPIDResult> fmpidt(trackListHandle, evt, fPIDModuleLabel);
  double max_MVA = -999;
  //Loop over the tracks to get the longest one NICE
  for (unsigned int i_track = 0; i_track < trackList.size(); i_track++){
    //Get the muon MVA value for this particular track
    std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(trackList[i_track].key());
    std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
    double muon_MVA = mvaOutMap["muon"]; 
    if (muon_MVA > max_MVA){
      max_MVA = muon_MVA;
      mytrack = trackList[i_track]; 
    }
  }
  return mytrack;
}

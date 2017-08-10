#include "LongestRecoTrackSelector.h"


art::Ptr<recob::Track> FDSelectionTools::LongestRecoTrackSelector::SelectTrack(art::Event const & evt){
  art::Ptr<recob::Track> mytrack;

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > trackList;
  if (evt.getByLabel(fTrackModuleLabel, trackListHandle)){
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  double longest_track_length = -999;
  //Loop over the tracks to get the longest one NICE
  for (unsigned int i_track = 0; i_track < trackList.size(); i_track++){
    double current_track_length = trackList[i_track]->Length();
    if (current_track_length > longest_track_length){
      longest_track_length = current_track_length;
      mytrack = trackList[i_track];
    }
  }

  return mytrack;
}

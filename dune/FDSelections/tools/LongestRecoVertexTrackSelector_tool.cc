#include "LongestRecoVertexTrackSelector.h"


art::Ptr<recob::Track> FDSelectionTools::LongestRecoVertexTrackSelector::SelectTrack(art::Event const & evt){
  art::Ptr<recob::Track> mytrack;

  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)){
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);
  }
  for (unsigned int i_pfp = 0; i_pfp <pfparticleList.size(); i_pfp++){
    art::Ptr<recob::PFParticle> pfparticle = pfparticleList[i_pfp];
    int pdgcode = pfparticle->PdgCode();
    if ((std::abs(pdgcode)==14 || std::abs(pdgcode)==12)){
      art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fTrackModuleLabel);
      const std::vector<art::Ptr<recob::Track> > primary_tracks = fmtpfp.at(pfparticle.key());

      double longest_track_length = -999;
      //Loop over the tracks to get the longest one NICE
      for (unsigned int i_track = 0; i_track < primary_tracks.size(); i_track++){
        double current_track_length = primary_tracks[i_track]->Length();
        if (current_track_length > longest_track_length){
          longest_track_length = current_track_length;
          mytrack = primary_tracks[i_track];
        }
      }
    }
    break; //no need to keep looping
  }


  return mytrack;
}

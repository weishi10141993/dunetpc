#include "HighestEnergyRecoShowerSelector.h"


art::Ptr<recob::Shower> FDSelectionTools::HighestEnergyRecoShowerSelector::SelectShower(art::Event const & evt){
  art::Ptr<recob::Shower> myshower;

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerList;
  if (evt.getByLabel(fShowerModuleLabel, showerListHandle)){
    art::fill_ptr_vector(showerList, showerListHandle);
  }
  else {
    return myshower;
  }

  //art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);

  // Get the highest energy shower
  int i_energyist_shower = -1;
  double energy_energyist_shower = -999.;
  art::ServiceHandle<geo::Geometry> geom;

  for (std::vector<art::Ptr<recob::Shower> >::const_iterator showerIt = showerList.begin(); showerIt != showerList.end(); ++showerIt) {
    //const std::vector<art::Ptr<recob::Hit> > showerHits = fmhs.at(showerIt->key());
    //std::map<int,double> showerEnergy;
    for (unsigned int plane = 0; plane < geom->MaxPlanes(); ++plane)
      showerEnergy[plane] = fShowerEnergyAlg.ShowerEnergy(showerHits, plane);
    int best_plane = -1;
    double highest_energy_plane = 0;
    for (std::map<int,double>::const_iterator showerEnergyIt = showerEnergy.begin(); showerEnergyIt != showerEnergy.end(); ++showerEnergyIt) {
      if (showerEnergyIt->second > highest_energy_plane) {
	highest_energy_plane = showerEnergyIt->second;
	best_plane = showerEnergyIt->first;
      }
    }
    if (best_plane < 0)
      return myshower;
    double current_energy = showerEnergy.at(best_plane);
    if (current_energy > energy_energyist_shower) {
      energy_energyist_shower = current_energy;
      i_energyist_shower = showerIt->key();
    }
  }
  if (i_energyist_shower < 0)
    return myshower;

  myshower = showerList.at(i_energyist_shower); 
  return myshower;
}

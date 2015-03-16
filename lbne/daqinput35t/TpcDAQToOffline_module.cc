////////////////////////////////////////////////////////////////////////
// Class:       TpcDAQToOffline
// Module Type: producer
// File:        TpcDAQToOffline_module.cc
//
// Generated at Mon Sep  1 10:00:30 2014 by Jonathan Davies using artmod
// from cetpkgsupport v1_06_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>

//lbne-artdaq includes
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragments.hh"

//larsoft includes
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "Geometry/Geometry.h"

#include "tpcFragmentToRawDigits.h"

#include "utilities/UnpackFragment.h"

namespace DAQToOffline {
  class TpcDAQToOffline;
}

class DAQToOffline::TpcDAQToOffline : public art::EDProducer {
public:
  explicit TpcDAQToOffline(fhicl::ParameterSet const & pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TpcDAQToOffline(TpcDAQToOffline const &) = delete;
  TpcDAQToOffline(TpcDAQToOffline &&) = delete;
  TpcDAQToOffline & operator = (TpcDAQToOffline const &) = delete;
  TpcDAQToOffline & operator = (TpcDAQToOffline &&) = delete;
  void produce(art::Event & evt) override;
  void reconfigure(const fhicl::ParameterSet &pset);
  void printParameterSet();

private:

  std::string fFragType;
  std::string fRawDataLabel;
  std::string fOutputDataLabel;
  bool fDebug;
  raw::Compress_t        fCompression;      ///< compression type to use
  unsigned int           fZeroThreshold;    ///< Zero suppression threshold


};


DAQToOffline::TpcDAQToOffline::TpcDAQToOffline(fhicl::ParameterSet const & pset)
{

  this->reconfigure(pset);

  produces< std::vector<raw::RawDigit> > (fOutputDataLabel);  

}

void DAQToOffline::TpcDAQToOffline::reconfigure(fhicl::ParameterSet const& pset){

  fFragType = pset.get<std::string>("FragType");
  fRawDataLabel = pset.get<std::string>("RawDataLabel");
  fOutputDataLabel = pset.get<std::string>("OutputDataLabel");
  fDebug = pset.get<bool>("Debug");

  fZeroThreshold=0;
  fCompression=raw::kNone;
  if(fDebug) printParameterSet();

}

void DAQToOffline::TpcDAQToOffline::printParameterSet(){

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "Parameter Set" << std::endl;
  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;

  std::cout << "fFragType: " << fFragType << std::endl;
  std::cout << "fRawDataLabel: " << fRawDataLabel << std::endl;
  std::cout << "fOutputDataLabel: " << fOutputDataLabel << std::endl;
  std::cout << "fDebug: ";
  if(fDebug) std::cout << "true" << std::endl;
  else std::cout << "false" << std::endl;

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;    


}

void DAQToOffline::TpcDAQToOffline::produce(art::Event & evt)
{
  art::Handle<artdaq::Fragments> rawFragments;
  evt.getByLabel(fRawDataLabel, fFragType, rawFragments);

  art::EventNumber_t eventNumber = evt.event();

  //Check that the data is valid
  if(!rawFragments.isValid()){
    std::cerr << "Run: " << evt.run()
	      << ", SubRun: " << evt.subRun()
	      << ", Event: " << eventNumber
	      << " is NOT VALID" << std::endl;
    throw cet::exception("rawFragments NOT VALID");
  }

  // Check if there is RCE data in this event
  // Don't crash code if not present, just don't save anything
  try { rawFragments->size(); }
  catch(std::exception e) {
    std::cout << "WARNING: Raw RCE data not found in event " << eventNumber << std::endl;
    return;
  }

  auto digits = tpcFragmentToRawDigits(evt.id(), *rawFragments, fDebug, fCompression, fZeroThreshold);

  evt.put(std::make_unique<decltype(digits)>(std::move(digits)), fOutputDataLabel);
}

DEFINE_ART_MODULE(DAQToOffline::TpcDAQToOffline)



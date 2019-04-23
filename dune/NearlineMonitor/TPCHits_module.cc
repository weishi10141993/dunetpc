////////////////////////////////////////////////////////////////////////
// Class:       TPCHits
// Module Type: analyzer
// File:        TPCHits_module.cc
//
// Generated at Wed Nov 18 16:12:15 2015 by Jonathan Davies using artmod
// from cetpkgsupport v1_08_07.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"


//lbne-artdaq includes
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragment.hh"


//larsoft
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"

//ROOT
#include "TH1I.h"
#include "TH2I.h"
#include "TMath.h"

#include <sstream>

namespace nearline {
  class TPCHits;
}

class nearline::TPCHits : public art::EDAnalyzer {
public:
  explicit TPCHits(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCHits(TPCHits const &) = delete;
  TPCHits(TPCHits &&) = delete;
  TPCHits & operator = (TPCHits const &) = delete;
  TPCHits & operator = (TPCHits &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void endJob();
  void makeHistograms();
  void printRawDigitTimingInfo(art::Event const & e);
  void printTpcFragTimingInfo(art::Event const & e);


private:
  art::InputTag fTPCHitTag;
  art::InputTag fRawDigitsTag;
  bool fPrintRawDigitTimingInfo;

  //Variables for comparing previous to next event
  art::EventID fEventIDLast;
  art::Timestamp fTimestampLast;
  art::EventNumber_t fEventLast;

  bool fPrintTpcFragTimingInfo;


  std::vector<unsigned int> fTotalNumHitsPerTPC;
  std::vector<TH1I*> fHistNumHitsPerTPC;
  std::vector<TH1I*> fHistNumTicksPerEvent;
  std::vector<TH1I*> fHistNumTicksPerEventPair;
  std::vector<TH2I*> fHistNumTicksPerEventVsLast;


  std::vector<unsigned int> fNumChansTPC;
  std::vector<unsigned int> fNumTicksTPC;
  std::vector<unsigned int> fNumChansTPCLastEvent;
  std::vector<unsigned int> fNumTicksTPCLastEvent;
  std::vector<unsigned int> fNumChansTPCTwoEvents;
  std::vector<unsigned int> fNumTicksTPCTwoEvents;


};


nearline::TPCHits::TPCHits(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  reconfigure(p);

  art::ServiceHandle<geo::Geometry> geom;
  fTotalNumHitsPerTPC = std::vector<unsigned int>(geom->NTPC());
  fNumChansTPC = std::vector<unsigned int>(geom->NTPC());
  fNumTicksTPC = std::vector<unsigned int>(geom->NTPC());
  fNumChansTPCLastEvent = std::vector<unsigned int>(geom->NTPC());
  fNumTicksTPCLastEvent = std::vector<unsigned int>(geom->NTPC());
  fNumChansTPCTwoEvents = std::vector<unsigned int>(geom->NTPC());
  fNumTicksTPCTwoEvents = std::vector<unsigned int>(geom->NTPC());


  makeHistograms();

}

void nearline::TPCHits::makeHistograms(){

  mf::LogInfo loginfo("TPCHits");
  loginfo << "===================================="   << "\n"
          << "Making Histograms"                      << "\n"
          << "===================================="   << "\n";
  
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<geo::Geometry> geom;
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  const double samplingRate = detprop->SamplingRate();
  const unsinged int millisliceSize = 5e6;
  const unsigned int maxNumTicksPerEvent = millisliceSize / samplingRate;

  fHistNumHitsPerTPC = std::vector<TH1I*> (geom->NTPC());
  fHistNumTicksPerEvent = std::vector<TH1I*>(geom->NTPC());
  fHistNumTicksPerEventPair = std::vector<TH1I*>(geom->NTPC());
  fHistNumTicksPerEventVsLast = std::vector<TH2I*>(geom->NTPC());  

  //Currently 5ms of data per pair of events
  //Sampling rate of the TPC is 2MHz (0.5 ns spacing) -> can get this from the det properties
  //That means that 5000 / 0.5 = 10,000 ticks per 2 events
  //Max number of hits is therefore == 10,000 * nwires per TPC
  
  loginfo << "Making Number of Hits per TPC Histograms\n";

  for(size_t tpc = 0; tpc < fHistNumHitsPerTPC.size(); tpc++){
    std::ostringstream histName;
    histName << "hist_num_hits_tpc_" << tpc;
    
    std::ostringstream histTitle;
    histTitle << "Number of hits per event TPC " << tpc << ";Log_{10}(Number of hits);Number of events";

    const unsigned int numPlanes = geom->Nplanes(tpc);
    unsigned int totalNumWires = 0;
    for(unsigned int plane=0;plane<numPlanes;plane++){
      unsigned int numWires = geom->Nwires(plane, tpc);
      totalNumWires += numWires;
    }//plane
    
    unsigned int maxNumHits = totalNumWires * maxNumTicksPerEvent;
    
    int numBins = 100;
    int binLow = -1;
    int binHigh = TMath::Log10(maxNumHits);
    TH1I* tempHist = tfs->make<TH1I>(histName.str().c_str(), histTitle.str().c_str(), numBins, binLow, binHigh);
    fHistNumHitsPerTPC.at(tpc) = tempHist;
    loginfo << "histName: "  << histName.str()      << "\n";
  }

  loginfo << "\n"
          << "Making Number of ticks / wire / event / TPC Histograms\n";

  for(size_t tpc = 0; tpc < fHistNumTicksPerEvent.size(); tpc++){
    std::ostringstream histName;
    histName << "hist_num_ticks_per_event_tpc_" << tpc;
    
    std::ostringstream histTitle;
    histTitle << "Number of ticks / wire / event in TPC " << tpc << ";Number of ticks;Number of events";
    
    int numBins = 100;
    int binLow = 0;
    int binHigh = maxNumTicksPerEvent;
    TH1I* tempHist = tfs->make<TH1I>(histName.str().c_str(), histTitle.str().c_str(), numBins, binLow, binHigh);
    fHistNumTicksPerEvent.at(tpc) = tempHist;
    loginfo << "histName: "  << histName.str()      << "\n";
  }//tpc

  loginfo << "\n"
          << "Making Number of ticks / wire / event pair / TPC Histograms\n";

  for(size_t tpc = 0; tpc < fHistNumTicksPerEventPair.size(); tpc++){
    std::ostringstream histName;
    histName << "hist_num_ticks_per_event_pair_tpc_" << tpc;
    
    std::ostringstream histTitle;
    histTitle << "Number of ticks / wire / event pair in TPC " << tpc << ";Number of ticks;Number of event pairs";
    
    int numBins = 100;
    int binLow = 0;
    int binHigh = maxNumTicksPerEvent;
    TH1I* tempHist = tfs->make<TH1I>(histName.str().c_str(), histTitle.str().c_str(), numBins, binLow, binHigh);
    fHistNumTicksPerEventPair.at(tpc) = tempHist;
    loginfo << "histName: "  << histName.str()      << "\n";
  }//tpc

  loginfo << "\n"
          << "Making Number of ticks / wire / event vs last / TPC Histograms\n";

  for(size_t tpc = 0; tpc < fHistNumTicksPerEventVsLast.size(); tpc++){
    std::ostringstream histName;
    histName << "hist_num_ticks_per_event_vs_last_tpc_" << tpc;
    
    std::ostringstream histTitle;
    histTitle << "Number of ticks / wire / event vs last event in TPC " << tpc << ";Number of ticks;Number of ticks (last)";
    
    int numBins = 100;
    int binLow = 0;
    int binHigh = maxNumTicksPerEvent;
    TH2I* tempHist = tfs->make<TH2I>(histName.str().c_str(), histTitle.str().c_str(), numBins, binLow, binHigh, numBins, binLow, binHigh);
    fHistNumTicksPerEventVsLast.at(tpc) = tempHist;
    loginfo << "histName: "  << histName.str()      << "\n";
  }//tpc

  loginfo << "===================================="   << "\n"
          << "Finished Making Histograms"             << "\n"
          << "===================================="   << "\n";
}

void nearline::TPCHits::reconfigure(fhicl::ParameterSet const & p ){
  
  fTPCHitTag = p.get<art::InputTag>("TPCHitTag", "a:b:c");
  fRawDigitsTag = p.get<art::InputTag>("RawDigitsTag", "a:b:c");
  fPrintRawDigitTimingInfo = p.get<bool>("fPrintRawDigitTimingInfo", true);
  fPrintTpcFragTimingInfo = p.get<bool>("fPrintTpcFragTimingInfo", true);

  mf::LogInfo("TPCHits") << "===================================="   << "\n"
                         << "Parameter Set"                          << "\n"
                         << "===================================="   << "\n"
                         << "fTPCHitTag:          " << fTPCHitTag            << "\n"
                         << "fRawDigitsTag:       " << fRawDigitsTag         << "\n"
                         << "fPrintRawDigitTimingInfo: " << fPrintRawDigitTimingInfo   << "\n"
                         << "fPrintTpcFragTimingInfo: " << fPrintTpcFragTimingInfo   << "\n"
                         << "===================================="   << "\n";
}

void nearline::TPCHits::printTpcFragTimingInfo(art::Event const & e){

  mf::LogInfo("TPCHits::printTpcFragTimingInfo") << "Hello World" << std::endl;

  art::Handle<artdaq::Fragments> rawFragments;
  art::InputTag tpcFragTag("daq:TPC");
  bool retVal =  e.getByLabel("daq:TPC", rawFragments);
  if(retVal==true) 
    mf::LogInfo("TPCHits::printTpcFragTimingInfo") << "Getting TPC Frag SUCCESS: " << tpcFragTag << std::endl;
  else{ 
    mf::LogWarning("TPCHits::printTpcFragTimingInfo") << "Getting TPC Frag FAIL: " << tpcFragTag << std::endl;
    return;
  }
  try { rawFragments->size(); }
  catch(std::exception e) {
    mf::LogError("TPCHits::printTpcFragTimingInfo") << "WARNING: Issue with rawFragments for TPC hits" << std::endl;
    return;
  }

  for(size_t fragIndex = 0; fragIndex < rawFragments->size(); fragIndex++){
    const artdaq::Fragment &singleFragment = rawFragments->at(fragIndex);
    lbne::TpcMilliSliceFragment msf(singleFragment);
    auto nMicroSlices = msf.microSliceCount();
    
    mf::LogInfo("TPCHits::printTpcFragTimingInfo") << "fragIndex: " << fragIndex 
                                                   << " fragmentID: " << singleFragment.fragmentID() 
                                                   << " ms counter: " << nMicroSlices
                                                   << std::endl;

  }//fragIndex
}



void nearline::TPCHits::printRawDigitTimingInfo(art::Event const & e){

  art::EventID eventID = e.id();
  art::Timestamp timestamp = e.time();
  art::EventNumber_t event = e.event();

  mf::LogInfo("TPCHits::printRawDigitTimingInfo") << "EventID: " << eventID
          << " time: " << timestamp.value()
          << " EventNum:  " << event
          << std::endl;
  
  if((timestamp.value() - fTimestampLast.value() < 0) || (event - fEventLast < 0)){
    mf::LogInfo("TPCHits::printRawDigitTimingInfo") << "ERROR - negative delta timestamp of eventnum" << std::endl;
  }

  fEventIDLast = eventID;
  fTimestampLast = timestamp;
  fEventLast = event;
}


void nearline::TPCHits::analyze(art::Event const & e)
{

  if(fPrintTpcFragTimingInfo) printTpcFragTimingInfo(e);
  if(fPrintRawDigitTimingInfo) printRawDigitTimingInfo(e);
  return;//FIXME

  art::ServiceHandle<geo::Geometry> geom;
  std::vector<unsigned int> thisNumHits(geom->NTPC());
  art::Handle<std::vector<recob::Hit> > hitHandle;
  art::Handle<std::vector<raw::RawDigit> > digitHandle;
                            
  bool retVal = e.getByLabel(fTPCHitTag, hitHandle);
  if(retVal==true) 
    ;
  //mf::LogInfo("TPCHits") << "Getting Hits SUCCESS: " << fTPCHitTag << std::endl;
  else{ 
    mf::LogWarning("TPCHits") << "Getting Hits FAIL: " << fTPCHitTag << std::endl;
    return;
  }
  
  try { hitHandle->size(); }
  catch(std::exception e) {
    mf::LogError("TPCHits") << "WARNING: Issue with hitHandle for TPC hits" << std::endl;
    return;
  }

  if(!hitHandle.isValid()){
    mf::LogError("TPCHits") << "Run: " << e.run()
                            << ", SubRun: " << e.subRun()
                            << ", Event: " << e.event()
                            << " is NOT VALID";
    throw cet::exception("hits NOT VALID");
    return;
  }


  retVal = e.getByLabel(fRawDigitsTag, digitHandle);
  if(retVal==true) 
    ;
  //mf::LogInfo("TPCHits") << "Getting RawDigits SUCCESS: " << fRawDigitsTag << std::endl;
  else{ 
    mf::LogWarning("TPCHits") << "Getting RawDigits FAIL: " << fRawDigitsTag << std::endl;
    return;
  }
  
  try { digitHandle->size(); }
  catch(std::exception e) {
    mf::LogError("TPCHits") << "WARNING: Issue with digitHandle for RawDigits" << std::endl;
    return;
  }

  if(!digitHandle.isValid()){
    mf::LogError("TPCHits") << "Run: " << e.run()
                            << ", SubRun: " << e.subRun()
                            << ", Event: " << e.event()
                            << " is NOT VALID";
    throw cet::exception("RawDigit NOT VALID");
    return;
  }

  std::vector<art::Ptr<recob::Hit> > hitlist;
  art::fill_ptr_vector(hitlist, hitHandle);
 
  //  size_t numHits = hitlist.size();
  //  mf::LogInfo("TPCHits") << "NumHits: " << numHits << " hits" << std::endl;
  for (size_t tpchit_index = 0; tpchit_index<hitlist.size(); ++tpchit_index){
    geo::WireID wireid = hitlist[tpchit_index]->WireID();
    thisNumHits.at(wireid.TPC)++;
    fTotalNumHitsPerTPC.at(wireid.TPC)++;
  }//hit

  for(size_t tpc=0; tpc < fHistNumHitsPerTPC.size();tpc++){
    TH1I* tempHist = fHistNumHitsPerTPC.at(tpc);
    if(thisNumHits.at(tpc)==0)
      tempHist->Fill(-1);
    else
      tempHist->Fill(TMath::Log10(thisNumHits.at(tpc)));
  }//tpc
  
  size_t numDigitChans = digitHandle->size();
  //  mf::LogInfo("TPCHits") << "numDigitChans: " << numDigitChans << std::endl;
  //  size_t numChans = geom->Nchannels();
  //  mf::LogInfo("TPCHits") << "numChans: " << numChans << std::endl;

  fNumChansTPC = std::vector<unsigned int> (geom->NTPC());
  fNumTicksTPC = std::vector<unsigned int> (geom->NTPC());
  for(size_t rdIter=0;rdIter<numDigitChans;rdIter++){
    art::Ptr<raw::RawDigit> digitVec(digitHandle, rdIter);
    auto channel =  digitVec->Channel();    
    auto numSamples = digitVec->Samples();
    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);

    if(wids.at(0).Plane !=2 ) continue;
    fNumChansTPC.at(wids.at(0).TPC)++;
    fNumTicksTPC.at(wids.at(0).TPC)+=numSamples;
    fNumChansTPCTwoEvents.at(wids.at(0).TPC)++;
    fNumTicksTPCTwoEvents.at(wids.at(0).TPC)+=numSamples;
    

  }//rdIter
  
  // mf::LogError loginfo("TPCHits");
  // loginfo << "Numbers of channels\n";
  // for(size_t tpc=0;tpc<fNumChansTPC.size();tpc++){
  //   if(tpc!=5) continue;
  //   loginfo << "tpc: " << tpc 
  //           << " numChannels: " << fNumChansTPC.at(tpc) 
  //           << " numSamples: " << fNumTicksTPC.at(tpc)
  //           << " numSamplesTwoEvents: " << fNumTicksTPCTwoEvents.at(tpc)
  //           << "\n";
    
  // }
  
  for(size_t tpc=0;tpc<fHistNumTicksPerEvent.size();tpc++){

    unsigned int numChans = fNumChansTPC.at(tpc);
    unsigned int numTicks = fNumTicksTPC.at(tpc);
    unsigned int numTicksLastEvent = fNumTicksTPCLastEvent.at(tpc);
    //    unsigned int numTicksEventPair = fNumTicksTPCTwoEvents.at(tpc);

    // if(numTicks>0){
    //   mf::LogInfo loginfo("TPCHits");
    //   loginfo << "tpc: " << tpc << "\n"
    //           << " numTicks: " << numTicks << "\n"
    //           << " numTicksEventPair: " << numTicksEventPair << "\n";
    // }

    if(numChans==0) continue;
    TH1I* tempHist = fHistNumTicksPerEvent.at(tpc);    
    tempHist->Fill(numTicks/numChans);

    tempHist = fHistNumTicksPerEventPair.at(tpc);    
    //    tempHist->Fill(numTicksEventPair/numChans);
    tempHist->Fill((numTicksLastEvent+numTicks)/numChans);

    TH2I* tempHist2D = fHistNumTicksPerEventVsLast.at(tpc);
    tempHist2D->Fill(numTicks/numChans, numTicksLastEvent/numChans);

  }//tpc

  fNumChansTPCLastEvent = fNumChansTPC;
  fNumTicksTPCLastEvent = fNumTicksTPC;
  fNumChansTPCTwoEvents = fNumChansTPC;
  fNumTicksTPCTwoEvents = fNumTicksTPC;
  
}

void nearline::TPCHits::endJob(){

  std::ostringstream os;

  os << "===================================="   << "\n"
     << "End Job"                                << "\n"
     << "===================================="   << "\n";
  
  for(size_t i=0; i<fTotalNumHitsPerTPC.size(); i++) os << "i: " << i << " contents: " << fTotalNumHitsPerTPC.at(i) << std::endl;

  os << "===================================="   << "\n";
  mf::LogInfo("TPCHits") << os.str();
  
  
}

DEFINE_ART_MODULE(nearline::TPCHits)

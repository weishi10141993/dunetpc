#ifndef BEAMDATA_PROTODUNEBEAMEVENT_H
#define BEAMDATA_PROTODUNEBEAMEVENT_H

#include <vector>
#include <bitset>
#include <iostream>
#include <map>
#include "lardataobj/RecoBase/Track.h"
#include "dune/DuneObj/ProtoDUNEBeamSpill.h"

namespace beam
{
  

  class ProtoDUNEBeamEvent{
    public:
      ProtoDUNEBeamEvent();
     ~ProtoDUNEBeamEvent(){};

      const std::pair< double, double >  & GetT0() const{ return t0;};
      const double & GetT0Sec()  const { return t0.first; };
      const double & GetT0Nano() const { return t0.second; };
      void                                 SetT0(std::pair< double, double > theT0){ t0 = theT0;};
     
      const FBM &       GetFBM(std::string) const;
      void              SetFBMTrigger(std::string, FBM); 

      void                       DecodeFibers(std::string);
      const short              & GetFiberStatus(std::string, size_t) const;
      const std::vector<short> & GetActiveFibers(std::string) const;
      const double             & GetFiberTime(std::string) const; 

      std::bitset<32>    toBinary(double num){return std::bitset<32>( (uint32_t(num)) );}; 

      void                       SetCKov0(CKov theCKov){ CKov0 = theCKov; }; 
      void                       SetCKov1(CKov theCKov){ CKov1 = theCKov; };
      const short  &             GetCKov0Status() const{ return CKov0.trigger; };
      const short  &             GetCKov1Status() const{ return CKov1.trigger; };
      const double &             GetCKov0Time() const{ return CKov0.timeStamp; };
      const double &             GetCKov1Time() const{ return CKov1.timeStamp; };
      const double &             GetCKov0Pressure() const{ return CKov0.pressure; };
      const double &             GetCKov1Pressure() const{ return CKov1.pressure; };


      void SetTOF0Trigger( std::pair<double,double> theT){ TOF0 = theT;};
      void SetTOF1Trigger( std::pair<double,double> theT){ TOF1 = theT;}; 
      void SetTOFChan    ( int theChan )                 { TOFChan = theChan; };

      const std::pair< double, double > & GetTOF0() const { return TOF0; };
      const double & GetTOF0Sec()  const { return TOF0.first; };
      const double & GetTOF0Nano() const { return TOF0.second; };

      const std::pair< double, double > & GetTOF1() const { return TOF1; };
      const double & GetTOF1Sec()  const { return TOF1.first; };
      const double & GetTOF1Nano() const { return TOF1.second; };

     
      void DecodeTOF(){ theTOF = ( (TOF1.first  - TOF0.first ) * 1.e9 
                                 + (TOF1.second - TOF0.second) ); };

      const double &            GetTOF() const { return theTOF; };
      const int &               GetTOFChan() const{ return TOFChan; };

      void                                AddBeamTrack(recob::Track theTrack){ Tracks.push_back(theTrack);};
      const recob::Track &                GetBeamTrack(size_t i) const{ return Tracks.at(i);};
      size_t                              GetNBeamTracks() {return Tracks.size();}
      const std::vector< recob::Track > & GetBeamTracks() const;
      void                                ClearBeamTracks(){ Tracks.clear(); };

      void              SetActiveTrigger(size_t theTrigger){ activeTrigger = theTrigger; isMatched = true; };
      const bool &      CheckIsMatched()const { return isMatched; };
      void              SetUnmatched(){ isMatched = false; };
      const size_t &    GetActiveTrigger() const{ return activeTrigger; };

      void                          AddRecoBeamMomentum( double theMomentum ){ RecoBeamMomenta.push_back( theMomentum ); };
      const std::vector< double > & GetRecoBeamMomenta() const               { return RecoBeamMomenta;};
      size_t                        GetNRecoBeamMomenta()                    { return RecoBeamMomenta.size(); };
      const double &                GetRecoBeamMomentum( size_t i ) const    { return RecoBeamMomenta.at(i); };  
      void                          ClearRecoBeamMomenta()                   { RecoBeamMomenta.clear(); }; 

      void              SetBITrigger(int theTrigger){ BITrigger = theTrigger; };
      const int &       GetBITrigger() const{ return BITrigger; };

      void              SetTimingTrigger(int theTrigger){ TimingTrigger = theTrigger; };
      const int &       GetTimingTrigger() const{ return TimingTrigger; };
      
      void              SetSpillStart(double theSpillStart){ SpillStart = theSpillStart; };
      const double &    GetSpillStart() const{ return SpillStart; };

      void              SetSpillOffset(double theSpillOffset){ SpillOffset = theSpillOffset; };
      const double &    GetSpillOffset() const{ return SpillOffset; };

      void              SetCTBTimestamp(long long theCTBTimestamp){ CTBTimestamp = theCTBTimestamp; };
      const long long & GetCTBTimestamp() const{ return CTBTimestamp; };

      void              SetRDTimestamp(long long theRDTimestamp){ RDTimestamp = theRDTimestamp; };
      const long long & GetRDTimestamp() const{ return RDTimestamp; };



    private:

      //Time of a coincidence between 2 TOFs
      //Signalling a good particle
      //Known as 'GeneralTrigger'
      //
      std::pair<double,double> t0;

      //Timestamp from the CTB signaling a 
      //Good particle signal was received
      //
      long long CTBTimestamp;
      long long RDTimestamp;

      //Set of FBMs
      //Indices: [Monitor in beam]
      std::map<std::string, FBM > fiberMonitors;

      //Set of TOF detectors
      //
      std::pair< double, double > TOF0;
      std::pair< double, double > TOF1;
      int TOFChan;
      double theTOF;

      //Set of Cerenkov detectors
      //
      CKov CKov0;
      CKov CKov1;

      std::vector<recob::Track> Tracks;
  
      size_t activeTrigger;
      bool   isMatched = false;

      std::vector< double > RecoBeamMomenta;

      int BITrigger;
      int TimingTrigger;
      double SpillStart;
      double SpillOffset;

  };

  inline ProtoDUNEBeamEvent::ProtoDUNEBeamEvent(){ 
    FBM dummy; 
    dummy.ID = -1;
    fiberMonitors["dummy"] = dummy;
  }

  inline const std::vector< recob::Track > & ProtoDUNEBeamEvent::GetBeamTracks() const { return Tracks; }
  

  ////////////Fiber Monitor Access
  inline const FBM & ProtoDUNEBeamEvent::GetFBM(std::string FBMName) const{
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
    
      std::cout << "Error FBM " << FBMName << " not found" << std::endl;
        
      for(auto itF = fiberMonitors.begin(); itF != fiberMonitors.end(); ++itF){
        std::cout << "\t" << itF->first << std::endl; 
      }

      return fiberMonitors.at("dummy");
    }

    return fiberMonitors.at(FBMName); 
  }

  inline void ProtoDUNEBeamEvent::SetFBMTrigger(std::string FBMName, FBM theFBM){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM not found" << std::endl;
      std::cout << "Setting New FBM " << FBMName << std::endl;
    }

    fiberMonitors[FBMName] = theFBM;
  }



  inline void ProtoDUNEBeamEvent::DecodeFibers(std::string FBMName){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " Not Found"  << std::endl;
      return;
    }
    //This always clears the currently active fibers in the FBM.
    fiberMonitors[FBMName].active.clear();
    
    for(int iSet = 0; iSet < 6; ++iSet){
      
      std::bitset<32> theseFibers = toBinary( fiberMonitors[FBMName].fiberData[iSet] );

      for(int  iFiber = 0; iFiber < 32; ++iFiber){      
        fiberMonitors[FBMName].fibers[iSet*32 + iFiber] = theseFibers[iFiber];
        if(theseFibers[iFiber]) fiberMonitors[FBMName].active.push_back(iSet*32 + iFiber);
      }
    }

    fiberMonitors[FBMName].decoded = true;
  }


  inline const short & ProtoDUNEBeamEvent::GetFiberStatus(std::string FBMName, size_t iFiber) const{
/*    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return -1;          
    }
    if( (iFiber > 191)){
      std::cout << "Please input fiber in range [0,191]" << std::endl;
      return -1;
    }
*/
    return fiberMonitors.at(FBMName).fibers.at(iFiber);
  }

  inline const std::vector<short> & ProtoDUNEBeamEvent::GetActiveFibers(std::string FBMName) const{
//    std::vector<short> active;

    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return fiberMonitors.at("dummy").active;          
    }
    
/*    for(size_t iF = 0; iF < 192; ++iF){
      if(fiberMonitors[FBMName].fibers[iF]) active.push_back(iF); 
    }
*/
    return fiberMonitors.at(FBMName).active;
  }

  inline const double & ProtoDUNEBeamEvent::GetFiberTime(std::string FBMName) const{
/*    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return -1;          
    }
*/
    return fiberMonitors.at(FBMName).timeStamp;
  }

  /////////////////////////////////


}


#endif

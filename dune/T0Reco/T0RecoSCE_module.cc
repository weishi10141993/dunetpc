////////////////////////////////////////////////////////////////////////
// Class:       T0RecoSCE
// Module Type: analyser
// File:        T0RecoSCE_module.cc
//
// Joshua Thompson - joshualt@fnal.gov
// developed from work by Hannah Rogers   - hannah.rogers@colostate.edu
// based on uboonecode modules by David Caratelli and Chris Barnes
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

// services etc...
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"

// ROOT
#include "TVector3.h"
#include <TTree.h>

// C++
#include <memory>
#include <iostream>
#include <utility>

class T0RecoSCE;

class T0RecoSCE : public art::EDAnalyzer {
public:
	explicit T0RecoSCE(fhicl::ParameterSet const & p);
	// The destructor generated by the compiler is fine for classes
	// without bare pointers or other resource use.
	
	// Plugins should not be copied or assigned.
	T0RecoSCE(T0RecoSCE const &) = delete;
	T0RecoSCE(T0RecoSCE &&) = delete;
	T0RecoSCE & operator = (T0RecoSCE const &) = delete;
	T0RecoSCE & operator = (T0RecoSCE &&) = delete;
	
	void beginJob() override;
	
	//Required functions.
	void analyze(art::Event const & evt) override;
	
private:
	// Functions to be used throughout module
  
	void   SortTrackPoints	(const recob::Track& track, std::vector<TVector3>& sorted_trk);

	void   SplitTrack	(const recob::Track& track, std::vector<TVector3>& sorted_trk);
	
	size_t FlashMatch	(const double reco_time, std::vector<double> op_times_v);
		
	// Declare member data here

	std::string fTrackProducer;
	std::string fFlashProducer;
	std::string fHitProducer;
	std::string	fTriggerProducer;
	std::string	fPFPProducer;

	std::string	fOpHitProducer;

	bool		fUseMC;

	double		fEdgeWidth; // [cm]
	
	bool 		fCathode;
	bool 		fAnode;

	bool 		fAllFlash;

	double	 	fMinPE;
	double 		fMinTrackLength;

	double 		fFlashScaleFactor;
	double 		fFlashTPCOffset;

	bool 		fUseOpHits;

	int 		fFirstOpChannel;
	int 		fLastOpChannel;

	bool 		fDebug;	

	double		fDriftVelocity; // [cm/us]
	unsigned int 	fReadoutWindow;

	bool		fAnodeT0Check;
	std::string fAnodeT0Producer;

	double		det_top;
	double 		det_bottom;
	double 		det_front;
	double 		det_back;
	double 		det_width; // [cm]
	
	std::vector<double> 	op_times;

	std::vector<size_t> 	flash_id_v;
	art::Handle<std::vector<recob::OpFlash> > 	flash_h;

	std::vector<size_t> 	op_id_v;
	art::Handle<std::vector<recob::OpHit> > 	op_hit_h;

	TTree *track_tree;
	// Track parameters
	double 	rc_time;
	double 	length;
	double 	rc_xs, rc_xe, rc_xs_corr, rc_xe_corr;
	double 	rc_ys, rc_ye;
	double 	rc_zs, rc_ze;

	// Track parameters for cathode crossing anode piercing tracks
	double	anode_rc_time;
	double	cathode_rc_time;
	double	dt_anode_cathode;

	// Flash parameters
	double 	matched_flash_pe;	
	double 	matched_flash_time;
	double 	matched_flash_time_width;
	double 	corrected_matched_flash_time;
	double	matched_flash_centre_y;
	double	matched_flash_centre_z;
	double	matched_flash_max_pe_det_x;
	double	matched_flash_max_pe_det_y;
	double	matched_flash_max_pe_det_z;
	double	matched_flash_width_y;
	double	matched_flash_width_z;

	// Reco results
	double 	dt_flash_reco;

	// Track info
	bool 	readout_edge;
	bool 	TPC_entering_candidate;
	bool 	TPC_exiting_candidate;
	bool 	anode_piercing_candidate;
	bool 	cathode_crossing_track;
	int 	sister_track;

	// MC info
	double 	mc_time;
	double 	mc_particle_ts;
	double 	mc_particle_te;
	double	mc_particle_x_anode;
	double	mc_particle_y_anode;
	double	mc_particle_z_anode;
	double 	dt_mc_reco;
	bool	true_anode_piercer;

	// Tree for flash info
	TTree 	*flash_tree;
	double 	flash_reco_time_diff;
	double 	flash_time;
	double 	flash_time_width;
	double 	corrected_flash_time;
	double 	flash_pe;
	double	flash_centre_y;
	double	flash_centre_z;
	double	flash_width_y;
	double	flash_width_z;

	// Tree for event info
	TTree 	*ev_tree;
	int 	run;
	int		event;
	int 	total_particle_ctr;
	int 	ev_ctr;
	
	};

T0RecoSCE::T0RecoSCE(fhicl::ParameterSet const & fcl)
	:
	EDAnalyzer(fcl) {

	fTrackProducer     	= fcl.get<std::string>	("TrackProducer"    	);
	fHitProducer		= fcl.get<std::string>	("HitProducer"      	);
	fFlashProducer     	= fcl.get<std::string>	("FlashProducer"    	);
	fTriggerProducer	= fcl.get<std::string>	("TriggerProducer"  	);
	fPFPProducer		= fcl.get<std::string>	("PFPProducer"  	);

	fOpHitProducer		= fcl.get<std::string>	("OpHitProducer"	);

	fUseMC            	= fcl.get<bool>			("UseMC"            	);

	fEdgeWidth			= fcl.get<double>		("EdgeWidth"       	);

	fMinPE				= fcl.get<double> 		("MinPE"           	);
	fMinTrackLength		= fcl.get<double>		("MinTrackLength"    	);

	fCathode			= fcl.get<bool>  		("CathodeCrossers"  	);
	fAnode				= fcl.get<bool>  		("AnodePiercers"    	);

	fDebug				= fcl.get<bool>  		("Debug"            	);

	fAllFlash			= fcl.get<bool>  		("AllFlashToTrackTimeDiffs");

	fFlashScaleFactor	= fcl.get<double>		("FlashScaleFactor" 	);
	fFlashTPCOffset		= fcl.get<double>		("FlashTPCOffset"	);
	
	fUseOpHits			= fcl.get<bool>			("UseOpHits"		);
	fFirstOpChannel		= fcl.get<int>			("FirstOpChannel"	);
	fLastOpChannel		= fcl.get<int>			("LastOpChannel"	);

	fAnodeT0Check		= fcl.get<bool>			("CheckAssignedAnodeT0"	);
	fAnodeT0Producer	= fcl.get<std::string>	("AnodeT0Producer"		);


	// get boundaries based on detector bounds
	auto const* geom = lar::providerFrom<geo::Geometry>();
  
	det_top = fEdgeWidth;
	det_bottom = fEdgeWidth;
	det_front = fEdgeWidth;
	det_back = fEdgeWidth;

	for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
		geo::TPCGeo const& TPC = geom->TPC(tID);
   
		if(TPC.DriftDistance() < 25.0) continue;
   
		double origin[3] = {0.};
		double center[3] = {0.};
		TPC.LocalToWorld(origin, center);
   
		double tpc_top = center[1] + TPC.HalfHeight();
		double tpc_bottom = center[1] - TPC.HalfHeight();
		double tpc_front = center[2] - TPC.HalfLength();
		double tpc_back = center[2] + TPC.HalfLength();
   
		if (tpc_top 	> det_top) 	det_top = tpc_top;
		if (tpc_bottom 	< det_bottom) 	det_bottom = tpc_bottom;
		if (tpc_front 	< det_front) 	det_front = tpc_front;
		if (tpc_back 	> det_back) 	det_back = tpc_back;  
   
		det_width = TPC.DriftDistance();
		}

  	// Use 'detp' to find 'efield' and 'temp'
	auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
	double efield = detp -> Efield();
	std::cout << "Nominal electric field is: " << efield*1000 << " V/cm" << std::endl;

	double temp   = detp -> Temperature();
	std::cout << "LAr temperature is: " << temp << " K" << std::endl;

	// Determine the drift velocity from 'efield' and 'temp'
	fDriftVelocity = detp -> DriftVelocity(efield,temp);
	std::cout << "Drift velocity is: " << fDriftVelocity << " cm/us" << std::endl;

	// Get Readout window length
	
	fReadoutWindow = detp->ReadOutWindowSize();
	std::cout << "Readout window is: " << fReadoutWindow << " ticks" << std::endl;
	}  
  
void T0RecoSCE::beginJob(){

	art::ServiceHandle<art::TFileService> tfs;

	track_tree = tfs->make<TTree>("track_tree","SCE calibrations variables");
	//Track parameters
	track_tree->Branch("rc_time",&rc_time,"rc_time/D");
	track_tree->Branch("length", &length, "length/D");
	track_tree->Branch("rc_xs",&rc_xs,"rc_xs/D");
	track_tree->Branch("rc_xs_corr",&rc_xs_corr,"rc_xs/D");
	track_tree->Branch("rc_ys",&rc_ys,"rc_ys/D");
	track_tree->Branch("rc_zs",&rc_zs,"rc_zs/D");
	track_tree->Branch("rc_xe",&rc_xe,"rc_xe/D");
	track_tree->Branch("rc_xe_corr",&rc_xe_corr,"rc_xe_corr/D");
	track_tree->Branch("rc_ye",&rc_ye,"rc_ye/D");
	track_tree->Branch("rc_ze",&rc_ze,"rc_ze/D");

	if(fCathode&&fAnode) {
		track_tree->Branch("cathode_rc_time",&cathode_rc_time,"cathode_rc_time/D");
		track_tree->Branch("anode_rc_time",&anode_rc_time,"anode_rc_time/D");
		track_tree->Branch("dt_anode_cathode",&dt_anode_cathode,"dt_anode_cathode/D");
		}

	// Information on whether track crosses anode or cathode
	track_tree->Branch("TPC_entering_candidate",&TPC_entering_candidate,"TPC_entering_candidate/B");
	track_tree->Branch("TPC_exiting_candidate",&TPC_exiting_candidate,"TPC_exiting_candidate/B");
	track_tree->Branch("cathode_crossing_track",&cathode_crossing_track,"cathode_crossing_track/B");
	
	// Flash parameters
	track_tree->Branch("matched_flash_time",&matched_flash_time,"matched_flash_time/D");
	track_tree->Branch("matched_flash_time_width",&matched_flash_time_width,"matched_flash_time_width/D");
	//track_tree->Branch("corrected_matched_flash_time",&corrected_matched_flash_time,"corrected_matched_flash_time/D");
	track_tree->Branch("matched_flash_pe",&matched_flash_pe,"matched_flash_pe/D");
	track_tree->Branch("matched_flash_centre_y",&matched_flash_centre_y,"matched_flash_centre_y/D");
	track_tree->Branch("matched_flash_centre_z",&matched_flash_centre_z,"matched_flash_centre_z/D");
	track_tree->Branch("matched_flash_max_pe_det_x",&matched_flash_max_pe_det_x,"matched_flash_max_pe_det_x/D");
	track_tree->Branch("matched_flash_max_pe_det_y",&matched_flash_max_pe_det_y,"matched_flash_max_pe_det_y/D");
	track_tree->Branch("matched_flash_max_pe_det_z",&matched_flash_max_pe_det_z,"matched_flash_max_pe_det_z/D");
	track_tree->Branch("matched_flash_width_y",&matched_flash_width_y,"matched_flash_width_y/D");
	track_tree->Branch("matched_flash_width_z",&matched_flash_width_z,"matched_flash_width_z/D");

	// Flash -> track time difference
	track_tree->Branch("dt_flash_reco",&dt_flash_reco,"dt_flash_reco/D");
	
	// Branches for MC truth info
	if(fUseMC) {
		track_tree->Branch("mc_time",&mc_time,"mc_time/D");
		track_tree->Branch("dt_mc_reco",&dt_mc_reco,"dt_mc_reco/D");

		track_tree->Branch("mc_particle_x_anode",&mc_particle_x_anode,"mc_particle_x_anode/D");
		track_tree->Branch("mc_particle_y_anode",&mc_particle_y_anode,"mc_particle_y_anode/D");
		track_tree->Branch("mc_particle_z_anode",&mc_particle_z_anode,"mc_particle_z_anode/D");
		track_tree->Branch("true_anode_piercer",&true_anode_piercer,"true_anode_piercer/B");
		}

	// Flash Tree
	if (fAllFlash) {
		flash_tree = tfs->make<TTree>("flash_tree","Flash properties and reco time differences");
		flash_tree->Branch("flash_time",&flash_time,"flash_time/D");
		flash_tree->Branch("flash_time_width",&flash_time_width,"flash_time_width/D");
		flash_tree->Branch("corrected_flash_time",&corrected_flash_time,"corrected_flash_time/D");
		flash_tree->Branch("flash_reco_time_diff",&flash_reco_time_diff,"flash_reco_time_diff/D");
		flash_tree->Branch("flash_pe",&flash_pe,"flash_pe/D");
		flash_tree->Branch("flash_centre_y",&flash_centre_y,"flash_centre_y/D");
		flash_tree->Branch("flash_centre_z",&flash_centre_z,"flash_centre_z/D");
		flash_tree->Branch("flash_width_y",&flash_width_y,"flash_width_y/D");
		flash_tree->Branch("flash_width_z",&flash_width_z,"flash_width_z/D");
		}

	// Event Tree
	ev_tree = tfs->make<TTree>("ev_tree","Event information");
	ev_tree->Branch("total_particle_ctr",&total_particle_ctr,"total_particle_ctr/I");
	ev_tree->Branch("ev_ctr",&ev_ctr,"ev_ctr/I");
	ev_tree->Branch("run",&run,"run/I");
	ev_tree->Branch("event",&event,"event/I");
	ev_ctr = 0;
	total_particle_ctr = 0;
	
	}
  
void T0RecoSCE::analyze(art::Event const & evt){
	
	event = evt.event();
	run = evt.run();
 	int track_number = 0;
	ev_ctr++;

	// Load detector clocks for later
	auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();

	double TPC_trigger_offset = 0.0;

	std::vector<std::vector<TLorentzVector>> mcpart_list;
	
	std::cout << "Event number: " << ev_ctr << std::endl;
	if(fDebug) std::cout << "Set event number: " << event << "\nTop: " << det_top 
	<< "\nBottom: " << det_bottom << "\nFront: " << det_front << "\nBack: " << det_back 
	<< "\nEdge width: " << fEdgeWidth << std::endl;  
  
	op_times.clear();
	flash_id_v.clear();
	flash_h.clear();

	op_id_v.clear();
	op_hit_h.clear();
	
	// load Flashes
	if (fDebug) std::cout << "Loading flash from producer " << fFlashProducer << std::endl;

	if(fAnode||fAllFlash){
		evt.getByLabel(fFlashProducer,flash_h);

		// make sure flashes look good
		if(!flash_h.isValid()) {
			std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    		throw std::exception();
  			}
  		}

    if(fUseOpHits) evt.getByLabel(fOpHitProducer, op_hit_h);

  	art::Handle<std::vector<recob::OpFlash> > trigger_h;
	double trigger_time = 0;

  	if(!fUseMC){
		if(fDebug) std::cout << "Loading trigger time from producer " 
		<< fTriggerProducer << std::endl;
  		evt.getByLabel(fTriggerProducer, trigger_h);
  		trigger_time = trigger_h->at(0).Time();
  	}

	// load PFParticles

	auto reco_particles_h = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFPProducer);
	
	if(!reco_particles_h.isValid()) {
		std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate PFParticles!"<<std::endl;
		throw std::exception(); 
		}

	// Utilities for PFParticles and tracks
	protoana::ProtoDUNEPFParticleUtils pfpUtil;
	protoana::ProtoDUNETrackUtils trackUtil;
	protoana::ProtoDUNETruthUtils truthUtil;

	// load MCParticles

  	art::Handle<std::vector<simb::MCParticle> > mcpart_h;
  	evt.getByLabel("largeant",mcpart_h);

  	// if we should use MCParticle
  	if (fUseMC){
		// make sure particles exist
		if(!mcpart_h.isValid()) {
			std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCParticle!"<<std::endl;
			throw std::exception(); 
			}
   	
		// ^ if use MCParticle
		}

	// Get trigger to TPC Offset
		TPC_trigger_offset = detclock->TriggerOffsetTPC();
		if(fDebug) std::cout << "TPC time offset from trigger: " 
			<< TPC_trigger_offset << " us" << std::endl;

	// Prepare a vector of optical flash times, if flash above some PE cut value

	if((fAnode||fAllFlash)&&!fUseOpHits) {
  		size_t flash_ctr = 0;
  		for (auto const& flash : *flash_h){
    		if (flash.TotalPE() > fMinPE){
				double op_flash_time = flash.Time() - trigger_time;
				if(fUseMC) op_flash_time = op_flash_time - TPC_trigger_offset;
      			op_times.push_back(op_flash_time);
      			flash_id_v.push_back(flash_ctr);
      			//if (fDebug) std::cout << "\t Flash: " << flash_ctr 
					//<< " has time : " << flash.Time() - trigger_time 
					//<< ", PE : " << flash.TotalPE() << std::endl;
   				}
    		flash_ctr++;
  			} // for all flashes

  		if(!fAnode) {
			auto const& output_flash = FlashMatch(0.0,op_times);
			if (fDebug) std::cout << "Output all flashes - closest flash to trigger time is: "
				<< output_flash << std::endl; }

   		if(fDebug) std::cout << "Selected a total of " << op_times.size() << " OpFlashes/OpHits" << std::endl;
		}

	if(fUseOpHits) {
		size_t op_ctr = 0;
		for (auto const& op_hit : *op_hit_h){
			int op_hit_channel = op_hit.OpChannel();
			if(op_hit_channel>=fFirstOpChannel&&op_hit_channel<=fLastOpChannel) {
				double op_hit_time = op_hit.PeakTime() - trigger_time;
				//if(fUseMC) op_hit_time = op_hit_time - TPC_trigger_offset;
				op_times.push_back(op_hit_time);
				op_id_v.push_back(op_ctr);
				//if(fDebug) std::cout << "OpHit channel: " << op_hit_channel << std::endl;
				}
			op_ctr++;
			}
		}
 
	// LOOP THROUGH RECONSTRUCTED PFPARTICLES

	size_t ev_particle_ctr = 0;

	for(unsigned int p = 0; p < reco_particles_h->size(); ++p){

		recob::PFParticle particle = (*reco_particles_h)[p];

    	// Only consider primary particles
    	if(!particle.IsPrimary()) continue;

		ev_particle_ctr++;
		total_particle_ctr++;

		const recob::Track* track = pfpUtil.GetPFParticleTrack(particle,evt,fPFPProducer,fTrackProducer);
		if(track == 0x0) { 
			if(fDebug) std::cout << "\tPFParticle " << ev_particle_ctr << " is not track like" << std::endl;
			continue; } 

		if (fDebug) std::cout << "\tLooping through reco PFParticle " << ev_particle_ctr << std::endl;  

		const std::vector<const recob::Hit*>& hit_v = trackUtil.GetRecoTrackHits(*track,evt,fTrackProducer);

		rc_time = fReadoutWindow;
		anode_rc_time = fReadoutWindow;
		cathode_rc_time = -fReadoutWindow;
		dt_anode_cathode = -999.;

		length = 0.;
		rc_xs = fReadoutWindow/10.;
		rc_xe = -fReadoutWindow/10.; 
		rc_xs_corr = 399.;
		rc_xe_corr = -399.;
		rc_ys = -99.;
		rc_ye = -99.; 
		rc_zs = -99.; 
		rc_ze = -99.;

		matched_flash_time = -fReadoutWindow;
		corrected_matched_flash_time = -fReadoutWindow;
		matched_flash_time_width = -1.;
		matched_flash_pe = 0.;
		matched_flash_centre_y = -99.;
		matched_flash_centre_z = -99.;
		matched_flash_width_y = -99.;
		matched_flash_width_z = -99.;
		matched_flash_max_pe_det_x = 0.;
		matched_flash_max_pe_det_y = -99.;
		matched_flash_max_pe_det_z = -99.;

		dt_flash_reco = 999;

		mc_time = 2*fReadoutWindow;

		anode_piercing_candidate = false;
		cathode_crossing_track = false;
		sister_track = 0;
  
 		// Get sorted points for the track object [assuming downwards going]
    	std::vector<TVector3> sorted_trk;
	    SortTrackPoints(*track,sorted_trk);

    	TVector3 track_start = sorted_trk.at(0);
   		TVector3 track_end = sorted_trk.at(sorted_trk.size() - 1);

	    if(fDebug) std::cout << "\t\tTrack goes from (" << track_start.X() << ", " 
			<< track_start.Y() << ", " << track_start.Z() << ") --> (" << 
			track_end.X() << ", " << track_end.Y() << ", " << track_end.Z() << 
			")" << std::endl;

		if(sqrt(pow(track_start.X() - track_end.X(),2.0) 
			+ pow(track_start.Y() - track_end.Y(),2.0)
			+ pow(track_start.Z() - track_end.Z(),2.0)) < fMinTrackLength )	{
      		if(fDebug) std::cout << "\t\t\tParticle track too short. Skipping." << std::endl;
      		continue;
	      	}

     	// Determine if the track crosses the cathode 
		auto const* geom = lar::providerFrom<geo::Geometry>();   
    	auto const* hit = hit_v.at(0);
	    const geo::WireID wireID = hit->WireID();
		const auto TPCGeoObject = geom->TPC(wireID.TPC,wireID.Cryostat);
		short int driftDir_start = TPCGeoObject.DetectDriftDirection();
		short int driftDir_end = 0;

	    for (size_t ii = 1; ii < hit_v.size(); ii++) {
    		const geo::WireID wireID2 = hit_v.at(ii)->WireID();
			const auto TPCGeoObject2 = geom->TPC(wireID2.TPC,wireID2.Cryostat);
			driftDir_end = TPCGeoObject2.DetectDriftDirection(); 
		
			if(driftDir_end + driftDir_start == 0){
				cathode_crossing_track = true;
				ii = hit_v.size();
				}
			}

		const simb::MCParticle* mc_particle = 0x0;
		int last_mc_point = 0;

		// if we should use MC info -> continue w/ MC validation
      	if (fUseMC){

			mc_particle = truthUtil.GetMCParticleFromRecoTrack(*track,evt,fTrackProducer);

			if(mc_particle==0x0) { 
				if(fDebug) std::cout << "\t\t\tNo MC particle matched to PFParticle " 
				<< ev_particle_ctr << std::endl;
				continue; 
				}

			mc_particle_ts = mc_particle->T(0);

			last_mc_point = mc_particle->NumberTrajectoryPoints()-1;
			mc_particle_te = mc_particle->T(last_mc_point)/1000;

			mc_time = detclock->G4ToElecTime(mc_particle_ts) + TPC_trigger_offset;

			if(fDebug&&fabs(mc_particle_te - mc_particle_ts)>1) std::cout << 
			"\t\t\tMC Particle end time: " << mc_particle_te << 
			" us is significantly different from start time: " << mc_time << 
			" us" << std::endl;
			}

		// -------------------------------------------------------------------------------
		// CATHODE CROSSERS

		if(fCathode){

			if(cathode_crossing_track) {
				//GET T0 FROM PFPARTICLE
				auto t0_v = pfpUtil.GetPFParticleT0(particle,evt,fPFPProducer);
				if(t0_v.size() == 0) continue; 
				auto t0 = t0_v.at(0);
				cathode_rc_time = t0.Time()/1000; 

				if (fDebug) std::cout << "\t\tTrack crossers cathode - PFParticle has t0: " 
					<< cathode_rc_time << " against MC particle t0: " << mc_time << std::endl;
			
				rc_time = cathode_rc_time;

				std::vector<TVector3> split_trk = sorted_trk;

    			SplitTrack(*track,split_trk);
    			std::vector<TVector3> top_trk = {split_trk.at(0), split_trk.at(1)};
	    		std::vector<TVector3> bottom_trk = {split_trk.at(2), split_trk.at(3)};
    	
		    	if(fDebug) std::cout << "\t\tCathode-crossing point: (" 
					<< top_trk.at(1).X() << ", " << top_trk.at(1).Y() << ", " 
					<< top_trk.at(1).Z() << ") --> (" << bottom_trk.at(0).X() 
					<< ", " << bottom_trk.at(0).Y() << ", " << bottom_trk.at(0).Z() 
					<< ")" << std::endl;
    	
    			// Top Track
		    	//if(fDebug) std::cout << "\tTop track" << std::endl;
    	
    			sister_track = track_number + 1;
		    	TVector3 top_track_start = top_trk.at(0);
	    		TVector3 top_track_end = top_trk.at(1);
		
				rc_xs = top_track_start.X();
				rc_xs_corr = rc_xs;
    			rc_ys = top_track_start.Y();
    			rc_zs = top_track_start.Z();
    			rc_xe = top_track_end.X();
    			rc_xe_corr = rc_xe;
    			rc_ye = top_track_end.Y();
    			rc_ze = top_track_end.Z();
    	
				length = sqrt(pow(rc_xs - rc_xe,2.0) + pow(rc_ys - rc_ye,2.0) + 
					pow(rc_zs - rc_ze,2.0));
		
				track_tree->Fill();
				track_number++;
		
				// Bottom Track
				//if(fDebug) std::cout << "\tBottom track" << std::endl;
				sister_track = track_number - 1;
		
				TVector3 bottom_track_start = bottom_trk.at(0);
		    	TVector3 bottom_track_end = bottom_trk.at(1);
		
				rc_xs = bottom_track_start.X();
				rc_xs_corr = rc_xs;
    			rc_ys = bottom_track_start.Y();
    			rc_zs = bottom_track_start.Z();
    			rc_xe = bottom_track_end.X();
    			rc_xe_corr = rc_xe;
    			rc_ye = bottom_track_end.Y();
    			rc_ze = bottom_track_end.Z();
    	
				length = sqrt(pow(rc_xs - rc_xe,2.0) + pow(rc_ys - rc_ye,2.0) + 
				pow(rc_zs - rc_ze,2.0));
		
				track_tree->Fill();
				track_number++;
				}
		
			else if(fDebug) std::cout << "\t\tTrack does not cross cathode."<< std::endl;
			}
		    
		// ------------------------------------------------------------------------------------
		// ANODE PIERCERS 
		if(fAnode){	
			// if(fDebug) std::cout << "\t\tThis track starts in TPC " << wireID.TPC <<
			//" which has a drift direction of " << driftDir_start << std::endl; 

			// create root trees variables
    
   			rc_xs = track_start.X();
		  	rc_ys = track_start.Y();
			rc_zs = track_start.Z();
	    	rc_xe = track_end.X();
	    	rc_ye = track_end.Y();
    		rc_ze = track_end.Z();
    		length = track->Length();

			double x_shift = 0;
			if(cathode_crossing_track&&cathode_rc_time>-fReadoutWindow) {
				x_shift = cathode_rc_time*driftDir_start*fDriftVelocity;
				}

			// Determine if track hits edge of readout window
			// NECESSARY TO REMOVE TRACKS THAT AREN'T FINISHED WHEN WINDOW CLOSES
			// AS WELL AS TRACKS THAT WERE BEING COLLECTED BEFORE THE WINDOW OPENED

    			readout_edge = false;
    			for (auto& hits : hit_v) {
    				auto hit_tick = hits->PeakTime();
					double hit_time = detclock->TPCTick2TrigTime(hit_tick);
    				//if(fDebug) std::cout << "\t\tHit from track " << trk_ctr << 
						//" at tick: " << hit_tick << ", in TPC " 
						//<< hits->WireID().TPC << ", plane " << hits->WireID().Plane 
						//<< " and wire " << hits->WireID().Wire << std::endl;
    				if(hit_tick < 20.0 || hit_tick > (fReadoutWindow - 20.0)){
    		 			readout_edge = true;
    		 			if(fDebug) std::cout << "\tTrack hits edge of readout window. "
							"Skipping." << std::endl;
						continue;
						}
					// If track within window, get reco time from earliest hit time
					if (hit_time < anode_rc_time) anode_rc_time = hit_time;
    				}

				if(readout_edge) continue;

				TPC_entering_candidate = false;
				TPC_exiting_candidate = false;
			
    			// Tracks which may enter TPC through an APA
    			if (rc_ys < (det_top - fEdgeWidth) && rc_zs > (det_front + fEdgeWidth) 
					&& rc_zs < (det_back - fEdgeWidth)) {

      				// reconstruct track T0 w.r.t. trigger time
					if( ( rc_xs > rc_xe && driftDir_start>0 ) || 
						( rc_xs < rc_xe && driftDir_start<0 ) ) {
						TPC_entering_candidate = true;
						if(fDebug) std::cout << "\t\tTrack may enter TPC through "
						"anode. Reco t0: " << anode_rc_time << " us" << std::endl;
     					}
    				}

    			// Tracks which may exit TPC through an APA
    			if (rc_ye > (det_bottom + fEdgeWidth) && rc_ze > (det_front + fEdgeWidth) 
					&& rc_ze < (det_back - fEdgeWidth)) {

      				// reconstruct track T0 w.r.t. trigger time	
					if( ( rc_xe > rc_xs && driftDir_end>0 ) || 
						( rc_xe < rc_xs && driftDir_end<0 ) ) {	
						TPC_exiting_candidate = true;
						if(fDebug) std::cout << "\t\tTrack may exit TPC through "
							"anode. Reco t0:" << anode_rc_time << " us" << std::endl;
      					}
  					}
			
			
				if(TPC_entering_candidate||TPC_exiting_candidate) 
					anode_piercing_candidate = true;

				if(!anode_piercing_candidate) {
					if(fDebug) std::cout << "\t\tTrack does not pierce anode." << std::endl;
					continue; 
					}

				if(TPC_entering_candidate&&TPC_exiting_candidate) {
					if(fDebug) std::cout << "\t\tTrack neither enters nor exits"
					" through non-anode TPC face. No useful end point for SCE" 
					" measurement." << std::endl;
					continue; 
					}

    			rc_xs_corr = rc_xs - x_shift + driftDir_start*anode_rc_time*fDriftVelocity;
    			rc_xe_corr = rc_xe + x_shift + driftDir_end*anode_rc_time*fDriftVelocity;

    			// FLASH MATCHING
    			size_t op_match_result = FlashMatch(anode_rc_time,op_times);

				if(op_match_result==99999) {
					if(fDebug) std::cout << "Unable to match flash to track." << std::endl;
					continue;
					}

    			if(!fUseOpHits) {	
					const art::Ptr<recob::OpFlash> flash_ptr(flash_h, op_match_result);

    				matched_flash_time = flash_ptr->Time() - trigger_time;
					if(!fUseMC) corrected_matched_flash_time = fFlashScaleFactor*matched_flash_time + fFlashTPCOffset;
					if(fUseMC) corrected_matched_flash_time = fFlashScaleFactor*matched_flash_time + fFlashTPCOffset - TPC_trigger_offset;
					matched_flash_time_width = flash_ptr->TimeWidth();

	    			dt_flash_reco = corrected_matched_flash_time - anode_rc_time;
	    
	    			matched_flash_pe = flash_ptr->TotalPE();
					matched_flash_centre_y = flash_ptr->YCenter();
					matched_flash_centre_z = flash_ptr->ZCenter();
					matched_flash_width_y = flash_ptr->YWidth();
					matched_flash_width_z = flash_ptr->ZWidth();

					unsigned int max_pe_channel = 9999;
					double max_pe = 0;
					unsigned int pd_ch;
					for(pd_ch = 0; pd_ch <= geom->MaxOpChannel(); pd_ch++) {
						double channel_i_pe = flash_ptr->PE(pd_ch);
						if(channel_i_pe > max_pe)  {
							max_pe = channel_i_pe;
							max_pe_channel = pd_ch;
							}
						}
	
					if(max_pe_channel==9999||max_pe==0) continue;			
	
					if(max_pe_channel<9999) {

						matched_flash_max_pe_det_x = -det_width;

						if(max_pe_channel>143) matched_flash_max_pe_det_x = det_width;
					
						/*double max_pe_det_v[3];
						geom->OpDetGeoFromOpChannel(max_pe_channel).GetCenter(max_pe_det_v);

						matched_flash_max_pe_det_x = max_pe_det_v[0];
						matched_flash_max_pe_det_y = max_pe_det_v[1]; 
						matched_flash_max_pe_det_z = max_pe_det_v[2];
						*/

						if(fDebug) std::cout << "\t\tOpChannel " << max_pe_channel << 
							" has maximum PE, and is located at: (" << 
							matched_flash_max_pe_det_x << ", " << 
							matched_flash_max_pe_det_y << ", " << 
							matched_flash_max_pe_det_z << ")" << std::endl;
						}

   					if(fDebug) std::cout << "\t\t Matched to flash w/ index " << op_match_result << " w/ PE " 
    					<< matched_flash_pe << ", corrected time " << corrected_matched_flash_time << 
						" us vs corrected reco time " << anode_rc_time << " us" << std::endl;
					}

				if(fUseOpHits) {
					const art::Ptr<recob::OpHit> op_ptr(op_hit_h, op_match_result);

					matched_flash_time = op_ptr->PeakTime() - trigger_time;

					if(fUseMC) matched_flash_time = matched_flash_time - TPC_trigger_offset;

	
					if(fDebug) std::cout << "\t\t Matched to op hit w/ index " << op_match_result << 
						" w/ time " << matched_flash_time << " us vs corrected reco time " 
						<< anode_rc_time << " us" << std::endl;

	    			dt_flash_reco = matched_flash_time - anode_rc_time;

					}

				dt_mc_reco = 999.;
				true_anode_piercer = false;

				mc_particle_x_anode = 0;
				mc_particle_y_anode = -99;
				mc_particle_z_anode = -99;

				// if we should use MC info -> continue w/ MC validation
      			if (fUseMC){

					for(int mc_hit_i=0; mc_hit_i<last_mc_point; mc_hit_i++) {

						double mc_particle_xi = mc_particle->Position(mc_hit_i).X();
						double mc_particle_yi = mc_particle->Position(mc_hit_i).Y();
						double mc_particle_zi = mc_particle->Position(mc_hit_i).Z();

						if(abs(mc_particle_xi) > (det_width - 1) && 
							abs(mc_particle_xi) < (det_width + 1) && 
							mc_particle_yi > det_bottom && mc_particle_yi < det_top && 
							mc_particle_zi > det_front && mc_particle_zi < det_back)
							true_anode_piercer = true;

						if(abs(abs(mc_particle_xi) - det_width) < 
							abs(abs(mc_particle_x_anode) - det_width)){
							mc_particle_x_anode = mc_particle_xi;
							mc_particle_y_anode = mc_particle_yi;
							mc_particle_z_anode = mc_particle_zi;
							}
						}

					if(fDebug&&true_anode_piercer) 
					std::cout << "\t\tMC Particle matched to track has t0: " 
						<< mc_time << " us against reconstructed t0: " 
						<< anode_rc_time << " us and passes through anode at: ("
						<< mc_particle_x_anode << ", " << mc_particle_y_anode << ", "
						<< mc_particle_z_anode << ")" << std::endl;

					} // ^ if we should use MCParticles

			// For verifying anode piercer T0 assigned by producer module
			
			if(fAnodeT0Check) {
				bool HasT0 = false;

    			unsigned int pIndex = particle.Self();

    			std::vector<anab::T0> t0_apt_v;

    			const art::FindManyP<anab::T0> findParticleT0s(reco_particles_h,evt,fAnodeT0Producer);
   				
				for(unsigned int p = 0; p < findParticleT0s.at(pIndex).size(); ++p){
      				t0_apt_v.push_back((*(findParticleT0s.at(pIndex)[p])));
    				}

				if(t0_apt_v.size() != 0) HasT0 = true;

				if(HasT0) {
					auto t0_apt = t0_apt_v.at(0);

					std::cout << "PFParticle has T0: " << (t0_apt.Time()/1000) << 
					" us against anode reco time: " << anode_rc_time << 
					" us with dt_flash_reco " << dt_flash_reco << " us." << std::endl;
					}
				}

			if (!fCathode||!cathode_crossing_track) rc_time = anode_rc_time;

			if (fUseMC) dt_mc_reco = mc_time - rc_time;

			if(anode_rc_time < fReadoutWindow && cathode_rc_time > -fReadoutWindow) 
			dt_anode_cathode = anode_rc_time - cathode_rc_time;

			track_tree->Fill();
			track_number++;
			
			}
    
    	}
    
	ev_tree->Fill();
	}

void T0RecoSCE::SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_trk)
{
	sorted_trk.clear();
		
	TVector3 track_start, track_end;	
	double start_y = det_bottom - fEdgeWidth;
	double end_y = det_top + fEdgeWidth;
	
	for (size_t ii = 0; ii < track.NumberTrajectoryPoints(); ii++){
		auto const& trk_loc = track.LocationAtPoint(ii);
		
		if ((trk_loc.X() < -998.)||(trk_loc.Y() < -998.)||(trk_loc.Z() < -998)) continue;
		
		if (trk_loc.Y() < end_y){
			end_y = trk_loc.Y();
			track_end = {trk_loc.X(), trk_loc.Y(), trk_loc.Z()};
		}
		if (trk_loc.Y() > start_y){
			start_y = trk_loc.Y();
			track_start = {trk_loc.X(), trk_loc.Y(), trk_loc.Z()};
		}
	}
	
	sorted_trk.push_back(track_start);
	sorted_trk.push_back(track_end);
}

void T0RecoSCE::SplitTrack(const recob::Track& track, std::vector<TVector3>& split_trk)
{
	split_trk.clear();
	
	TVector3 track_neg, track_neg_c, track_pos_c,track_pos;
	double neg_x = 0.0;
	double pos_x = 0.0;
	double neg_c = -2.0*fEdgeWidth;
	double pos_c = 2.0*fEdgeWidth;
	
	for (size_t ii = 0; ii < track.NumberTrajectoryPoints(); ii++){
		auto const& trk_loc = track.LocationAtPoint(ii);
		
		if ((trk_loc.X() < -998.)||(trk_loc.Y() < -998.)||(trk_loc.Z() < -998)) continue;
		
		if (trk_loc.X() < neg_x){
			neg_x = trk_loc.X();
			track_neg = {trk_loc.X(), trk_loc.Y(), trk_loc.Z()};
		}
		if (trk_loc.X() > pos_x){
			pos_x = trk_loc.X();
			track_pos = {trk_loc.X(), trk_loc.Y(), trk_loc.Z()};
		}
		if ((trk_loc.X() < 0.0) && (trk_loc.X() > neg_c)){
			neg_c = trk_loc.X();
			track_neg_c = {trk_loc.X(), trk_loc.Y(), trk_loc.Z()};
		}
		if ((trk_loc.X() > 0.0) && (trk_loc.X() < pos_c)){
			pos_c = trk_loc.X();
			track_pos_c = {trk_loc.X(), trk_loc.Y(), trk_loc.Z()};
		}
	}
	
	if( track_neg.Y() > track_pos.Y()){
		split_trk.push_back(track_neg);
		split_trk.push_back(track_neg_c);
		split_trk.push_back(track_pos_c);
		split_trk.push_back(track_pos);
	} else {
		split_trk.push_back(track_pos);
		split_trk.push_back(track_pos_c);
		split_trk.push_back(track_neg_c);
		split_trk.push_back(track_neg);
	}
	
}
  
size_t  T0RecoSCE::FlashMatch(const double reco_time, std::vector<double> op_times_v)
{
  	// loop through all reco'd flash times and see if one matches
  	// the time from the track/particle
  	double dt_min = 9999999.; 
  	size_t matched_op_id = 99999;

	flash_time = fReadoutWindow;
	corrected_flash_time = fReadoutWindow;

	flash_pe = -9;
	flash_centre_y = -9;
	flash_centre_z = -9;
	flash_width_y = -9;
	flash_width_z = -9;
	flash_time_width = -9;	

  	for (size_t i=0; i < op_times_v.size(); i++){
		double op_time_i = op_times_v[i];
    	double corrected_op_time_i = op_time_i*fFlashScaleFactor + fFlashTPCOffset;
    	flash_reco_time_diff = corrected_op_time_i - reco_time;
    	if (fabs(flash_reco_time_diff) < dt_min){
      		dt_min  = fabs(flash_reco_time_diff);
      		if(!fUseOpHits) matched_op_id = flash_id_v[i];
			if(fUseOpHits) matched_op_id = op_id_v[i];
    		}

		if (fAllFlash) {
			const art::Ptr<recob::OpFlash> flash_i_ptr(flash_h, flash_id_v[i]);
    		flash_time = op_time_i;
			corrected_flash_time = corrected_op_time_i;
			flash_pe = flash_i_ptr->TotalPE();
			flash_centre_y = flash_i_ptr->YCenter();
			flash_centre_z = flash_i_ptr->ZCenter();
			flash_width_y = flash_i_ptr->YWidth();
			flash_width_z = flash_i_ptr->ZWidth();
			flash_time_width = flash_i_ptr->TimeWidth();
			flash_tree->Fill();
			}
  		}
  	return matched_op_id;
	}

DEFINE_ART_MODULE(T0RecoSCE) 
  

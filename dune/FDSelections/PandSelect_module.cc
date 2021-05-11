////////////////////////////////////////////////////////////////////////
// Class:       PandSelect
// Plugin Type: producer (art v3_05_01)
// File:        PandSelect_module.cc
//
// Generated at Mon May 10 05:56:39 2021 by Isobel Mawby using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

// ROOT
#include "TMVA/Reader.h"

// ART
#include "art/Utilities/make_tool.h" 

// LARSOFT
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/Vertex.h"

// ADDED
#include "dune/FDSelections/FDSelectionData/PandSelectParams.h"
#include "dune/FDSelections/pandizzle/PandizzleAlg.h"
#include "dune/FDSelections/pandrizzle/PandrizzleAlg.h"
#include "tools/RecoShowerSelector.h"
#include "tools/RecoTrackSelector.h"

constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

class PandSelect;

class PandSelect : public art::EDProducer {
public:
  explicit PandSelect(fhicl::ParameterSet const& p);
  PandSelect(PandSelect const&) = delete;
  PandSelect(PandSelect&&) = delete;
  PandSelect& operator=(PandSelect const&) = delete;
  PandSelect& operator=(PandSelect&&) = delete;
  void produce(art::Event& e) override;
  void beginJob() override;
  void endJob() override;

  double GetPandizzleScore(art::Event& e);
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & e);
  double GetPandrizzleScore(art::Event& e);
  void FillVertexInformation(art::Event const & e);

private:
  // Module labels
  std::string fPFParticleModuleLabel;
  std::string fShowerModuleLabel;
  std::string fTrackModuleLabel;

  // Tools
  std::unique_ptr<FDSelectionTools::RecoShowerSelector> fRecoShowerSelector;
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;

  // Algs
  FDSelection::PandizzleAlg fPandizzleAlg;
  FDSelection::PandrizzleAlg fPandrizzleAlg;

  //Pandizzle TMVAStuff
  TMVA::Reader fPandizzleReader;
  float fTMVAPFPMichelNHits;
  float fTMVAPFPMichelElectronMVA;
  float fTMVAPFPMichelRecoEnergyPlane2;
  float fTMVAPFPTrackDeflecAngleSD;
  float fTMVAPFPTrackLength;
  float fTMVAPFPTrackEvalRatio;
  float fTMVAPFPTrackConcentration;
  float fTMVAPFPTrackCoreHaloRatio;
  float fTMVAPFPTrackConicalness;
  float fTMVAPFPTrackdEdxStart;
  float fTMVAPFPTrackdEdxEnd;
  float fTMVAPFPTrackdEdxEndRatio;
  float fTMVAPFPTrackPIDA;

  float f_TMVAPFPMichelNHits;
  float f_TMVAPFPMichelElectronMVA;
  float f_TMVAPFPMichelRecoEnergyPlane2;
  float f_TMVAPFPTrackDeflecAngleSD;
  float f_TMVAPFPTrackLength;
  float f_TMVAPFPTrackEvalRatio;
  float f_TMVAPFPTrackConcentration;
  float f_TMVAPFPTrackCoreHaloRatio;
  float f_TMVAPFPTrackConicalness;
  float f_TMVAPFPTrackdEdxStart;
  float f_TMVAPFPTrackdEdxEnd;
  float f_TMVAPFPTrackdEdxEndRatio;
  float f_TMVAPFPTrackPIDA;

  // Pandrizzle Stuff
  double fRecoNuVtxX;
  double fRecoNuVtxY;
  double fRecoNuVtxZ;
};

///////////////////////////////////////////////////////////////////////////////////////////////////

PandSelect::PandSelect(fhicl::ParameterSet const& p) : 
  EDProducer{p},
  fShowerModuleLabel(p.get< std::string >("ModuleLabels.ShowerModuleLabel")),
  fTrackModuleLabel(p.get< std::string >("ModuleLabels.TrackModuleLabel")),
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(p.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))},
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(p.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fPandizzleAlg(p),
  fPandrizzleAlg(p),
  fPandizzleReader("")
{
  produces<pandselect::PandSelectParams>();

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
  std::string weight_file_name = "Pandizzle_TMVAClassification_BDTG.weights.xml";
  std::string weight_file_path;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(weight_file_name, weight_file_path);
  fPandizzleReader.BookMVA("BDTG",weight_file_path);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandSelect::produce(art::Event& e)
{
  std::unique_ptr<pandselect::PandSelectParams> outputPandSelectParams = std::make_unique<pandselect::PandSelectParams>();
  outputPandSelectParams->selTrackPandizzleScore = GetPandizzleScore(e);
  outputPandSelectParams->selShowerPandrizzleScore = GetPandrizzleScore(e);
  e.put(std::move(outputPandSelectParams));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

double PandSelect::GetPandizzleScore(art::Event& e)
{
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(e.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return kDefDoub;
  }
  //Get the selected track
  art::Ptr<recob::Track> sel_track = fRecoTrackSelector->FindSelectedTrack(e);

  //If we didn't find a selected track then what's the point?
  if (!(sel_track.isAvailable())) {
    std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no track returned from selection" << std::endl; 
    return kDefDoub;
  }

  art::Ptr<recob::PFParticle> track_pfp = GetPFParticleMatchedToTrack(sel_track, e);

  fPandizzleAlg.ProcessPFParticle(track_pfp, e);

  f_TMVAPFPMichelNHits = (float)(fPandizzleAlg.GetIntVar("PFPMichelNHits"));
  f_TMVAPFPMichelElectronMVA = fPandizzleAlg.GetFloatVar("PFPMichelElectronMVA");
  f_TMVAPFPMichelRecoEnergyPlane2 = fPandizzleAlg.GetFloatVar("PFPMichelRecoEnergyPlane2");
  f_TMVAPFPTrackDeflecAngleSD = fPandizzleAlg.GetFloatVar("PFPTrackDeflecAngleSD");
  f_TMVAPFPTrackLength = fPandizzleAlg.GetFloatVar("PFPTrackLength");
  f_TMVAPFPTrackEvalRatio = fPandizzleAlg.GetFloatVar("PFPTrackEvalRatio");
  f_TMVAPFPTrackConcentration = fPandizzleAlg.GetFloatVar("PFPTrackConcentration");
  f_TMVAPFPTrackCoreHaloRatio = fPandizzleAlg.GetFloatVar("PFPTrackCoreHaloRatio");
  f_TMVAPFPTrackConicalness = fPandizzleAlg.GetFloatVar("PFPTrackConicalness");
  f_TMVAPFPTrackdEdxStart = fPandizzleAlg.GetFloatVar("PFPTrackdEdxStart");
  f_TMVAPFPTrackdEdxEnd = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEnd");
  f_TMVAPFPTrackdEdxEndRatio = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEndRatio");
  f_TMVAPFPTrackPIDA = fPandizzleAlg.GetFloatVar("PFPTrackPIDA");

  fTMVAPFPMichelNHits = f_TMVAPFPMichelNHits;
  fTMVAPFPMichelElectronMVA = f_TMVAPFPMichelElectronMVA;
  fTMVAPFPMichelRecoEnergyPlane2 = f_TMVAPFPMichelRecoEnergyPlane2;
  fTMVAPFPTrackDeflecAngleSD =  f_TMVAPFPTrackDeflecAngleSD;
  fTMVAPFPTrackLength = f_TMVAPFPTrackLength;
  fTMVAPFPTrackEvalRatio = f_TMVAPFPTrackEvalRatio;
  fTMVAPFPTrackConcentration = f_TMVAPFPTrackConcentration;
  fTMVAPFPTrackCoreHaloRatio = f_TMVAPFPTrackCoreHaloRatio;
  fTMVAPFPTrackConicalness = f_TMVAPFPTrackConicalness;
  fTMVAPFPTrackdEdxStart = f_TMVAPFPTrackdEdxStart;
  fTMVAPFPTrackdEdxEnd = f_TMVAPFPTrackdEdxEnd;
  fTMVAPFPTrackdEdxEndRatio = f_TMVAPFPTrackdEdxEndRatio;

  return fPandizzleReader.EvaluateMVA("BDTG");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::PFParticle> PandSelect::GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & e){
  art::Ptr<recob::PFParticle> matched_pfp;
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(e.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToTrack Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return matched_pfp;
  }
  art::FindManyP<recob::PFParticle> fmpfpt(trackListHandle, e, fTrackModuleLabel);
  const std::vector<art::Ptr<recob::PFParticle> > sel_track_pfps = fmpfpt.at(track.key());
  if (sel_track_pfps.size() != 1){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToTrack NUMBER OF PFP MATCHED TO A TRACK DOES NOT EQUAL 1: " << sel_track_pfps.size() << std::endl;
    return matched_pfp;
  }
  matched_pfp = sel_track_pfps[0];
  return matched_pfp;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

double PandSelect::GetPandrizzleScore(art::Event& e)
{
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  if (!(e.getByLabel(fShowerModuleLabel, showerListHandle))){
    std::cout<<"Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
    return kDefDoub;
  }

  //Get the selected shower
  art::Ptr<recob::Shower> sel_shower = fRecoShowerSelector->FindSelectedShower(e);

  //If we didn't find a selected track then what's the point?
  if (!(sel_shower.isAvailable())) {
    std::cout<<"FDSelection::CCNuSelection::RunShowerSelection - no shower selected by tool"<<std::endl;
    return kDefDoub;
  }

  FillVertexInformation(e);

  FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(sel_shower, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), e));
  return pandrizzleRecord.GetMVAScore();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandSelect::FillVertexInformation(art::Event const & e){

  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(e.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }

  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);
  //Loop over the neutrinos
  if (nu_pfps.size() != 1){
    std::cout<<"FDSelection::CCNuSelection: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return; //do nothing
  }
  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, e, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());

  if (sel_pfp_vertices.size() == 0){ //Nothing to do
    return;
  }
  else if (sel_pfp_vertices.size() > 1){
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandSelect::beginJob()
{
  fTMVAPFPMichelNHits            = kDefDoub;
  fTMVAPFPMichelElectronMVA      = kDefDoub;
  fTMVAPFPMichelRecoEnergyPlane2 = kDefDoub;
  fTMVAPFPTrackDeflecAngleSD     = kDefDoub;
  fTMVAPFPTrackLength            = kDefDoub;
  fTMVAPFPTrackEvalRatio         = kDefDoub;
  fTMVAPFPTrackConcentration     = kDefDoub;
  fTMVAPFPTrackCoreHaloRatio     = kDefDoub;
  fTMVAPFPTrackConicalness       = kDefDoub;
  fTMVAPFPTrackdEdxStart         = kDefDoub;
  fTMVAPFPTrackdEdxEnd           = kDefDoub;
  fTMVAPFPTrackdEdxEndRatio      = kDefDoub;
  fTMVAPFPTrackPIDA              = kDefDoub;

  f_TMVAPFPMichelNHits            = kDefDoub;
  f_TMVAPFPMichelElectronMVA      = kDefDoub;
  f_TMVAPFPMichelRecoEnergyPlane2 = kDefDoub;
  f_TMVAPFPTrackDeflecAngleSD     = kDefDoub;
  f_TMVAPFPTrackLength            = kDefDoub;
  f_TMVAPFPTrackEvalRatio         = kDefDoub;
  f_TMVAPFPTrackConcentration     = kDefDoub;
  f_TMVAPFPTrackCoreHaloRatio     = kDefDoub;
  f_TMVAPFPTrackConicalness       = kDefDoub;
  f_TMVAPFPTrackdEdxStart         = kDefDoub;
  f_TMVAPFPTrackdEdxEnd           = kDefDoub;
  f_TMVAPFPTrackdEdxEndRatio      = kDefDoub;
  f_TMVAPFPTrackPIDA              = kDefDoub;

  fRecoNuVtxX = kDefDoub;
  fRecoNuVtxY = kDefDoub;
  fRecoNuVtxZ = kDefDoub;
}

void PandSelect::endJob()
{
}

DEFINE_ART_MODULE(PandSelect)

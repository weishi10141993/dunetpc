#ifndef PANDRIZZLEALG_H_SEEN
#define PANDRIZZLEALG_H_SEEN

///////////////////////////////////////////////
// PandrizzleAlg.h
//
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//#include "larreco/RecoAlg/ShowerEnergyAlg.h"

// c++
#include <map>
#include <memory>
#include <variant>
#include <vector>

// ROOT
#include "TTree.h"

//Custom
//#include "FDSelectionUtils.h"

//constexpr int kMaxObjects = 999;

namespace FDSelection {
    class PandrizzleAlg;
}

class FDSelection::PandrizzleAlg {
    public:

        PandrizzleAlg(const fhicl::ParameterSet& pset);
        void Run(const art::Event& evt);

        enum Vars{
            EvalRatio = 0,
            Concentration,
            CoreHaloRatio,
            Conicalness,
            dEdxBestPlane,
            Displacement,
            DCA,
            Wideness,
            EnergyDensity,
            TerminatingValue

        };
    private:

        std::map<Vars, std::unique_ptr<Float_t> >  fInputs;

};

#endif

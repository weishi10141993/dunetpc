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
#include <functional>
#include <map>
#include <memory>
#include <variant>
#include <vector>

// ROOT
#include "TMVA/Reader.h"
#include "TTree.h"

//Custom
//#include "FDSelectionUtils.h"

//constexpr int kMaxObjects = 999;

namespace FDSelection {
    class PandrizzleAlg;
}

class FDSelection::PandrizzleAlg {
    public:

        enum Vars{
            kEvalRatio = 0,
            kConcentration,
            kCoreHaloRatio,
            kConicalness,
            kdEdxBestPlane,
            kDisplacement,
            kDCA,
            kWideness,
            kEnergyDensity,
            kTerminatingValue
        };

        using InputVarsToReader = std::map<Vars, std::unique_ptr<Float_t> >;

        class Record {
            public:
                Record(const InputVarsToReader &inputVars);

                Float_t GetVar(const FDSelection::PandrizzleAlg::Vars var);

            private:
                using InputVars = std::map<Vars, Float_t>;
                InputVars fInputs;
        };

        PandrizzleAlg(const fhicl::ParameterSet& pset);

        Record Run(const art::Ptr<recob::Shower> pShower, const TVector3 &nuVertex, const art::Event& evt);

    private:
        Float_t* GetVarPtr(const FDSelection::PandrizzleAlg::Vars var);
        void SetVar(const FDSelection::PandrizzleAlg::Vars var, const Float_t value);

        std::string fPIDModuleLabel;

        TMVA::Reader fReader;
        InputVarsToReader fInputsToReader;

        Record ReturnEmptyRecord();


};

Float_t FDSelection::PandrizzleAlg::Record::GetVar(const FDSelection::PandrizzleAlg::Vars var)
{
    return (fInputs.at(var));
}

Float_t* FDSelection::PandrizzleAlg::GetVarPtr(const FDSelection::PandrizzleAlg::Vars var)
{
    return fInputsToReader.at(var).get();
}

void FDSelection::PandrizzleAlg::SetVar(const FDSelection::PandrizzleAlg::Vars var, const Float_t value)
{
    std::map<Vars, std::unique_ptr<Float_t> >::iterator itr(fInputsToReader.find(var));
    if (itr != fInputsToReader.end())
    {
        *(itr->second.get()) = value;
        std::cout<<"The value for " << var << " is " << *(itr->second.get()) << std::endl;
    }
    return;
}

#endif

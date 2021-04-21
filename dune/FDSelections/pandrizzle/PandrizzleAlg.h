#ifndef PANDRIZZLEALG_H_SEEN
#define PANDRIZZLEALG_H_SEEN

///////////////////////////////////////////////
// PandrizzleAlg.h (D. Brailsford)
//
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 

// LArSoft
#include "lardataobj/RecoBase/Shower.h"

// c++
#include <map>
#include <memory>
#include <vector>

// ROOT
#include "TMVA/Reader.h"

namespace FDSelection
{
    class PandrizzleAlg 
    {
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
                kTerminatingValue //terminates the enum and not an actual variable
            };

            using InputVarsToReader = std::map<Vars, std::unique_ptr<Float_t> >;

            class Record {
                public:
                    Record(const InputVarsToReader &inputVars, const Float_t mvaScore, const bool isFilled);

                    Float_t GetVar(const FDSelection::PandrizzleAlg::Vars var);

                    bool IsFilled();

                    Float_t GetMVAScore();

                private:
                    using InputVars = std::map<Vars, Float_t>;
                    InputVars fInputs;
                    Float_t fMVAScore;
                    bool fIsFilled;
            };

            PandrizzleAlg(const fhicl::ParameterSet& pset);

            Record RunPID(const art::Ptr<recob::Shower> pShower, const TVector3 &nuVertex, const art::Event& evt);

        private:
            Float_t* GetVarPtr(const FDSelection::PandrizzleAlg::Vars var);
            void SetVar(const FDSelection::PandrizzleAlg::Vars var, const Float_t value);

            std::string fPIDModuleLabel;

            TMVA::Reader fReader;
            InputVarsToReader fInputsToReader;

            Record ReturnEmptyRecord();


    };
}

Float_t FDSelection::PandrizzleAlg::Record::GetVar(const FDSelection::PandrizzleAlg::Vars var)
{
    return (fInputs.at(var));
}

bool FDSelection::PandrizzleAlg::Record::IsFilled()
{
    return fIsFilled;
}

Float_t FDSelection::PandrizzleAlg::Record::GetMVAScore()
{
    return fMVAScore;
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
    }
    return;
}

#endif

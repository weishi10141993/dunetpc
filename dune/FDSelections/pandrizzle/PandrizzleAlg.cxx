///////////////////////////////////////////////
// PandrizzleAlg.cxx
//
// D Brailsford
///////////////////////////////////////////////

//STL
#include <limits>

//ART
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "lardataobj/AnalysisBase/MVAPIDResult.h"

//Custom
#include "PandrizzleAlg.h"


namespace
{
    constexpr Float_t kDefValue(std::numeric_limits<Float_t>::lowest());

    using namespace FDSelection;
    void Reset(PandrizzleAlg::InputVarsToReader &inputVarsToReader)
    {
        for (PandrizzleAlg::Vars var = PandrizzleAlg::kEvalRatio; var < PandrizzleAlg::kTerminatingValue; var=static_cast<PandrizzleAlg::Vars>(static_cast<int>(var)+1)) 
        {
            auto [itr, inserted] = inputVarsToReader.try_emplace(var, std::make_unique<Float_t>(kDefValue));
            if (!inserted)
                *(itr->second.get()) = kDefValue;
        }
    }
}

FDSelection::PandrizzleAlg::Record::Record(const InputVarsToReader &inputVarsToReader, const Float_t mvaScore, const bool isFilled) :
    fMVAScore(mvaScore),
    fIsFilled(isFilled)
{
    for (PandrizzleAlg::Vars var = PandrizzleAlg::kEvalRatio; var < PandrizzleAlg::kTerminatingValue; var=static_cast<PandrizzleAlg::Vars>(static_cast<int>(var)+1)) 
    {
        fInputs.try_emplace(var, *(inputVarsToReader.at(var)));
    }
}

FDSelection::PandrizzleAlg::PandrizzleAlg(const fhicl::ParameterSet& pset) :
    fPIDModuleLabel(pset.get<std::string>("ModuleLabels.PIDModuleLabel")),
    fPandrizzleWeightFileName        (pset.get< std::string > ("PandrizzleWeightFileName")),
    fReader("",0)
{
    Reset(fInputsToReader);
    fReader.AddVariable("EvalRatio",GetVarPtr(kEvalRatio));
    fReader.AddVariable("Concentration",GetVarPtr(kConcentration));
    fReader.AddVariable("CoreHaloRatio",GetVarPtr(kCoreHaloRatio));
    fReader.AddVariable("Conicalness",GetVarPtr(kConicalness));
    fReader.AddVariable("dEdxBestPlane",GetVarPtr(kdEdxBestPlane));
    fReader.AddVariable("Displacement",GetVarPtr(kDisplacement));
    fReader.AddVariable("DCA",GetVarPtr(kDCA));
    fReader.AddVariable("Wideness",GetVarPtr(kWideness));
    fReader.AddVariable("EnergyDensity",GetVarPtr(kEnergyDensity));
    //const std::string weightFileName("Pandrizzle_TMVAClassification_BDTG.weights.xml");
    const std::string weightFileName(fPandrizzleWeightFileName);
    std::string weightFilePath;
    cet::search_path sP("FW_SEARCH_PATH");
    sP.find_file(weightFileName, weightFilePath);
    sP.find_file(weightFileName, weightFilePath);

    fReader.BookMVA("BDTG",weightFilePath);
}


FDSelection::PandrizzleAlg::Record FDSelection::PandrizzleAlg::RunPID(const art::Ptr<recob::Shower> pShower, const TVector3& nuVertex, const art::Event& evt) 
{


    const std::string weightFileName(fPandrizzleWeightFileName);
    std::string weightFilePath;
    cet::search_path sP("FW_SEARCH_PATH");
    sP.find_file(weightFileName, weightFilePath);
    sP.find_file(weightFileName, weightFilePath);
    std::cout << "pandrizzle file path: " << weightFilePath << std::endl;

    art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower> >{pShower}, evt, fPIDModuleLabel);
    art::Ptr<anab::MVAPIDResult> mvaPIDResult(findPIDResult.at(0));
    //MVAPID vars
    if (mvaPIDResult.isAvailable())
    {
        if (isnan(mvaPIDResult->evalRatio))
            SetVar(kEvalRatio, -0.5f);
        else
            SetVar(kEvalRatio, static_cast<Float_t>(mvaPIDResult->evalRatio));

        if (isnan(mvaPIDResult->concentration))
            SetVar(kConcentration, -2.f);
        else
            SetVar(kConcentration, std::min(static_cast<Float_t>(mvaPIDResult->concentration), 50.f));
        SetVar(kCoreHaloRatio, static_cast<Float_t>(mvaPIDResult->coreHaloRatio));
        SetVar(kConicalness, std::min(static_cast<Float_t>(mvaPIDResult->conicalness), 100.f));
    }
    else 
        ReturnEmptyRecord();

    //dEdx
    if (pShower->dEdx().size() > 0)
    {
        SetVar(kdEdxBestPlane, std::max(std::min(static_cast<Float_t>(pShower->dEdx().at(pShower->best_plane())), 20.f), -2.f));
    }
    else
        ReturnEmptyRecord();

    //Displacement
    SetVar(kDisplacement, std::min(static_cast<Float_t>((pShower->ShowerStart() - nuVertex).Mag()), 100.f));

    //Distance of closest approach
    double alpha((pShower->ShowerStart() - nuVertex).Dot(pShower->Direction()));
    TVector3 r(pShower->ShowerStart() + alpha*pShower->Direction());
    SetVar(kDCA, std::min(static_cast<Float_t>((r-nuVertex).Mag()), 50.f));

    //Wideness
    Float_t wideness(static_cast<Float_t>(pShower->OpenAngle()/pShower->Length()));
    if (isnan(wideness))
        SetVar(kWideness, -0.01f);
    else
        SetVar(kWideness, std::min( wideness, 0.1f));

    //Energy density
    if (pShower->Energy().size() > 0)
    {
        Float_t volume(static_cast<Float_t>((M_PI * pShower->Length() * pShower->Length() * pShower->Length() * std::tan(pShower->OpenAngle()))/3.));
        Float_t energyDensity(std::min(std::max(static_cast<Float_t>(pShower->Energy().at(2))/volume, -0.1f), 5.f));
        if (isnan(energyDensity))
            energyDensity = -0.1f;
        SetVar(kEnergyDensity, energyDensity);
    }
    else
        ReturnEmptyRecord();

    return Record(fInputsToReader, fReader.EvaluateMVA("BDTG"), true);
}

FDSelection::PandrizzleAlg::Record FDSelection::PandrizzleAlg::ReturnEmptyRecord()
{
    Reset(fInputsToReader);
    return Record(fInputsToReader, kDefValue, false);
}

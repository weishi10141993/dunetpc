///////////////////////////////////////////////
// PandrizzleAlg.cxx
//
// Reco and true PID stuff up
// D Brailsford & M Wallbank, June 2017
///////////////////////////////////////////////

#include <limits>

#include "PandrizzleAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


constexpr Float_t kDefValue(std::numeric_limits<Float_t>::lowest());

FDSelection::PandrizzleAlg::PandrizzleAlg(const fhicl::ParameterSet& pset) :
    fReader("",0)
//fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg"))
{
    for ( Vars var = kEvalRatio; var < kTerminatingValue; var=static_cast<Vars>(static_cast<int>(var)+1)) 
    {
        fInputs.try_emplace(var, std::make_unique<Float_t>(kDefValue));
    }
    fReader.AddVariable("EvalRatio",GetVarPtr(kEvalRatio));
    fReader.AddVariable("Concentration",GetVarPtr(kConcentration));
    fReader.AddVariable("CoreHaloRatio",GetVarPtr(kCoreHaloRatio));
    fReader.AddVariable("Conicalness",GetVarPtr(kConicalness));
    fReader.AddVariable("dEdxBestPlane",GetVarPtr(kdEdxBestPlane));
    fReader.AddVariable("Displacement",GetVarPtr(kDisplacement));
    fReader.AddVariable("DCA",GetVarPtr(kDCA));
    fReader.AddVariable("Wideness",GetVarPtr(kWideness));
    fReader.AddVariable("EnergyDensity",GetVarPtr(kEnergyDensity));
    const std::string weightFileName("Pandrizzle_TMVAClassification_BDTG.weights.xml");
    std::string weightFilePath;
    cet::search_path sP("FW_SEARCH_PATH");
    sP.find_file(weightFileName, weightFilePath);
    std::cout<<"The found path is " << weightFilePath << std::endl;
    sP.find_file(weightFileName, weightFilePath);

    fReader.BookMVA("BDTG",weightFilePath);
}


void FDSelection::PandrizzleAlg::Run(const art::Event& evt) 
{
    std::cout<<"starting"<<std::endl;
    std::cout<<"value is: " << fReader.EvaluateMVA("BDTG")<<std::endl;;
    std::cout<<"done"<<std::endl;
}

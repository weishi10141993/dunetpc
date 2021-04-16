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

#include "TMVA/Reader.h"

constexpr Float_t kDefValue(std::numeric_limits<Float_t>::lowest());

FDSelection::PandrizzleAlg::PandrizzleAlg(const fhicl::ParameterSet& pset) //:
  //fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg"))
{
    for ( const auto e : { Vars::EvalRatio, Vars::TerminatingValue } ) 
    {
        fInputs.try_emplace(e, std::make_unique<Float_t>(kDefValue));
    }

}


void FDSelection::PandrizzleAlg::Run(const art::Event& evt) 
{
    /*
    std::unique_ptr<float> var1 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var2 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var3 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var4 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var5 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var6 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var7 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var8 = std::unique_ptr<float>(new float(0.2f));
    std::unique_ptr<float> var9 = std::unique_ptr<float>(new float(0.2f));

    fInputs[Var1] = std::move(var1);
    fInputs[Var2] = std::move(var2);
    fInputs[Var3] = std::move(var3);
    fInputs[Var4] = std::move(var4);
    fInputs[Var5] = std::move(var5);
    fInputs[Var6] = std::move(var6);
    fInputs[Var7] = std::move(var7);
    fInputs[Var8] = std::move(var8);
    fInputs[Var9] = std::move(var9);
    */


 //   std::cout<<"starting"<<std::endl;
 // TMVA::Reader fReader("", 1);
 // fReader.AddVariable("EvalRatio",&*(std::get<std::unique_ptr<float> >(fInputs[Var1]).get()));
 // fReader.AddVariable("Concentration",&*(std::get<std::unique_ptr<float> >(fInputs[Var2]).get()));
 // fReader.AddVariable("CoreHaloRatio",&*(std::get<std::unique_ptr<float> >(fInputs[Var3]).get()));
 // fReader.AddVariable("Conicalness",&*(std::get<std::unique_ptr<float> >(fInputs[Var4]).get()));
 // fReader.AddVariable("dEdxBestPlane",&*(std::get<std::unique_ptr<float> >(fInputs[Var5]).get()));
 // fReader.AddVariable("Displacement",&*(std::get<std::unique_ptr<float> >(fInputs[Var6]).get()));
 // fReader.AddVariable("DCA",&*(std::get<std::unique_ptr<float> >(fInputs[Var7]).get()));
 // fReader.AddVariable("Wideness",&*(std::get<std::unique_ptr<float> >(fInputs[Var8]).get()));
 // fReader.AddVariable("EnergyDensity",&*(std::get<std::unique_ptr<float> >(fInputs[Var9]).get()));
 // fReader.BookMVA("BDTG","/dune/app/users/dbrailsf/oscillation/nu_mu/cutsel/soft/srcs/dunetpc/dune/FDSelections/weights/Pandrizzle_TMVAClassification_BDTG.weights.xml");
 // //std::cout<<"value is: " << fReader.EvaluateMVA("BDTG")<<std::endl;;
 // //std::cout<<"done"<<std::endl;

}

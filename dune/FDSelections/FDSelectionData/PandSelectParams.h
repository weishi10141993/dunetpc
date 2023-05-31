#ifndef PAND_SELECT_PARAMS_H
#define PAND_SELECT_PARAMS_H

#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

namespace pandselect
{
  class PandSelectParams
  {
  public:

      double selTrackPandizzleScore;
      double selShowerPandrizzleScore;
      double selShowerJamPandrizzleScore;

      dune::EnergyRecoOutput energyRecoNumu;
      dune::EnergyRecoOutput energyRecoNue;
  };
}

#endif

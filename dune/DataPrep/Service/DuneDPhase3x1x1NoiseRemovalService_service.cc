// DuneDPhase3x1x1NoiseRemovalService_service.cc

#include "DuneDPhase3x1x1NoiseRemovalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

//**********************************************************************

DuneDPhase3x1x1NoiseRemovalService::
DuneDPhase3x1x1NoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) :
    fRoiThreshold( pset.get<float>("RoiThreshold") ),
    fRoiPadLow( pset.get<int>("RoiPadLow") ),
    fRoiPadHigh( pset.get<int>("RoiPadHigh") ),
    fGeometry( &*art::ServiceHandle<geo::Geometry>() )
{
}
//**********************************************************************

int DuneDPhase3x1x1NoiseRemovalService::update(AdcChannelDataMap& datamap) const {
  const std::string myname = "DuneDPhase3x1x1NoiseRemovalService:update: ";
  if ( datamap.size() == 0 ) {
    std::cout << myname << "WARNING: No channels found." << std::endl;
    return 1;
  }

  unsigned int nsig = 0;
  bool first = true;
  for (const AdcChannelDataMap::value_type& ent : datamap)
  {
    const AdcChannelData& data = ent.second;
    if (first) { nsig = data.samples.size(); first = false; }
    else if (data.samples.size() != nsig)
    {
      std::cout << myname << "WARNING: Channels have inconsistent sample counts." << std::endl;
      return 2;
    }
  }
  
  if (nsig == 0)
  {
    std::cout << myname << "WARNING: No ADC samples found." << std::endl;
    return 3;
  }

  std::cout << myname << "Processing noise removal..." << std::endl;

  auto ch_groups = makeDaqGroups(32);
  removeCoherent(ch_groups, datamap);

  std::cout << myname << "...done." << std::endl;

  return 0;
}
//**********************************************************************

void DuneDPhase3x1x1NoiseRemovalService::removeCoherent(const GroupChannelMap & ch_groups, AdcChannelDataMap& datamap) const
{
  if (datamap.empty()) return;

  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector<size_t> ch_averaged(n_samples);
  std::vector<double> correction(n_samples);

  for (const auto & entry : ch_groups)
  {
    const auto & channels = entry.second;
    std::fill(ch_averaged.begin(), ch_averaged.end(), 0);
    std::fill(correction.begin(), correction.end(), 0);
    for (unsigned int ch : channels)
    {
        //std::cout << "  ch:" << ch;
        auto iacd = datamap.find(ch);
        if (iacd == datamap.end()) continue;

        const AdcChannelData & adc = iacd->second;
        auto mask = roiMask(adc, fRoiThreshold);

        for (size_t s = 0; s < n_samples; ++s)
        {
            if (!mask[s]) { continue; }

            AdcFlag flag = adc.flags.size() ? adc.flags[s] : AdcGood;
            if (flag != AdcGood) { continue; }

            correction[s] += adc.samples[s];
            ch_averaged[s]++;
        }
    }
    //std::cout << std::endl;
    for (size_t s = 0; s < n_samples; ++s)
    {
        if (ch_averaged[s] > 0) { correction[s] /= ch_averaged[s]; }
        //std::cout << " " << correction[s];
    }
    //std::cout << std::endl << std::endl;
    for (unsigned int ch : channels)
    {
        auto iacd = datamap.find(ch);
        if (iacd == datamap.end()) continue;

        AdcChannelData & adc = iacd->second;
        for (size_t s = 0; s < n_samples; ++s)
        {
            if (ch_averaged[s] > 0)
            {
                adc.samples[s] -= correction[s];
            }
        }
    }
  }
}
//**********************************************************************

std::vector<bool> DuneDPhase3x1x1NoiseRemovalService::roiMask(const AdcChannelData & adc, AdcSignal thr) const
{
  std::vector<bool> mask(adc.samples.size(), true);

  bool inroi = false;
  for (int i = 0; i < (int)adc.samples.size(); ++i)
  {
    AdcSignal sig = adc.samples[i];
    if (inroi)
    {
      if (sig > thr) { mask[i] = false; }
      else
      {
        for (int p = 1; p <= fRoiPadHigh; ++p) { if (i + p < (int)mask.size()) { mask[i + p] = false; } }
        inroi = false;
      }
    }
    else
    {
      if (sig > thr )
      {
        for (int p = 0; p < fRoiPadLow; ++p) { if (i - p >= 0) { mask[i + p] = false; } }
        mask[i] = false;
        inroi = true;
      }
    }
  }
  return mask;
}

GroupChannelMap DuneDPhase3x1x1NoiseRemovalService::makeDaqGroups(size_t gsize) const
{
  GroupChannelMap groups;

  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

  const unsigned int nchan = fGeometry->Nchannels();
  for (unsigned int ch = 0; ch < nchan; ++ch)
  {
    if (chStatus.IsPresent(ch) && !chStatus.IsNoisy(ch)) { groups[get311Chan(ch) / gsize].push_back(ch); }
    else { std::cout << "skip channel " << ch << std::endl; } // can remove once verified it is working ok
  }

  return groups;
}
//**********************************************************************

GroupChannelMap DuneDPhase3x1x1NoiseRemovalService::makeGroups(size_t gsize) const
{
  GroupChannelMap groups;

  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

  const unsigned int nchan = fGeometry->Nchannels();
  for (unsigned int ch = 0; ch < nchan; ++ch)
  {
    if (chStatus.IsPresent(ch) && !chStatus.IsNoisy(ch)) { groups[ch / gsize].push_back(ch); }
    else { std::cout << "skip channel " << ch << std::endl; } // can remove once verified it is working ok
  }

  return groups;
}
//**********************************************************************

size_t DuneDPhase3x1x1NoiseRemovalService::get311Chan(size_t LAr_chan)
{
  size_t crate = LAr_chan / 320;
  size_t Chan311;

  LAr_chan = 8*(LAr_chan/8+1)-LAr_chan%8 -1;

  if(crate == 0)
  {
    LAr_chan = 32*(LAr_chan/32+1)-LAr_chan%32 -1;
    size_t card = 4 - ((LAr_chan / 32) % 5);
    if(LAr_chan > 159)
    {
        size_t shift = 31 - (LAr_chan % 32);
        Chan311 = (2*card)*32 + shift;
    }
    else
    {
       size_t shift = 31 - (LAr_chan % 32);
       Chan311 = (2*card + 1)*32 + shift;
    }
  }
  else
  {
     size_t new_LAr_chan = LAr_chan - crate*320;
     size_t card = ((new_LAr_chan / 32) % 5);
     if(new_LAr_chan > 159)
     {
        size_t shift = new_LAr_chan % 32;
        Chan311 = (2*card)*32 + shift;
     }
     else
     {
       size_t shift = new_LAr_chan % 32;
       Chan311 = (2*card + 1)*32 + shift;
     }
     Chan311 = Chan311 + crate*320;
  } // end of if/else statementi

  return Chan311;
}
//**********************************************************************

std::ostream& DuneDPhase3x1x1NoiseRemovalService::print(std::ostream& out, std::string prefix) const
{
  out << prefix << "DuneDPhase3x1x1NoiseRemovalService:  ...info" << std::endl;
  return out;
}
//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneDPhase3x1x1NoiseRemovalService, AdcNoiseRemovalService)

//**********************************************************************

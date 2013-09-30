////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter35.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/GeoObjectSorter35.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

namespace geo{

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in APA configuration
  //   same as standard
  static bool sortCryo35(const CryostatGeo* c1, const CryostatGeo* c2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.}; 
    c1->LocalToWorld(local, xyz1);
    c2->LocalToWorld(local, xyz2);

    return xyz1[0] < xyz2[0];   
  }


  //----------------------------------------------------------------------------
  // Define sort order for tpcs in APA configuration.
  static bool sortTPC35(const TPCGeo* t1, const TPCGeo* t2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    // The goal is to number TPCs first in the x direction so that,
    // in the case of APA configuration, TPCs 2c and 2c+1 make up APA c.
    // then numbering will go in y then in z direction.

    // First sort all TPCs into same-z groups
    if(xyz1[2] < xyz2[2]) return true;
 
    // Within a same-z group, sort TPCs into same-y groups
    if(xyz1[2] == xyz2[2] && xyz1[1] < xyz2[1]) return true;
 
    // Within a same-z, same-y group, sort TPCs according to x
    if(xyz1[2] == xyz2[2] && xyz1[1] == xyz2[1] && xyz1[0] < xyz2[0]) return true;      

    // none of those are true, so return false
    return false;
  }


  //----------------------------------------------------------------------------
  // Define sort order for planes in APA configuration
  //   same as standard, but implemented differently
  static bool sortPlane35(const PlaneGeo* p1, const PlaneGeo* p2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    p1->LocalToWorld(local, xyz1);
    p2->LocalToWorld(local, xyz2);

    return xyz1[0] > xyz2[0];
  }

  //----------------------------------------------------------------------------
  // we want the wires to be sorted such that the smallest corner wire
  // on the readout end of a plane is wire zero, with wire number
  // increasing away from that wire.
  
  // Since 35t has an APA which is both above and below the world origin,
  // we cannot use the APA trick. If we could ask where wire 0 was, we could
  // still do this in a single implimentation, but we aren't sure what wire
  // center we will be getting, so this reversed sorting must be handled
  // at the plane level where there is one vertical center.
  // If the plane center is above, count from top down (the top stacked and
  // largest APAs) If the plane is below (bottom stacked APA) count bottom up
  struct sortWire35{

    std::string detVersion;

    sortWire35(std::string detv)
      : detVersion(detv)
    {}

    bool operator()(WireGeo* w1, WireGeo* w2){
      double xyz1[3] = {0.};
      double xyz2[3] = {0.};
      
      w1->GetCenter(xyz1); w2->GetCenter(xyz2);
      
      // immedieately take care of vertical wires regardless of which TPC
      // vertical wires should always have same y, and always increase in z direction
      if( xyz1[1]==xyz2[1] && xyz1[2]<xyz2[2] ) return true;
      
      ///////////////////////////////////////////////////////////
      // Hard code a number to tell sorting when to look
      // for top/bottom APAs and when to look for only one
      bool InVertSplitRegion;
      if(detVersion=="lbne35t")     InVertSplitRegion = xyz1[2] > 76.35;   // the old
      if(detVersion=="lbne35t4apa") InVertSplitRegion = ((51 < xyz1[2])    // the new
							 && (xyz1[2] < 102));
      ///////////////////////////////////////////////////////////

      // we want the wires to be sorted such that the smallest corner wire
      // on the readout end of a plane is wire zero, with wire number
      // increasing away from that wire. 
      
      if( InVertSplitRegion ){
        
        // if a bottom TPC, count from bottom up
        if( xyz1[1]<0 && xyz1[1] < xyz2[1] ) return true; 
        
        // if a top TPC, count from top down
        if( xyz1[1]>0 && xyz1[1] > xyz2[1] ) return true;
        
      } 
      else {
        
        if( xyz1[1] > xyz2[1] ) return true;
        
      }

      return false;
    }
  };

  //----------------------------------------------------------------------------
  GeoObjectSorter35::GeoObjectSorter35(fhicl::ParameterSet const& p)
    : fDetVersion(p.get< std::string >("DetectorVersion", "lbne35t4apa"))
  {
  }

  //----------------------------------------------------------------------------
  GeoObjectSorter35::~GeoObjectSorter35()
  {
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorter35::SortCryostats(std::vector<geo::CryostatGeo*> & cgeo) const
  {
    std::sort(cgeo.begin(), cgeo.end(), sortCryo35);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorter35::SortTPCs(std::vector<geo::TPCGeo*>  & tgeo) const
  {
    
    std::sort(tgeo.begin(), tgeo.end(), sortTPC35);

    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorter35::SortPlanes(std::vector<geo::PlaneGeo*> & pgeo,
				     geo::DriftDirection_t  const& driftDir) const
  {
    // sort the planes to increase in drift direction
    // The drift direction has to be set before this method is called.  It is set when
    // the CryostatGeo objects are sorted by the CryostatGeo::SortSubVolumes method
    if     (driftDir == geo::kPosX) std::sort(pgeo.rbegin(), pgeo.rend(), sortPlane35);
    else if(driftDir == geo::kNegX) std::sort(pgeo.begin(),  pgeo.end(),  sortPlane35);
    else if(driftDir == geo::kUnknownDrift)
      throw cet::exception("TPCGeo") << "Drift direction is unknown, can't sort the planes";

    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorter35::SortWires(std::vector<geo::WireGeo*> & wgeo) const
  {
    sortWire35 sw35(fDetVersion);

    std::sort(wgeo.begin(), wgeo.end(), sw35);

    return;
  }

}

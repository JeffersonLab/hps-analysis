//
// Created by Maurik Holtrop on 12/31/22.
//

#ifndef HPS_ANALYSIS_PARTICLES_H
#define HPS_ANALYSIS_PARTICLES_H

#include "TObject.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
using namespace ROOT;
using namespace ROOT::VecOps;

class Particles {

public:
   RVec<bool> fiducial_cut(RVec<int> ix, RVec<int> iy){
      RVec<bool> out;
      for(size_t i=0;i< ix.size();++i){
         if(
               !(ix[i] <= -23 || ix[i] >= 23) && /* Cut out the left and right side */
               !(iy[i] <= -6 || iy[i] >= 6)   && /* Cut out the top and bottom row */
               !(iy[i] >= -1 && iy[i] <= 1)   && /* Cut out the first row around the gap */
               !(iy[i] >= -2 && iy[i] <= 2 && ix[i] >= -11 && ix[i] <= 1)
               )
         {
            out.push_back(true);
         }else{
            out.push_back(false);
         }
      }
      return out;
   }

ClassDef(Particles,1)
};


#endif //HPS_ANALYSIS_PARTICLES_H

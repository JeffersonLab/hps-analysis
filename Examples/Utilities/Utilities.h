//
// Created by Maurik Holtrop on 12/31/22.
//

#ifndef HPS_ANALYSIS_UTILITIES_H
#define HPS_ANALYSIS_UTILITIES_H

#include "TObject.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
using namespace ROOT;
using namespace ROOT::VecOps;

class Utilities {

public:
   static RVec<bool> fiducial_cut(RVec<int> ix, RVec<int> iy){
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

   static RVec<bool> fiducial_cut_extended(RVec<int> ix, RVec<int> iy, std::vector< std::pair<int,int> > exclude){
      RVec<bool> out;
      for(size_t i=0;i< ix.size();++i){
         if(
               !(ix[i] <= -23 || ix[i] >= 23) && /* Cut out the left and right side */
               !(iy[i] <= -6 || iy[i] >= 6)   && /* Cut out the top and bottom row */
               !(iy[i] >= -1 && iy[i] <= 1)   && /* Cut out the first row around the gap */
               !(iy[i] >= -2 && iy[i] <= 2 && ix[i] >= -11 && ix[i] <= 1)
               )
         {
            bool do_we_include = true;
            for(const auto& cc : exclude) {
               if (ix[i] >= cc.first - 1 && ix[i] <= cc.first + 1 && iy[i] >= cc.second - 1 && iy[i] <= cc.second + 1) {
                  // ix,iy is within the 3x3 exclusion zone for this crystal.
                  do_we_include = false;
                  break;
               }
               // This is where the decision to skip ix==0 bites. If first==1 then first-2 should be excluded too,
               // same for first==-1, then first+2 should be excluded.
               if (cc.first == 1 && ix[i] == -1 && iy[i] >= cc.second - 1 && iy[i] <= cc.second + 1) {
                  do_we_include = false;
                  break;
               }
               if (cc.first == -1 && ix[i] == 1 && iy[i] >= cc.second - 1 && iy[i] <= cc.second + 1) {
                  do_we_include = false;
                  break;
               }
            }
            out.push_back(do_we_include);
         }else{
            out.push_back(false);
         }
      }
      return out;
   }


ClassDef(Utilities, 1)
};


#endif //HPS_ANALYSIS_UTILITIES_H

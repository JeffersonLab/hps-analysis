//
//  FEESelector.hpp
//  BaseAna
//
//

#ifndef FEESelector_hpp
#define FEESelector_hpp

#include <iostream>
#include <TObject.h>

#include "BaseAna.h"

class FEESelector: public TObject {
  
public:
  int    hitcut;                 ///< Cut on the minimum number of required hits in a cluster. Set to 0 for no cut. [#]
  double min_seed_hit_energy;    ///< Minimum energy for the seed hit. [GeV]
  double min_cluster_energy;     ///< Minimum energy for the cluster. [GeV]
  double min_seed_cluster_ratio; ///< Minimum ratio for SeedEnergy/ClusterEnergy [#.#]
  double min_time;               ///< Minimum Time for the cluster [ns].
  double max_time;               ///< Maxumum Time for the cluster [ns].
  
private:
  BaseAna *base;                 ///< Store a pointer to the BaseAna derived class.
  
public:
  
  FEESelector(){ cerr << "ERROR -- Class FEESelector should not be created without a pointer to a BaseAna class\n"; };
  FEESelector(BaseAna *ba);
  void Set_Engineering2015(void);
  
  bool Select(int np);
  bool Select(HpsParticle *part);
  
  ~FEESelector(){};
  
  
  
  ClassDef(FEESelector,0);
};

#endif /* FEESelector_hpp */

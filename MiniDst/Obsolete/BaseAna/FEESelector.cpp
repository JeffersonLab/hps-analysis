/*!
 * \class  FEESelector
 * \brief  The base analysis class for hps-dst based analysis projects.
 *
 * This is a simple selector class. Initialize it with a pointer to a BaseAna class or "this" for the
 * Process() method in a BaseAna derived class. Then the Select() method will return true if the event
 * passed the selection criteria.
 */

#include "FEESelector.h"

ClassImp(FEESelector);

///
/// Standard constuctor.
/// This constructor should *always* be called with a BasaAna derived pointer.
///
/// Default parameters will be set with Set_Engineering2015().
///
FEESelector::FEESelector(BaseAna *ba): base(ba){
  Set_Engineering2015();
};


///
/// Set the defaults for the Engineering run 2015 for clean (calibration) FEE.
/// These are the default values used for this class.
///
void FEESelector::Set_Engineering2015(void){
  hitcut = 0;                  // default is not to require a min hit
  min_seed_hit_energy = 0.4;   // Default for clean FEE for 1.05 GeV data
  min_cluster_energy  = 0.6;   // Default for clean FEE for
  min_seed_cluster_ratio = 0.4;
  min_time  = 30.;
  max_time  = 70.;
};


///
/// Select will return true if the event seems to be an FEE
/// The class must be inialized with a BaseAna pointer. Currently
/// it does not check if this pointer is OK and has an event!
///
/// @param[in] np The number of the particle in the event.
/// \return{  True if particle is considered an FEE, false otherwise }
///
bool FEESelector::Select(int np){
  HpsParticle *part=base->GetFSParticle(np);
  return(Select(part));
}

///
/// Select will return true if the event seems to be an FEE
/// The class must be inialized with a BaseAna pointer. Currently
/// it does not check if this pointer is OK and has an event!
///
/// @param[in] part A particle in the event to test if it is an FEE.
/// \return{  True if particle is considered an FEE, false otherwise }
///
/// TODO: Arrange test in most likely to fail to least likely to fail for speed optimization.
///
bool FEESelector::Select(HpsParticle *part){

  HpsEvent *event =base->event;
  
  if( !event->isSingle0Trigger() && !event->isSingle1Trigger() ) return false; // Throw out all but Single0 or Single1 Triggers.
  
  TRefArray *clus_ref = part->getClusters();

  if(clus_ref->GetEntries()!= 1) return false;  // Must have one and only one cluster.
  EcalCluster *clus =(EcalCluster *)clus_ref->At(0);

  TRefArray *hits = clus->getEcalHits();
  int nhits = hits->GetEntries();
  
  if(nhits <= hitcut) return false;             // Must have at least hitcut hits.
  
  double cluster_energy = clus->getEnergy();
  
  if(cluster_energy < min_cluster_energy) return false;  // Cluster must have at least min_cluster_energy
  
  EcalHit *seed = clus->getSeed();
  double seed_energy = seed->getEnergy();
  
  if( seed_energy < min_seed_hit_energy) return false;                    // Must have at least min_seed_hit_energy
  if( seed_energy < cluster_energy*min_seed_cluster_ratio) return false;  // seed_energy/cluster_energy > min_seed_cluter_ratio.
                                                                          // Note that muliplication is faster and safer then division.
  
  double time = clus->getClusterTime();
  
  if( time < min_time || time > max_time ) return false;
  
  return true;  // All tests passed, so I think this was an FEE.
  
}

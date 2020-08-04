/**
 * @file: EcalCluster.cxx
 * @author: Omar Moreno <omoreno1@ucsc.edu>
 * @section Institution \n
 *          Santa Cruz Institute for Particle Physics
 *          University of California, Santa Cruz
 * @date: February 19, 2013
 */

#include "EcalCluster.h"

ClassImp(EcalCluster)

EcalCluster::EcalCluster()
    : TObject(),
      ecal_hits(new TRefArray()),
      n_ecal_hits(0),
      x(0),
      y(0),
      z(0),
      energy(0),
      cluster_time(0) {
}

EcalCluster::EcalCluster(const EcalCluster &ecal_cluster_obj)
    : TObject(),
      ecal_hits(new TRefArray()), 
      seed_hit(ecal_cluster_obj.seed_hit), 
      n_ecal_hits(ecal_cluster_obj.n_ecal_hits),
      x(ecal_cluster_obj.x), 
      y(ecal_cluster_obj.y), 
      z(ecal_cluster_obj.z),
      energy(ecal_cluster_obj.energy),
      cluster_time(ecal_cluster_obj.cluster_time) {

    *ecal_hits = *ecal_cluster_obj.ecal_hits; 
}

EcalCluster::~EcalCluster() {
    Clear();
    delete ecal_hits;   
}

EcalCluster &EcalCluster::operator=(const EcalCluster &ecal_cluster_obj) {
    
    // Check for self-assignment
    if(this == &ecal_cluster_obj) return *this;

    TObject::operator=(ecal_cluster_obj);
    Clear();

    this->n_ecal_hits = ecal_cluster_obj.n_ecal_hits;
    this->x = ecal_cluster_obj.x;
    this->y = ecal_cluster_obj.y;
    this->z = ecal_cluster_obj.z;
    this->energy = ecal_cluster_obj.energy;
    this->cluster_time = ecal_cluster_obj.cluster_time;
    
    // These will be commented out until they are part of the variables
    // in the reconstruction.
    //this->m2 = ecal_cluster_obj.m2; 
    //this->m3 = ecal_cluster_obj.m3; 

    ecal_hits = new TRefArray(); 
    *ecal_hits = *ecal_cluster_obj.ecal_hits;

    seed_hit = ecal_cluster_obj.seed_hit;     

    return *this;
}

void EcalCluster::Clear(Option_t* /*option*/) {
    TObject::Clear();
    ecal_hits->Delete();
    seed_hit = NULL;
    n_ecal_hits = 0;
}

void EcalCluster::setPosition(const float* position) {
    x = position[0];
    y = position[1];
    z = position[2];
}

void EcalCluster::addHit(EcalHit* hit) { 
    ++n_ecal_hits;

    if(seed_hit.GetObject() == NULL
            || (static_cast<EcalHit*>(seed_hit.GetObject())->getEnergy() < hit->getEnergy())) {
        seed_hit = hit;
        cluster_time = hit->getTime();
    }

    ecal_hits->Add(hit);
}

std::vector<double> EcalCluster::getPosition() const {
    std::vector<double> position(3, 0); 
    position[0] = x; 
    position[1] = y; 
    position[2] = z; 
    
    return position;    
}

EcalHit* EcalCluster::getSeed() const {
    return static_cast<EcalHit*>(seed_hit.GetObject());
}

TRefArray* EcalCluster::getEcalHits() const {
    return ecal_hits;
}


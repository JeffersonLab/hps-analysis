/**
 * @file MCParticle.cxx
 * @brief Class which implements an MC particle that stores information about
 *        tracks from the simulation.
 * @author: Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "MCParticle.h"

ClassImp(MCParticle)

MCParticle::MCParticle()
    : TObject() {
}

MCParticle::~MCParticle() {
    TObject::Clear();
    delete daughters_;
    delete parents_;
}

void MCParticle::Clear(Option_t *) {

    TObject::Clear();

    daughters_->Delete();
    parents_->Delete();

    energy_ = 0;
    pdg_id_ = 0;
    gen_status_ = -1;
    time_ = 0;
    x_ = 0;
    y_ = 0;
    z_ = 0;
    end_x_ = 0;
    end_y_ = 0;
    end_z_ = 0;
    px_ = 0;
    py_ = 0;
    pz_ = 0;
    mass_ = 0;
    charge_ = 0;
}

void MCParticle::Print(Option_t *) const {
    std::cout << "[ MCParticle ]: " 
              << "\t Energy: " << energy_ << "\n" 
              << "\t PDG ID: " << pdg_id_ << "\n" 
              << "\t Generator Status: " << gen_status_ << "\n" 
              << "\t Time: " << time_ << "\n" 
              << "\t Vertex: [ " << x_ << ", " << y_ << ", " << z_ << " ]\n " 
              << "\t EndPoint: [ " << end_x_ << ", " << end_y_ << ", " << end_z_ << " ]\n " 
              << "\t Momentum: [ " << px_ << ", " << py_ << ", " << pz_ << " ]\n "
              << "\t Mass: " << mass_ << "\n" 
              << "\t Daughter count: " << daughters_->GetEntries() << "\n"
              << "\t Parent count: " << parents_->GetEntries() << "\n"
              << std::endl;
}

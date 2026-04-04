/**
 * @file: EcalHit.cxx
 * @author: Omar Moreno <omoreno1@ucsc.edu>
 * @section Institution \n
 *          Santa Cruz Institute for Particle Physics
 *          University of California, Santa Cruz
 * @date: January 28, 2014
 */

#include "EcalHit.h"

ClassImp(EcalHit)

EcalHit::EcalHit() 
    : TObject(),
      index_x(0),
      index_y(0),
      energy(0), 
      hit_time(0) {
}

EcalHit::~EcalHit() {
    Clear(); 
}

void EcalHit::Clear(Option_t* /* options */) {
    TObject::Clear(); 
}

void EcalHit::setCrystalIndices(int x_index, int y_index) {
    index_x = x_index;
    index_y = y_index;
}

// These will be commented out until a need for them is demonstrated.
/*void EcalHit::setPosition(const double* position) {
    x = position[0];
    y = position[1];
    z = position[2];
}*/

/*std::vector<double> EcalHit::getPosition() const {
    std::vector<double> position(3,0);
    position[0] = x;
    position[1] = y;
    position[2] = z;
    return position;
}*/

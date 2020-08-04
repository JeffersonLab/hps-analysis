/**
 * @file LegacyGblTrack.cxx
 * @brief Class used to encapsulate GBL track refit information.  This class
 *        was used when the GBL refit was done at the DST level and is now 
 *        deprecated.
 * @author Per Hansson Adrian <phansson@slac.stanford.edu>
 *         SLAC
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date February 3, 2014
 *
 */

#include "LegacyGblTrack.h"

ClassImp(LegacyGblTrack)

LegacyGblTrack::LegacyGblTrack() 
    : TObject(),
      seed_track(NULL), 
      cov_matrix(5, 5),
      d0(0),
      phi0(0),
      omega(0),
      tan_lambda(0),
      z0(0),      
      chi_squared(-1.), 
      ndof(0), 
      px(0), 
      py(0), 
      pz(0) {
}

LegacyGblTrack::~LegacyGblTrack() {
    Clear(); 
}

void LegacyGblTrack::Clear(Option_t *option) { 
    TObject::Clear(); 
}

void LegacyGblTrack::setTrackParameters(const double d0, 
        const double phi0,
        const double omega,
        const double tan_lambda,
        const double z0) {

    this->d0 = d0; 
    this->phi0 = phi0; 
    this->omega = omega; 
    this->tan_lambda = tan_lambda;
    this->z0 = z0;  
}

void LegacyGblTrack::setMomentumVector(double x, double y, double z) {
    px = x;
    py = y;
    pz = z;
}

std::vector<double> LegacyGblTrack::getMomentum() {
    std::vector<double> p(3, 0);
    p[0] = px; 
    p[1] = py; 
    p[2] = pz;
    return p; 
}

void LegacyGblTrack::toString() {
    std::cout << "LegacyGblTrack: " << std::endl;
    std::cout << "    (omega, tan_lambda, phi0, d0, z0): " 
        << getOmega() 
        << "," << getTanLambda() 
        << "," << getPhi0() 
        << "," << getD0() 
        << "," << getZ0() << std::endl;
}

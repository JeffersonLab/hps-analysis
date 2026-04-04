/**
 *
 * @file SvtTrack.h
 * @brief Class used to describe an HPS SVT track.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date February 19, 2013
 *
 */

#include "SvtTrack.h"
#include "SvtHit.h"

ClassImp(SvtTrack)

SvtTrack::SvtTrack(){}

SvtTrack::SvtTrack(const SvtTrack &svtTrackObj)
: n_hits(svtTrackObj.n_hits),
track_volume(svtTrackObj.track_volume),
type(svtTrackObj.type),
d0(svtTrackObj.d0),
phi0(svtTrackObj.phi0),
omega(svtTrackObj.omega),
tan_lambda(svtTrackObj.tan_lambda),
z0(svtTrackObj.z0),
chi_squared(svtTrackObj.chi_squared),
track_time(svtTrackObj.track_time),
x_at_ecal(svtTrackObj.x_at_ecal),
y_at_ecal(svtTrackObj.x_at_ecal),
z_at_ecal(svtTrackObj.x_at_ecal) {
    
    *svt_hits = *svtTrackObj.svt_hits;
    fs_particle = svtTrackObj.fs_particle;
    gbl_track = svtTrackObj.gbl_track;
    memcpy(&isolation, svtTrackObj.isolation, sizeof(isolation));
    memcpy(&covmatrix, svtTrackObj.covmatrix, sizeof(covmatrix));
}


SvtTrack &SvtTrack::operator=(const SvtTrack &svtTrackObj) {
    
    // Check for self-assignment
    if(this == &svtTrackObj) return *this;
    
    TObject::operator=(svtTrackObj);
    Clear();
    delete svt_hits;
    
    this->n_hits = svtTrackObj.n_hits;
    this->track_volume = svtTrackObj.track_volume;
    this->type = svtTrackObj.type;
    this->d0 = svtTrackObj.d0;
    this->phi0 = svtTrackObj.phi0;
    this->omega = svtTrackObj.omega;this->tan_lambda = svtTrackObj.tan_lambda;this->z0 = svtTrackObj.z0;
    this->chi_squared = svtTrackObj.chi_squared;
    this->track_time = svtTrackObj.track_time;
    this->x_at_ecal = svtTrackObj.x_at_ecal;
    this->y_at_ecal = svtTrackObj.y_at_ecal;
    this->z_at_ecal = svtTrackObj.z_at_ecal;
    
    svt_hits = new TRefArray();
    *svt_hits = *svtTrackObj.svt_hits;
    fs_particle = svtTrackObj.fs_particle;
    gbl_track = svtTrackObj.gbl_track;
    memcpy(&isolation, svtTrackObj.isolation, sizeof(isolation));
    memcpy(&covmatrix, svtTrackObj.covmatrix, sizeof(covmatrix));
    return *this;
}

SvtTrack::~SvtTrack() {
    Clear();
    delete svt_hits;
}

void SvtTrack::Clear(Option_t* /* option */) {
    TObject::Clear();
    svt_hits->Delete();
    memset(isolation, 0, sizeof(isolation));
    memset(covmatrix, 0, sizeof(covmatrix));
    n_hits = 0;
}

void SvtTrack::setTrackParameters(double D0, double Phi0, double Omega,
                                  double Tan_lambda, double Z0) {
    d0         = D0;
    phi0       = Phi0;
    omega      = Omega;
    tan_lambda = Tan_lambda;
    z0         = Z0;
}

void SvtTrack::setPositionAtEcal(const double* position) { 
    x_at_ecal = position[0];
    y_at_ecal = position[1];
    z_at_ecal = position[2];
}

int SvtTrack::getCharge() const{
    // From ReconParticleDriver line 464:
    // Derive the charge of the particle from the track.
    int charge = ((omega) > 0) ? 1 : (((omega) < 0) ? -1 : 0); // Compute the charge for By>0
    return( -charge); // In HPS the By field is down, so we need to flip the charge.
}

std::vector<double> SvtTrack::getMomentum(double bfield) const{
    // If the bfield is given, compute the momentum from the curvature
    // and the 3-vector from the angles.
    // For the 2016 run the bfield was 0.5234000000000001
    // See org.lcsim.base.BaseTrack
    // Compute momentum from parameters and magnetic field.
    if( abs(bfield)>0.00000001){
        if( abs(omega)<0.0000001) return {0,0,0};  // Return zero instead of ~infinite.
        double Pt = abs(1./omega)*bfield* 2.99792458E-4; //
        double b0 = Pt*cos(phi0);
        double b1 = Pt*sin(phi0);
        double b2 = Pt*tan_lambda;
        return {b1,b2,b0};  // Coordinate transformation.
    }else{
        if (fs_particle == NULL) return {0, 0, 0};
        HpsParticle *part = static_cast<HpsParticle*>(fs_particle.GetObject());
        if( part == nullptr) return {0,0,0};
        return part->getMomentum();
    }
}

std::vector<double> SvtTrack::getPositionAtEcal() const{ 
    std::vector<double> position = {
        x_at_ecal,
        y_at_ecal,
        z_at_ecal
    };
    return position;
}

void SvtTrack::addHit(SvtHit* hit) {
    ++n_hits;
    svt_hits->Add(static_cast<TObject*>(hit));
}

TRefArray* SvtTrack::getSvtHits() const {
    return svt_hits;
}

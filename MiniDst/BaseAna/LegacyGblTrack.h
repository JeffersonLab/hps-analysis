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

#ifndef __LEGACY_GBL_TRACK_H__
#define __LEGACY_GBL_TRACK_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <iostream>

//------------//
//--- ROOT ---//
//------------//
#include <TObject.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TRef.h>
#include <TMatrixD.h>

// Forward declarations
class SvtTrack; 

class LegacyGblTrack : public TObject {

    // TODO: Add more documentation
    // TODO: Add copy constructor
    // TODO: Add copy assignment

    public:
        
        /**
         * Default constructor
         */
        LegacyGblTrack();

        /**
         * Destructor
         */
        ~LegacyGblTrack();

        /**
         *
         */
        void Clear(Option_t *option=""); 

        /**
         * Set the GBL track parameters
         *
         * @param d0: 
         * @param phi0:
         * @param omega:
         * @param tan_lambda:
         * @param z0:
         */
        void setTrackParameters(const double d0, 
                const double phi0,
                const double omega,
                const double tan_lambda,
                const double z0);

        /**
         * Set the chi^2 of the fit to the track.
         *
         * @param chi_squared : The chi^2 of the fit to the track.
         */        
        void setChi2(const double Chi_squared) { chi_squared = Chi_squared; }

        /**
         *
         */
        void setNdf(const double Ndof) { ndof = Ndof; }

        /**
         *
         */
        void setMomentumVector(const double px, const double py, const double pz);

        /**
         *
         */
        void setCov(const TMatrixD& Cov_matrix) { cov_matrix = Cov_matrix; }

        /**
         *
         */
        void setSeedTrack(SvtTrack* Seed_track) { seed_track = static_cast<TObject*>(Seed_track); }

        /**
         *
         */
        std::vector<double> getMomentum(); 

        /**
         *
         */
        double getD0() const { return d0; }

        /**
         *
         */
        double getPhi0() const { return phi0; }

        /**
         *
         */
        double getOmega() const { return omega; }

        /**
         *
         */
        double getTanLambda() const { return tan_lambda; }

        /**
         *
         */
        double getZ0() const { return z0; }

        /**
         *
         */
        double getChi2() const { return chi_squared; }

        /**
         *
         */
        double getNdf() const { return ndof; }

        /**
         *
         */
        TRef getSeedTrack() const { return seed_track; }

        /**
         *
         */
        void toString();

        ClassDef(LegacyGblTrack, 1) //Track class for use with GBL

    public:

            TRef seed_track; 
            TMatrixD cov_matrix;

            double d0;
            double phi0;
            double omega;
            double tan_lambda;
            double z0;
            double chi_squared;
            double ndof;
            double px;
            double py;
            double pz;

}; // LegacyGblTrack

#endif // _LEGACY_GBL_TRACK_H_

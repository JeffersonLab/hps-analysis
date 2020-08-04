/**
 * 
 * @file GblTrack.h
 * @brief Class used to describe an HPS SVT GBL track.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date December 11, 2015
 * 
 */

#ifndef __GBL_TRACK_H__
#define __GBL_TRACK_H__

//-----------------//
//--- HPS Event ---//
//-----------------//
#include "SvtTrack.h"

class GblTrack : public SvtTrack { 

    public: 
        
        /** Constructor */
        GblTrack();

        /**
         * Set the lambda kink of the given layer.
         *
         * @param layer Layer number associated with the given lambda kink.
         * @param lambda_kink The lambda kink value.
         */
        void setLambdaKink(const int layer, const double lambda_kink) { lambda_kinks[layer] = lambda_kink; }

        /**
         * Set the phi kink of the given layer.
         *
         * @param layer Layer number associated with the given phi kink.
         * @param phi_kink The phi kink value.
         */
        void setPhiKink(const int layer, const double phi_kink) { phi_kinks[layer] = phi_kink; }

        /**
         * Set the seed track that was refit to create this GBL track.
         *
         * @param Seed_track The seed track (SvtTrack object)
         */
        void setSeedTrack(SvtTrack* Seed_track) { seed_track = static_cast<TObject*>(Seed_track); }

        /**
         * Get the lambda kink value of the given layer.
         *
         * @param layer The SVT layer of interest.
         * @return The lambda kink value of the given layer.
         */
        double getLambdaKink(const int layer) const { return lambda_kinks[layer]; }

        /**
         * Get the phi kink value of the given layer.
         *
         * @param layer The SVT layer of interest.
         * @return The phi kink value of the given layer.
         */
        double getPhiKink(const int layer) const { return phi_kinks[layer]; }

        /**
         * Get the seed track that was refit to create this GBL track.
         *
         * @return The seed track.
         */ 
        SvtTrack* getSeedTrack() const { return static_cast<SvtTrack*>(seed_track.GetObject()); }

        ClassDef(GblTrack, 2)

    public:

        /** Reference to corresponding seed track */
        TRef seed_track; 

        /** Array used to store the lambda kinks for each of the sensor layers. */
        double lambda_kinks[14];
        
        /** Array used to store the phi kinks for each of the sensor layers. */
        double phi_kinks[14];
};

#endif // __GBL_TRACK_H__

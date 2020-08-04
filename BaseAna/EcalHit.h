/**
 * @file: EcalHit.h
 * @author: Omar Moreno <omoreno1@ucsc.edu>
 * @section Institution \n
 *          Santa Cruz Institute for Particle Physics
 *          University of California, Santa Cruz
 * @date: January 28, 2014
 */

#ifndef __ECAL_HIT_H__
#define __ECAL_HIT_H__

//------------//
//--- ROOT ---//
//------------//
#include <TObject.h>

//-----------------//
//--- HPS Event ---//
//-----------------//
#include "CalorimeterHit.h"

class EcalHit : public TObject, public CalorimeterHit {

    public:

        /**
         * Default Constructor
         */
        EcalHit(); 

        /**
         * Destructor
         */
        ~EcalHit();

        /**
         *
         */
        void Clear(Option_t* option=""); 

        /**
         * Set the energy of the hit in GeV.
         *
         * @param Energy : The energy of the hit
         */
        void setEnergy(const double Energy){ energy = Energy; };
       
        /**
         * Set the time of the hit
         *
         * @param Hit_time : The time of the hit
         */
        void setTime(const double Hit_time) { hit_time = Hit_time; };

        /**
         * Set the indices of the crystal
         *
         * @param index_x : The index along x
         * @param index_y : The index along y
         */
        void setCrystalIndices(int index_x, int index_y);

        /**
         * Get the energy of the hit in GeV.
         *
         * @return The energy of the hit
         */
        double getEnergy() const { return energy; };
       
        /**
         * Get the time of the hit.
         *
         * @return The time of the hit
         */
        double getTime() const { return hit_time; };

        /**
         * Get the crystal index along x.
         *
         * @return The index along x
         */
        int getXCrystalIndex() const { return index_x; };
        
        /**
         * Get the crystal index along y.
         *
         * @return The index along y
         */
        int getYCrystalIndex() const { return index_y; };

        ClassDef(EcalHit, 1);

        // Comment these out until a need for it is demonstrated
        //void setPosition(const double*);
        //std::vector<double> getPosition() const;

    public:

        int index_x;
        int index_y;

        //int x;
        //int y;
        //int z;

        double energy;
        double hit_time;

}; // EcalHit

#endif // _ECAL_HIT_H_

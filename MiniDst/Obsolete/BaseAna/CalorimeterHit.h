/**
 * @file: CalorimeterHit.h
 * @author: Omar Moreno <omoreno1@ucsc.edu>
 * @section Institution \n
 *			Santa Cruz Institute for Particle Physics
 *			University of California, Santa Cruz
 * @date: January 28, 2014
 *
 */

#ifndef __CALORIMETER_HIT_H__
#define __CALORIMETER_HIT_H__

class CalorimeterHit { 

	public: 

		virtual ~CalorimeterHit() {}; 

    /**
     * Get the energy of the hit in GeV.
     *
     * @return The energy of the hit
     */
		virtual double getEnergy() const = 0;
		
    /**
     * Set the energy of the hit in GeV.
     *
     * @param energy : The energy of the hit
     */
		virtual void setEnergy(const double energy) = 0;

        // Comments these out for now until a need for them is demonstrated
		//virtual void setPosition(const double*) = 0;
		//virtual std::vector<double> getPosition() const = 0;

}; // CalorimeterHit

#endif // __CALORIMETER_HIT_H__

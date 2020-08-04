/**
 *	@author:	Omar Moreno <omoreno1@ucsc.edu>
 *	@section institution
 *				Santa Cruz Institute for Particle Physics
 *				University of California, Santa Cruz
 *	@version:	v 1.0
 *	@date:		February 11, 2014
 *
 */

#ifndef _HPS_MC_PARTICLE_H_
#define _HPS_MC_PARTICLE_H_

//-- ROOT ---//
//-----------//
#include <TObject.h>

class HpsMCParticle : public TObject {

	public:

		HpsMCParticle();
		HpsMCParticle(const HpsMCParticle &);
		~HpsMCParticle();
		HpsMCParticle &operator=(const HpsMCParticle &);

		void Clear(Option_t *option="");

		void setPDG(const int Pdg){ pdg = Pdg; };
		void setCharge(const int Charge){ charge = Charge; };
		void setGeneratorStatus(const int Generator_status){ generator_status = Generator_status; };
		void setEnergy(const double Energy){ energy = Energy; };
		void setMass(const double Mass){ mass = Mass; };
		void setMomentum(const double*);
		void setEndpoint(const double*);

		int getPDG() const { return pdg; };
		int getCharge() const { return charge; };
		int getGeneratorStatus() const { return generator_status; };
		double getEnergy() const { return energy; };
		double getMass() const { return mass; };
		std::vector<double> getMomentum() const;
		std::vector<double> getEndpoint() const;

		ClassDef(HpsMCParticle, 2)

	public:

		int pdg;
		int charge;
		int generator_status;

		double energy;
		double mass;
		double px;
		double py;
		double pz;
		double endpt_x;
		double endpt_y;
		double endpt_z;

}; // HpsMCParticle

#endif // _HPS_MC_PARTICLE_H_

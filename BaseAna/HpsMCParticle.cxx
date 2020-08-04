/**
 *	@author:	Omar Moreno <omoreno1@ucsc.edu>
 *	@section institution
 *				Santa Cruz Institute for Particle Physics
 *				University of California, Santa Cruz
 *	@version:	v 1.0
 *	@date:		February 11, 2014
 *
 */

//--- HPS Event ---//
//-----------------//
#include "HpsMCParticle.h"

ClassImp(HpsMCParticle)

HpsMCParticle::HpsMCParticle()
	: TObject(), pdg(0), charge(0), generator_status(0),
	  energy(0), mass(0), px(0), py(0), pz(0), endpt_x(0), endpt_y(0), endpt_z(0)
{}

HpsMCParticle::HpsMCParticle(const HpsMCParticle &mcParticleObj)
	:TObject(), pdg(mcParticleObj.pdg), charge(mcParticleObj.charge),
	 generator_status(mcParticleObj.generator_status),
	 energy(mcParticleObj.energy), mass(mcParticleObj.mass),
	 px(mcParticleObj.px), py(mcParticleObj.py), pz(mcParticleObj.pz),
	endpt_x(mcParticleObj.endpt_x), endpt_y(mcParticleObj.endpt_y), endpt_z(mcParticleObj.endpt_z)
{
}

HpsMCParticle::~HpsMCParticle()
{
	Clear();
}

HpsMCParticle &HpsMCParticle::operator=(const HpsMCParticle &mcParticleObj)
{
	// Check for self-assignment
	if(this == &mcParticleObj) return *this;

	TObject::operator=(mcParticleObj);
	Clear();

	this->pdg = mcParticleObj.pdg;
	this->charge = mcParticleObj.charge;
	this->generator_status = mcParticleObj.generator_status;
	this->energy = mcParticleObj.energy;
	this->mass = mcParticleObj.mass;
	this->px = mcParticleObj.px;
	this->py = mcParticleObj.py;
	this->pz = mcParticleObj.pz;
	this->endpt_x = mcParticleObj.endpt_x;
	this->endpt_y = mcParticleObj.endpt_y;
	this->endpt_z = mcParticleObj.endpt_z;

	return *this;
}

void HpsMCParticle::Clear(Option_t* /*option*/)
{
	TObject::Clear();
}

void HpsMCParticle::setMomentum(const double* momentum)
{
	px = momentum[0];
	py = momentum[1];
	pz = momentum[2];
}

void HpsMCParticle::setEndpoint(const double* end_point)
{
	endpt_x = end_point[0];
	endpt_y = end_point[1];
	endpt_z = end_point[2];
}

std::vector<double> HpsMCParticle::getMomentum() const
{
	std::vector<double> momentum(3, 0);
	momentum[0] = px;
	momentum[1] = py;
	momentum[2] = pz;
	return momentum;
 }

std::vector<double> HpsMCParticle::getEndpoint() const
{
	std::vector<double> end_point(3, 0);
	end_point[0] = endpt_x;
	end_point[1] = endpt_y;
	end_point[2] = endpt_z;
	return end_point;
}

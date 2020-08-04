/**
 *	@author:	Omar Moreno <omoreno1@ucsc.edu>
 *	@section institution
 *				Santa Cruz Institute for Particle Physics
 *				University of California, Santa Cruz
 *	@version:	v 1.0
 *	@date:		January 28, 2014
 *
 */


#ifndef _CLUSTER_H_
#define _CLUSTER_H_

class Cluster {

	public:
		
		virtual ~Cluster(){}; 
		
		virtual std::vector<double> getPosition() const = 0; 
		virtual double getEnergy() const = 0; 

		virtual void setPosition(const float*) = 0;
		virtual void setEnergy(const double) = 0;

};

#endif // _CLUSTER_H_

/**
 *
 * @author: 	Omar Moreno <omoreno1@ucsc.edu>
 * @section institution
 * 				Santa Cruz Institute for Particle Physics
 * 				University of California, Santa Cruz
 * @version:    v 0.1
 * @date:       February 19, 2013
 */

#include "SvtHit.h"

ClassImp(SvtHit)

SvtHit::SvtHit()
	: TObject(), layer(0), x(0), y(0), z(0),
	  cxx(0), cxy(0), cxz(0), cyy(0), cyz(0), czz(0), time(0)
{}

SvtHit::~SvtHit()
{
    Clear(); 
}

void SvtHit::Clear(Option_t* /* options */)
{
    TObject::Clear();  
}

void SvtHit::setPosition(const double* position)
{
    x = position[0]; 
    y = position[1];
    z = position[2]; 
}

void SvtHit::setCovarianceMatrix(const std::vector<float> covariance_matrix)
{
	cxx = covariance_matrix[0];
	cxy = covariance_matrix[1];
	cxz = covariance_matrix[2];
	cyy = covariance_matrix[3];
	cyz = covariance_matrix[4];
	czz = covariance_matrix[5];
}

std::vector<double> SvtHit::getPosition() const
{
	std::vector<double> position(3, 0);
	position[0] = x;
	position[1] = y;
	position[2] = z;

	return position;
}

std::vector<double> SvtHit::getCovarianceMatrix() const
{
	std::vector<double> covariance_matrix(6, 0);
	covariance_matrix[0] = cxx;
	covariance_matrix[1] = cxy;
	covariance_matrix[2] = cxz;
	covariance_matrix[3] = cyy;
	covariance_matrix[4] = cyz;
	covariance_matrix[5] = czz;

	return covariance_matrix;
}

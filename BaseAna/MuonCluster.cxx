/**
 *
 * @author: 	Omar Moreno <omoreno1@ucsc.edu>
 * @section institution
 * 				Santa Cruz Institute for Particle Physics
 * 				University of California, Santa Cruz
 * @version:    v 0.1
 * @date:       April 09, 2013
 */

#include "MuonCluster.h"

ClassImp(MuonCluster)

MuonCluster::MuonCluster()
{}

MuonCluster::MuonCluster(const MuonCluster &muonClusterObj)
{}

MuonCluster::~MuonCluster()
{
	Clear();
}

MuonCluster &MuonCluster::operator=(const MuonCluster &muonClusterObj)
{
	// Check for self-assignment
	if(this == &muonClusterObj) return *this;

	TObject::operator=(muonClusterObj);
	Clear();

	return *this;
}

void MuonCluster::Clear(Option_t* /* option */)
{
	TObject::Clear();
}

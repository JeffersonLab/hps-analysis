/**
 *
 * @author: 	Omar Moreno <omoreno1@ucsc.edu>
 * @section institution
 * 				Santa Cruz Institute for Particle Physics
 * 				University of California, Santa Cruz
 * @version:    v 0.1
 * @date:       April 09, 2013
 */

#ifndef _MUON_CLUSTER_H_
#define _MUON_CLUSTER_H_

//--- C++ ---//
#include <iostream>

//--- ROOT ---//
#include <TObject.h>
#include <TClonesArray.h>

class MuonCluster : public TObject {

	public:

		MuonCluster();
		MuonCluster(const MuonCluster &muonClusterObj);
		virtual ~MuonCluster();
		MuonCluster &operator=(const MuonCluster &muonClusterObj);

		void Clear(Option_t *option = "");

		ClassDef(MuonCluster, 1); 

}; // MuonCluster

#endif // _MUON_CLUSTER_H_

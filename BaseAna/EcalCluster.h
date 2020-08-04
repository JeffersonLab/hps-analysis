/**
 * @file: EcalCluster.h
 * @author: Omar Moreno <omoreno1@ucsc.edu>
 * @section Institution \n
 *          Santa Cruz Institute for Particle Physics
 *          University of California, Santa Cruz
 * @date: February 19, 2013
 */

#ifndef __ECAL_CLUSTER_H__
#define __ECAL_CLUSTER_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <iostream>

//-----------//
//-- ROOT ---//
//-----------//
#include <TObject.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TRef.h>

//-----------------//
//--- HPS Event ---//
//-----------------//
#include "Cluster.h"
#include "EcalHit.h"

class EcalCluster : public TObject, public Cluster {

    public:

        /**
         * Default Destructor
         */
        EcalCluster();
        
        /**
         * Copy Constructor
         *
         * @param ecal_cluster_obj : EcalCluster that will be cloned 
         */
        EcalCluster(const EcalCluster &ecal_cluster_obj);

        /**
         * Destructor
         */
        ~EcalCluster();

        /**
         *
         */
        EcalCluster &operator=(const EcalCluster &);

        /**
         *
         */
        void Clear(Option_t *option="");

        /**
         * Add a reference to an Ecal hit composing with this cluster.
         *
         * @param hit : Ecal hit composing with this cluster
         */
		void addHit(EcalHit* hit); 

        /**
         * Set the position of the Ecal cluster.
         *
         * @param position : The position of the Ecal cluster
         */
        void setPosition(const float* position);
        
        /**
         * Set the energy of the Ecal cluster.
         *
         * @param Energy : The energy of the Ecal cluster.
         */
        void setEnergy(const double Energy) { energy = Energy; };
        
        /**
         * Set the time of the cluster.  This currently comes from the hit time
         * of the seed hit.
         *
         * @param cluster_time : The time of the cluster
         */
        void setClusterTime(const double cluster_time) { this->cluster_time = cluster_time; };
        
        // These will be commented out until they are part of the variables
        // in the reconstruction.
        //void setM2(const double m2){ this->m2 = m2; };
        //void setM3(const double m3){ this->m3 = m3; };

        /**
         * Get the position of the Ecal cluster.
         *
         * @return The position of the Ecal cluster
         */
		std::vector<double> getPosition() const;  
        
        /**
         * Get the energy of the Ecal cluster.
         *
         * @return The energy of the Ecal cluster.
         */
        double getEnergy() const { return energy; };
       
        /**
         * Get the time of the cluster.  This currently comes from the hit time
         * of the seed hit.
         *
         * @return The time of the cluster
         */
        double getClusterTime() const { return cluster_time; };
        

        /**
         * Get the seed hit of the cluster.
         *
         * @return The seed hit of the cluster
         */
        EcalHit* getSeed() const;
        
        /**
         * Get the Ecal hits composing this cluster.
         *
         * @return An array of references to the Ecal hits composing this cluster
         */
        TRefArray* getEcalHits() const;
		
        //double getM2() const { return m2; };
        //double getM3() const { return m3; };

        ClassDef(EcalCluster, 1);	

    public:
		
		TRefArray* ecal_hits; 
		TRef seed_hit; 

        int n_ecal_hits; 

        double x; 
        double y; 
        double z;
        double energy; 
        double cluster_time;
        
        //double m2;
        //double m3;

}; // EcalCluster

#endif // _ECAL_CLUSTER_H_

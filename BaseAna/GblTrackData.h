/**
 * @section purpose: Stores track information needed by GBL
 * @author:     Per Hansson Adrian <phansson@slac.stanford.edu>
 * @date:       February 3, 2014
 *
 */

#ifndef _GBL_TRACK_DATA_H_
#define _GBL_TRACK_DATA_H_

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <iostream>
#include <sstream>
#include <math.h>

//------------//
//--- ROOT ---//
//------------//
#include <TObject.h>
#include <TClonesArray.h>
#include <TMatrixD.h>
#include <TRefArray.h>
#include <TRef.h>
#include <TMath.h>

//-----------------//
//--- HPS Event ---//
//-----------------//
#include "GblStripData.h"
#include "SvtTrack.h"

// TODO: Add documentation
class GblTrackData : public TObject {

    public:

        /**
         * Default construtor
         */    
        GblTrackData();

        /**
         * Destructor
         */
        ~GblTrackData();
       
        /**
         *
         */ 
        void Clear(Option_t *option="");
        
        /**
         *
         */
        void addStrip(GblStripData* strip);

        /**
         *
         */
        void setTrack(SvtTrack* Seed_track){ seed_track = Seed_track; }; 
        
        /**
         *
         */
        void setPrjPerToCl(const unsigned int, const unsigned int, const double);

        /**
         *
         */
        int getNStrips() const { return n_gbl_strip_hits; }
        
        /**
         *
         */
        GblStripData* getStrip(const int& i) const {
            return static_cast<GblStripData*>(m_gbl_strip_hits->At(i));
        }
        
        /**
         *
         */
        double getKappa() const;
        
        /**
         *
         */
        double getTheta() const;
        
        /**
         *
         */
        double getPhi()   const;
        
        /**
         *
         */
        double getD0()    const;
        
        /**
         *
         */
        double getZ0()    const;
        
        /**
         *
         */
        TMatrixD getPrjPerToCl() const { return m_prjPerToCl; }
        
        /**
         *
         */
        std::string toString() const;

        ClassDef(GblTrackData,1) //Track information needed by GBL
    
    public:

        TRef seed_track; 
        
        TRefArray* m_gbl_strip_hits;
        TMatrixD m_prjPerToCl;
        
        int n_gbl_strip_hits;
        
        double m_kappa;
        double m_theta;
        double m_phi;
        double m_d0;
        double m_z0;


}; // GblTrackData

#endif // _GBL_TRACK_DATA_H_

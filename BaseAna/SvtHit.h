/**
 *
 * @author: 	Omar Moreno <omoreno1@ucsc.edu>
 * @section institution
 * 				Santa Cruz Institute for Particle Physics
 * 				University of California, Santa Cruz
 * @version:    v 0.1
 * @date:       February 19, 2013
 */

#ifndef _SVT_HIT_H_
#define _SVT_HIT_H_

//--- C++ ---//
#include <iostream>

//-- ROOT ---//
#include <TObject.h>
#include <TClonesArray.h>

class SvtHit : public TObject { 

    public: 

		SvtHit();
		virtual ~SvtHit();
        void Clear(Option_t *option="");
        
        void setLayer(const int Layer){ layer = Layer; };
        void setPosition(const double*);
        void setCovarianceMatrix(std::vector<float>);
        void setTime(const double Time) { time = Time; };

        double getLayer() const { return layer; };
        std::vector<double> getPosition() const;
        std::vector<double> getCovarianceMatrix() const;
        double getTime() const { return time; };

        ClassDef(SvtHit, 1);	
    
    public:

        int layer;

        double x; 
        double y; 
        double z;
        double cxx;
        double cxy;
        double cxz;
        double cyy;
        double cyz;
        double czz;
        double time;

};

#endif

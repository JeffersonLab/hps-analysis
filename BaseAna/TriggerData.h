/**
 * @file: TriggerData.h
 * @author: Omar Moreno <omoreno1@ucsc.edu>
 * @section Institution \n
 *          Santa Cruz Institute for Particle Physics
 *          University of California, Santa Cruz
 * @date:   May 21, 2015
 */

#ifndef __TRIGGER_DATA_H__
#define __TRIGGER_DATA_H__

//------------//
//--- LCIO ---//
//------------//
#include "EVENT/LCGenericObject.h"

class TriggerData { 
    
    public: 

        /**
         * Constructor
         *
         * @param trigger_data : The LCGenericObeject that is being used to 
         *                       store the data from the TI
         */
        TriggerData(EVENT::LCGenericObject* trigger_data);

        /**
         *
         */
        double getTime() const { return time_stamp; };

        /**
         *
         */
        bool isSingle0Trigger() const { return single0; };

        /**
         *
         */
        bool isSingle1Trigger() const { return single1; };

        /**
         *
         */
        bool isPair0Trigger() const { return pair0; };
        
        /**
         *
         */
        bool isPair1Trigger() const { return pair1; };

        /**
         *
         */
        bool isPulserTrigger() const { return pulser; };

    public:

        void parseTriggerData(EVENT::LCGenericObject* trigger_data);  

        long time_stamp; 
        bool single0;
        bool single1;
        bool pair0;
        bool pair1; 
        bool pulser; 
};

#endif // __TRIGGER_DATA_H__

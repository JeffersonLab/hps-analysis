/**
 *
 * @file HpsEvent.h
 * @brief Event class used to encapsulate event information and physics 
 *        collections.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date February 19, 2013
 *
 */

#ifndef __HPS_EVENT_H__
#define __HPS_EVENT_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <stdexcept>

//----------//
//   ROOT   //
//----------//
#include <TObject.h>
#include <TClonesArray.h>

//---------------//
//   HPS Event   //
//---------------//
#include "EcalCluster.h"
#include "Cluster.h"
#include "SvtTrack.h"
#include "GblTrack.h"
#include "SvtHit.h"
#include "EcalHit.h"
#include "HpsParticle.h"
#include "HpsMCParticle.h"
#include "MCParticle.h"


class HpsEvent : public TObject { 

    // TODO: Add documentation

    public:

        /** Constructor */
        HpsEvent();

        /** 
         * Copy Constructor
         *
         * @param hpsEventObj An HpsEvent object 
         */
        HpsEvent(const HpsEvent &hpsEventObj);

        /** Destructor */
        ~HpsEvent();    
        
        /**
         * Copy assignment operator 
         *
         * @param hpsEventObj An HpsEvent object
         */ 
        HpsEvent &operator=(const HpsEvent &hpsEventObj);

        /** 
         */
        void Clear(Option_t *option="");
        
        /** 
         * Add a track, {@link SvtTrack} object, to the event.
         *
         * @return A pointer to the {@link SvtTrack} object that was added.
         */
        SvtTrack*       addTrack();
       
        /** 
         * Add a GBL track, {@link GblTrack} object, to the event.
         *
         * @return A pointer to the {@link GblTrack} object that was added.
         */
        GblTrack*       addGblTrack();

        /** 
         * Add an SVT hit, {@link SvtHit} object, to the event.
         *
         * @return A pointer to the {@link SvtHit} object that was added.
         */
        SvtHit*         addSvtHit();
        
        /** 
         * Add an Ecal cluster, {@link EcalCluster} object, to the event.
         *
         * @return A pointer to the {@link EcalCluster} object that was added.
         */
        EcalCluster*    addEcalCluster();
        
        /** 
         * Add a Ecal hit, {@link EcalHit} object, to the event.
         *
         * @return A pointer to the {@link EcalHit} object that was added.
         */
        EcalHit*        addEcalHit();
        
        /**
         * Create an {@link HpsParticle} object with the given 
         * {@link HpsParticle::ParticleType} and add it to the
         * event.
         *
         * @param type The type of particle that is being requested e.g. 
         *             HpsParticle::FINAL_STATE_PARTICLE.
         * @return A pointer to the newly created HpsParticle object. 
         */
        HpsParticle*    addParticle(HpsParticle::ParticleType type); 
        
        /** */
        MCParticle*  addMCParticle();
        
        //--- Setters ---//
        //---------------//

        /** Set the event number. */
        void setEventNumber(int Event_number){ event_number = Event_number; };

        /**
         * Set the event time stamp. The event time is currently the Unix time 
         * stamp associated with the event.
         *
         * @param Event_time The Unix time stamp of the event.
         */
        void setEventTime(const long Event_time) { event_time = Event_time; };

        /**
         *
         */
        void setPair0Trigger(const int Pair0_trigger) { pair0_trigger = Pair0_trigger; };
        
        /**
         *
         */
        void setPair1Trigger(const int Pair1_trigger) { pair1_trigger = Pair1_trigger; };
        
        /**
         *
         */
        void setPulserTrigger(const int Pulser_trigger) { pulser_trigger = Pulser_trigger; };

        /**
         * Set the event RF time.
         *
         * @param channel The channel from which the RF time was retrieved.
         * @param rf_time The event RF time. 
         */
        void setRfTime(const int channel, const double rf_time) { rf_times[channel] = rf_time; }; 

        /**
         *
         */
        void setRunNumber(int Run_number){ run_number = Run_number; };

        /**
         *
         */
        void setSingle0Trigger(const int Single0_trigger) { single0_trigger = Single0_trigger; };

        /**
         *
         */
        void setSingle1Trigger(const int Single1_trigger) { single1_trigger = Single1_trigger; };
       
        /**
         * Set the state of the SVT bias during the event i.e. was it on or 
         * off? 
         *
         * @param Svt_bias_state The state of the SVT bias. It's set to 0 if
         *                       the bias was off or 1 if it was on.
         */ 
        void setSvtBiasState(const int Svt_bias_state) { svt_bias_state = Svt_bias_state; };

        /**
         * Set the flag indicating whether the event was affected by SVT burst
         * noise.
         *
         * @param Svt_burstmode_noise Flag indicating whether an event was affected
         *                        by SVT burst noise.  It's set to 0 if it was
         *                        or 1 if it wasn't.
         */
        void setSvtBurstModeNoise(const int Svt_burstmode_noise) { svt_burstmode_noise = Svt_burstmode_noise; };

        /**
         * Set the flag indicating whether the SVT headers had errors.
         *
         * @param Svt_event_header_state Flag indicating whether the SVT event
         *                               headers had errors.
         *
         */
        void setSvtEventHeaderState(const int Svt_event_header_state) { svt_event_header_state
                                                                            = Svt_event_header_state; };
        
        /**
         * Set the flag indicating whether the SVT latency was correct
         * during an event.
         *
         * @param Svt_latency_state Flag indicating whether the SVT latency
         *                          was correct during an event.
         */
        void setSvtLatencyState(const int Svt_latency_state) { svt_latency_state = Svt_latency_state; };


        /**
         * Set the state of indicating whether the SVT was open or closed 
         * during an event. 
         *
         * @param Svt_position_state The state indicating whether the SVT was
         *                           open or closed. It's set to 0 if the SVT 
         *                           was open or 1 if it was closed. 
         */ 
        void setSvtPositionState(const int Svt_position_state) { svt_position_state = Svt_position_state; };
  
        //--- Getters ---//
        //---------------//

        /** */
        EcalCluster*   getEcalCluster(int) const;
        
        /** */
        EcalHit*       getEcalHit(int);
        
        /** */
        MCParticle* getMCParticle(int);
        
        /**
         *
         */ 
        int getEventNumber() const  { return event_number; };

        /**
         * Get the event time stamp i.e. the Unix time stamp associated
         * with the event.
         *
         * @return The event time stamp.  
         */
        long getEventTime() const { return event_time; };

        /** */
        int getNumberOfEcalHits() const { return n_ecal_hits; };
  
        /** */
        int getNumberOfEcalClusters()   const  { return n_ecal_clusters; };
       
  
        /**
         * Get the number of MC particles ({@link MCParticle}) objects in the event.
         */
        int getNumberOfMCParticles() const { return n_mc_particles; };
  
        /**
         * Get the number of particles ({@link HpsParticle} objects) of the
         * given {@link HpsParticle::ParticleType} in the event.
         *
         * @param type The type of particle that is being requested e.g. 
         *             HpsParticle::FINAL_STATE_PARTICLE.
         * @return The number of particles of the given type in the event.
         */
        int getNumberOfParticles(HpsParticle::ParticleType type) const; 
        
        /** */
        int getNumberOfSvtHits()        const  { return n_svt_hits; };

        /** */
        int getNumberOfTracks()         const  { return n_tracks; };
    
        /** */
        int getNumberOfGblTracks() const { return n_gbl_tracks; }; 

        /**
         * Get the particle object ({@link HpsParticle}) of the given type and
         * at the given position in the container from the event.  See the class
         * {@link HpsParticle} for the type of particles that are available.
         *
         * @param type The type of particle that is being requested e.g. 
         *             HpsParticle::FINAL_STATE_PARTICLE.
         * @param particle_index The position of the particle in the container
         * @return An HpsParticle object of the given type and at the given 
         *         position in the container
         * @throws std::runtime_error if the given type is invalid
         *
         */
        HpsParticle*   getParticle(HpsParticle::ParticleType type, int particle_index); 
      
        /**
         * Get the RF time.
         *
         * @param channel The channel associated with the RF time.
         * @return The RF time. 
         */ 
        double getRfTime(const int channel) const { return rf_times[channel]; }; 

        /**
         *
         */ 
        int getRunNumber() const  { return run_number; };

        /** */
        SvtHit*        getSvtHit(int);

        /** */
        SvtTrack*      getTrack(int);
       
        /** */
        GblTrack*      getGblTrack(int); 

        //---//
        
        /**
         * Indicate whether a pair0 trigger was registered for the event.
         *
         * @return Returns true if a pair0 trigger was registered for the,
         *         false otherwise.
         */
        bool isPair0Trigger() const { return pair0_trigger == 1; };

        /**
         * Indicate whether a pair1 trigger was registered for the event.
         *
         * @return Returns true if a pair1 trigger was registered for the,
         *         false otherwise.
         */      
        bool isPair1Trigger() const { return pair1_trigger == 1; };

        /**
         * Indicate whether a pulser (random) trigger was registered for the 
         * event.
         *
         * @return Returns true if a pulser trigger was registered for the,
         *         false otherwise.
         */
        bool isPulserTrigger() const { return pulser_trigger == 1; };

        /**
         * Indicate whether a single0 trigger was registered for the event.
         *
         * @return Returns true if a single0 trigger was registered for the,
         *         false otherwise.
         */
        bool isSingle0Trigger() const { return single0_trigger == 1; };

        /**
         * Indicate whether a single1 trigger was registered for the event.
         *
         * @return Returns true if a single1 trigger was registered for the,
         *         false otherwise.
         */
        bool isSingle1Trigger() const { return single1_trigger == 1; };

        /**
         * Indicate whether the SVT bias was on during the event.
         *
         * @return Returns true if the bias was one, false otherwise.
         */
        bool isSvtBiasOn() const { return svt_bias_state == 1; };

        /**
         * Indicates whether the SVT was open or closed during an event.
         *
         * @return Returns true if the SVT was closed, false otherwise.
         */ 
        bool isSvtClosed() const { return svt_position_state == 1; }; 

        /**
         * Indicate whether the SVT latency was correct during an event.
         *
         * @return Returns true if the SVT latency was correct, false 
         *         otherwise.
         */
        bool isSvtLatencyGood() const { return svt_latency_state == 1; };

        /**
         * Indicates whether the event was affected by SVT burst noise.
         *
         * @return Returns true if the event has SVT burst noise, false 
         *         otherwise. 
         */
        bool hasSvtBurstModeNoise() const { return svt_burstmode_noise == 0; };

        /**
         * Indicates whether the SVT event headers had errors.
         *
         * @return Returns true if the SVT event headers had an error, 
         *         false otherwise.
         */
        bool hasSvtEventHeaderErrors() const { return svt_event_header_state == 0; }; 

        ClassDef(HpsEvent, 3);  

    public:

        /** Collection of beam spot constrained Moller candidates */
        TClonesArray* bsc_moller_candidates; //->
        /** Collection of beam spot constrained v0 candidates */
        TClonesArray* bsc_v0_candidates;     //->
        /** Collection of Ecal clusters */
        TClonesArray* ecal_clusters;         //->
        /** Collection of Ecal hits */
        TClonesArray* ecal_hits;             //->
        /** Collection of final state particles */
        TClonesArray* fs_particles;          //->
        /** Collection of GBL tracks */
        TClonesArray* gbl_tracks;            //-> 
        /** Collection of Monte Carlo particles */
        TClonesArray* mc_particles;          //->
        /** Collection of target constrained Moller candidates */
        TClonesArray* tc_moller_candidates;  //->
        /** Collection of target constrained v0 candidates */
        TClonesArray* tc_v0_candidates;      //->
        /** Collection of SVT tracks */ 
        TClonesArray* tracks;                //->
        /** Collection of SVT 3D hits */
        TClonesArray* svt_hits;              //->
        /** Collection of unconstrained Moller candidates */
        TClonesArray* uc_moller_candidates;  //->
        /** Collection of unconstrained v0 candidates */
        TClonesArray* uc_v0_candidates;      //->
        /** Collection of unconstrained v0 candidates in same side of detector. */
        TClonesArray* uc_vc_candidates;      //->
        /** Collection of other electrons */
        TClonesArray* other_electrons;          //->

        //-- Event information --//
        //-----------------------//

        /** The RF time */
        double rf_times[2]; 

        /** Event number. */
        int event_number;
        
        /** The event time. */
        long event_time;

        /** 
         * Flag indicating that a pair0 trigger was registered. It's 
         * set to 1 if it was registered or 0 if it wasn't.
         */
        int pair0_trigger;
        
        /** 
         * Flag indicating that a pair1 trigger was registered. It's 
         * set to 1 if it was registered or 0 if it wasn't.
         */
        int pair1_trigger;
        
        /** 
         * Flag indicating that a pulser (random) trigger was registered. It's 
         * set to 1 if it was registered or 0 if it wasn't.
         */
        int pulser_trigger;
        
        /** Run number */
        int run_number;
        
        /** 
         * Flag indicating that a singles0 trigger was registered. It's 
         * set to 1 if it was registered or 0 if it wasn't.
         */
        int single0_trigger;
        
        /** 
         * Flag indicating that a singles1 trigger was registered. It's 
         * set to 1 if it was registered or 0 if it wasn't.
         */
        int single1_trigger;
        
        /** 
         * Flag indicating the state of the SVT bias. It's set to 0 if the bias
         * was off or 1 if it was on.
         */
        int svt_bias_state;
        
        /**
         * Flag indicating whether the event was affected by SVT burst noise. 
         * It's set to 0 if the event saw burst noise or 1 if it was fine.
         */ 
        int svt_burstmode_noise; 

        /**
         * Flag indicating whether the SVT event headers had an error. It's 
         * set to 0 if the event headers had and error, of 1 if it was errorless. 
         */
        int svt_event_header_state; 

        /**
         * Flag indicating whether the SVT latency was correct during an event. 
         * It's set to 0 if the latency was incorrect, 1 if it was fine.
         */
        int svt_latency_state; 

        /** 
         * Flag indicating whether the SVT was open or closed.  It's set to 0 if 
         * the SVT was open or 1 if it was closed.
         */ 
        int svt_position_state;
 
        int n_tracks;
        int n_gbl_tracks;
        int n_svt_hits;
        int n_ecal_clusters;
        int n_ecal_hits;
        int n_fs_particles;
        int n_uc_v0_candidates;
        int n_uc_vc_candidates; 
        int n_uc_moller_candidates;
        int n_bsc_v0_candidates;
        int n_bsc_moller_candidates;
        int n_tc_v0_candidates;
        int n_tc_moller_candidates; 
        int n_mc_particles;
        int n_other_electrons;
};

#endif



///
/// This class reads the LCIO files and parses them for the MiniDst root files.
///
/// In principle, you could also derive a new class from this class and write analysis code to directly
/// analyze the LCIO file without writing out the MiniDST root file. To do so, your Start() function would need
/// to disable the md_output_tree by setting it to nullptr, so it does not write out...
///
#include "LcioReader.h"

void LcioReader::Start(){
    MiniDst::Start();
};
long LcioReader::Run(int nevent){

    for(string file: input_file_list){
        data_type_is_known=false;
        is_2016_data = false;
        is_2019_data = false;
        lcio_reader->open(file);

        while( (lcio_event = lcio_reader->readNextEvent())){

            clear();  // Clear all the vectors that contain data, so we can push_back on them again.

            /////////////////////////////////////////////////////////////////////////////////////////////
            ///
            /// Event Header
            ///
            /////////////////////////////////////////////////////////////////////////////////////////////

            run_number = lcio_event->getRunNumber();

            if( !data_type_is_known){ // Determine the data type by looking at the collections
                if(run_number < 9000){ // The 2015 or 2016 engineering runs.
                    is_2016_data = true;
                    if(md_Debug & kInfo) cout << "LCIO -> This is 2015/2016 data. \n";
                }else if( run_number <= 10750 ){ // 2019 physics run
                    is_2019_data = true;
                    if(md_Debug & kInfo) cout << "LCIO -> This is 2019 data. \n";
                }

                col_names = lcio_event->getCollectionNames();
                if(md_Debug & kDebug_L1){
                    cout << "LCIO Collections in " << file << ":\n";
                    for(string s: *col_names){
                        cout << s << endl;
                    }
                }
                for(string s: *col_names){
// In the current reco of data, 2019 LCIO files have both the TriggerBank and the TSBank, so this cannot
// help us determine what type of data we are looking at.
//                    if( s == "TriggerBank"){
//                        is_2016_data = true;
//                        if(md_Debug & kInfo) cout << "LCIO -> This is 2015/2016 data. \n";
//                    }
//                    if( s == "TSBank"){
//                        is_2019_data = true;
//                        is_2016_data = false;
//                        if(md_Debug & kInfo) cout << "LCIO -> This is 2019 data. \n";
//                    }
                    if( s == "MCParticles"){
                        is_MC_data = true;
                        if(md_Debug & kInfo) cout << "LCIO -> This is Monte Carlo data. \n";
                    }
                }
                if( is_2016_data && is_2019_data) cout << "WOA - a file that is both 2016 and 2019 data!!!\n";
                data_type_is_known = true;
            }

            event_number = lcio_event->getEventNumber();
            if( (++evt_count)%Counter_Freq == 0) {
                printf("i: %10lu   event: %10d  run: %5d\n", evt_count, event_number, run_number);
            }
            time_stamp = lcio_event->getTimeStamp();

            // The following quantities from the LCIO header were set for 2016 data, but not consistently.
            // Currently for 2019 data they are all 0. These should actually be booleans not ints.
            int svt_bias_good = lcio_event->getParameters().getIntVal("svt_bias_good");
            int svt_burstmode_noise_good = lcio_event->getParameters().getIntVal("svt_burstmode_noise_good");
            int svt_latency_good = lcio_event->getParameters().getIntVal("svt_latency_good");
            int svt_position_good = lcio_event->getParameters().getIntVal("svt_position_good");
            int svt_event_header_good = lcio_event->getParameters().getIntVal("svt_event_header_good");
            // Pack them into a single unsigned int.
            svt_status = svt_bias_good +
                    (svt_burstmode_noise_good<<1) +
                    (svt_latency_good<<2) +
                    (svt_position_good<<3) +
                    (svt_event_header_good<<4);

            // Get the LCIO GenericObject collection containing the RF times
            EVENT::LCCollection* rf_hits
                    = static_cast<EVENT::LCCollection*>(lcio_event->getCollection("RFHits"));

            if(rf_hits->getNumberOfElements() != 1){
                if(md_Debug & kDebug_L1 )
                    cout << "Issue with rf_hits for event " << event_number << " == No rf hits found. \n";
                rf_time1 = -99.;
                rf_time2 = -99.;
            }else {
                EVENT::LCGenericObject *rf_times = static_cast<EVENT::LCGenericObject *>(rf_hits->getElementAt(0));
                if (rf_times->getNDouble() != 2) {
                    cout << "ERROR in LCIO - wrong number of hits. \n";
                } else {
                    rf_time1 = rf_times->getDoubleVal(0);
                    rf_time2 = rf_times->getDoubleVal(1);
                }
            }
            //
            // For trigger bit parsing, see EvioTool::TSBank.h
            // Note that LCIO has the type of bank in location 0, so compared to the EVIO data
            // all the words are shifted +1.
            //
            if(is_2016_data) {
                EVENT::LCCollection *lcio_triggers
                        = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("TriggerBank"));
                int n_trigger_banks = lcio_triggers->getNumberOfElements();
                if( n_trigger_banks != 3 ){
                    // There should be 3 banks from EVIO: 0xe10c (57612) 0xe10a (57610) and 0xe10f (57615)
                    cout << "Holy cow, there were inconsistent number of triggers in this event!";
                };
                for(int i=0; i< n_trigger_banks; ++i) {
                    EVENT::LCGenericObject *lcio_trig
                            = static_cast<EVENT::LCGenericObject *>(lcio_triggers->getElementAt(i));
                    int nvalues = lcio_trig->getNInt();
                    int test_val = lcio_trig->getIntVal(0);
                    if(md_Debug & kDebug_L1) printf("Trigger int 0 = 0x%04x (%6d) has %2d int values\n",
                            test_val,test_val,nvalues);
                    if( test_val == 57610 ){ // Trigger bits are here.
                        unsigned int trigger_bits = lcio_trig->getIntVal(1);
                        // We need to re-encode the bits so they are consistent with 2019 data. See MiniDst.h
                        TriggerBits_int_t trig_code{0};
                        trig_code.bits.Single_0_Top =  (trigger_bits & (1<<24)); // Note: sets to 0 or 1.
                        trig_code.bits.Single_0_Bot = trig_code.bits.Single_0_Top;
                        trig_code.bits.Single_1_Top =  (trigger_bits & (1<<25));
                        trig_code.bits.Single_1_Bot = trig_code.bits.Single_1_Top;
                        trig_code.bits.Pair_0  =  (trigger_bits & (1<<26));
                        trig_code.bits.Pair_1  =  (trigger_bits & (1<<27));
                        trig_code.bits.Pulser  =  (trigger_bits & (1<<29));
                        trigger = trig_code.intval; // Finally, set the trigger to be stored.

                        time_stamp = lcio_trig->getIntVal(4) & 0xFFFF; // Time stamp high word.
                        time_stamp <<= 32;
                        time_stamp += lcio_trig->getIntVal(3); // Time stamp low word.
                    }
                } // for each trigger_bank
            } // is_2016_data
            if(is_2019_data){
                EVENT::LCCollection *tsbank_data
                        = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("TSBank"));
                EVENT::LCGenericObject *lcio_ts_bank
                        = static_cast<EVENT::LCGenericObject*>(tsbank_data->getElementAt(0));


                time_stamp = lcio_ts_bank->getIntVal(4) & 0xFFFF; // Time stamp high word.
                time_stamp <<= 32;
                time_stamp += lcio_ts_bank->getIntVal(3); // Time stamp low word.

                unsigned long trigger_number = lcio_ts_bank->getIntVal(4) & 0xFFFF0000;
                trigger_number <<= 16;
                trigger_number += lcio_ts_bank->getIntVal(2);

                // Prescaled is in word 5, Ext is in word 6
                trigger = lcio_ts_bank->getIntVal(5);
                ext_trigger = lcio_ts_bank->getIntVal(6);
                // unsigned int no_prescale_trigger = ext_trigger = lcio_ts_bank->getIntVal(7);

//                EVENT::LCCollection *vtpbank_data
//                        = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("VTPBank"));
//                EVENT::LCGenericObject* vtpbank
//                        = static_cast<EVENT::LCGenericObject*>(vtpbank_data->getElementAt(0));

// Todo: Add the VTP Parsing when needed.
            } // End of 2019 specific header parsing.

            ////////////////////////////////////////////////////////////////////////////////////////////////
            ///
            /// Ecal
            ///
            /// Parse the "EcalCalHits" and "EcalClustersCorr" LCIO collections.
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////////

            map<IMPL::CalorimeterHitImpl*, int> hit_pointer_to_index; // Store a map to speed up hit to cluster linking.
            UTIL::BitField64 ecal_hit_field_decoder("system:6,layer:2,ix:-8,iy:-6");

            if(write_ecal_hits || write_ecal_cluster){ // Need the hits to fill the clusters!
                EVENT::LCCollection* ecal_hits =
                        static_cast<EVENT::LCCollection*>(lcio_event->getCollection("EcalCalHits"));

                for (int ihit=0; ihit < ecal_hits->getNumberOfElements(); ++ihit) {
                    IMPL::CalorimeterHitImpl* lcio_hit
                            = static_cast<IMPL::CalorimeterHitImpl*>(ecal_hits->getElementAt(ihit));

                    hit_pointer_to_index[lcio_hit] = ihit;
                    // Gets the CellID which identifies the crystal.
                    // int id0 = lcio_hit->getCellID0();
                    // 0.1 ns resolution is sufficient to distinguish any 2 hits on the same crystal.
                    // int id1 = static_cast<int>(10.0*lcio_hit->getTime());

                    if(write_ecal_hits) {  // Only store the results if we actually want them.
                        ecal_hit_energy.push_back(lcio_hit->getEnergy());
                        ecal_hit_time.push_back(lcio_hit->getTime());

                        Long64_t value = Long64_t(lcio_hit->getCellID0() & 0xffffffff) |
                                         (Long64_t(lcio_hit->getCellID1()) << 32);
                        ecal_hit_field_decoder.setValue(value);

                        ecal_hit_index_x.push_back(ecal_hit_field_decoder["ix"]);
                        ecal_hit_index_x.push_back(ecal_hit_field_decoder["iy"]);
                    }
                }
            }

            if(write_ecal_cluster){
                EVENT::LCCollection* clusters =
                        static_cast<EVENT::LCCollection*>(lcio_event->getCollection("EcalClustersCorr"));
                for(int iclus=0;  iclus< clusters->getNumberOfElements(); ++iclus){
                    IMPL::ClusterImpl* lcio_clus = static_cast<IMPL::ClusterImpl*>(clusters->getElementAt(iclus));

                    ecal_cluster_energy.push_back(lcio_clus->getEnergy());
                    const float *position = lcio_clus->getPosition();
                    ecal_cluster_x.push_back(position[0]);
                    ecal_cluster_y.push_back(position[1]);
                    ecal_cluster_z.push_back(position[2]);

                    // Link the hits to the cluster.
                    EVENT::CalorimeterHitVec clus_hits = lcio_clus->getCalorimeterHits();
                    double seed_energy{-99.};
                    double seed_time{-99};
                    int    seed_index{-99};
                    Long64_t seed_cellid0{-99};
                    int    seed_ix{-99};
                    int    seed_iy{-99};
                    vector<int> clus_hit_indexes;
                    for( int j_hit=0; j_hit< clus_hits.size(); ++j_hit){
                        IMPL::CalorimeterHitImpl* hit = static_cast<IMPL::CalorimeterHitImpl*>(clus_hits[j_hit]);
                        if( hit->getEnergy() >  seed_energy ){
                            seed_energy = hit->getEnergy();
                            seed_time = hit->getTime();
                            seed_index = hit_pointer_to_index[hit];
                            seed_cellid0 = Long64_t(hit->getCellID0() & 0xffffffff) |
                                             (Long64_t(hit->getCellID1()) << 32);

                        }
                        int hit_index = hit_pointer_to_index[hit];
                        clus_hit_indexes.push_back(hit_index);
                    }
                    ecal_cluster_time.push_back(seed_time);        // The cluster time = seed hit time.
                    ecal_cluster_hits.push_back(clus_hit_indexes);
                    ecal_cluster_seed_index.push_back(seed_index);
                    ecal_cluster_seed_energy.push_back(seed_energy);
                    ecal_hit_field_decoder.setValue(seed_cellid0);
                    ecal_cluster_seed_ix.push_back(ecal_hit_field_decoder["ix"]);
                    ecal_cluster_seed_iy.push_back(ecal_hit_field_decoder["iy"]);
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////
            ///
            /// SVT
            ///
            /// Parse the "RotatedHelicalTrackHits" and "GBLTracks" and "MatchedTracks" LCIO collections.
            ///
            /// Todo: Add the raw tracker hits if needed.
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////////

            std::map<IMPL::TrackerHitImpl*, int> svt_hit_map;
            if(write_svt_hits || write_tracks){  // Need the svt hitmap either way.
                UTIL::BitField64 raw_svt_hit_decoder("system:6,barrel:3,layer:4,module:12,sensor:1,side:32:-2,strip:12");
                EVENT::LCCollection* tracker_hits = lcio_event->getCollection("RotatedHelicalTrackHits");
                for (int i_svt_hit = 0; i_svt_hit< tracker_hits->getNumberOfElements(); ++i_svt_hit) {
                    IMPL::TrackerHitImpl* lcio_svt_hit =
                            static_cast<IMPL::TrackerHitImpl*>(tracker_hits->getElementAt(i_svt_hit));
                    svt_hit_map[lcio_svt_hit] = i_svt_hit;

                    if(write_svt_hits){
                        svt_hit_time.push_back(lcio_svt_hit->getTime());
                        const double *pos = lcio_svt_hit->getPosition();
                        svt_hit_x.push_back(pos[0]);
                        svt_hit_y.push_back(pos[1]);
                        svt_hit_z.push_back(pos[2]);
                        svt_hit_edep.push_back(lcio_svt_hit->getEDep()); // Not in 2016 data.

                        EVENT::LCObjectVec raw_hits = lcio_svt_hit->getRawHits();
                        EVENT::TrackerRawData *lcio_raw_hit = static_cast<EVENT::TrackerRawData *>(raw_hits.at(0));
                        ULong64_t value = EVENT::long64(lcio_raw_hit->getCellID0() & 0xffffffff) |
                                          (EVENT::long64(lcio_raw_hit->getCellID1()) << 32);
                        raw_svt_hit_decoder.setValue(value);

//                            int system = raw_svt_hit_decoder["system"];
//                            int barrel = raw_svt_hit_decoder["barrel"];
                        int layer = raw_svt_hit_decoder["layer"];
//                            int module = raw_svt_hit_decoder["module"];
//                            int sensor  = raw_svt_hit_decoder["sensor"];
//                            int side  = raw_svt_hit_decoder["side"];
//                            int strip = raw_svt_hit_decoder["strip"];
                        svt_hit_layer.push_back(layer);
                    }
                }
            }

            if( write_tracks) {
                EVENT::LCCollection* gbl_tracks = lcio_event->getCollection("GBLTracks");
                EVENT::LCCollection* matched_tracks{nullptr};
                int n_matched_tracks=0;
                if(!write_only_gbl_tracks) {
                    matched_tracks = lcio_event->getCollection("MatchedTracks");
                    n_matched_tracks = matched_tracks->getNumberOfElements();
                }

                map<EVENT::Track*, int> gbl_track_to_index_map;
                map<EVENT::Track*, int> matched_track_to_index_map;

                int n_gbl_track = gbl_tracks->getNumberOfElements();
                int n_total_tracks = n_gbl_track + n_matched_tracks;
                EVENT::Track* lcio_track{nullptr};
                for(int track_number =0; track_number < n_total_tracks; ++track_number){
                    // We go through all the GBL tracks first, and if asked for, then go through all the
                    // matched tracks.
                    if(track_number < n_gbl_track) {
                        lcio_track = static_cast<EVENT::Track *>(gbl_tracks->getElementAt(track_number));
                        gbl_track_to_index_map[lcio_track] = track_number;
                        track_gbl_ref.push_back(track_number); // GBL Track points to itself.
                        track_ref.push_back(-99); // Pointer to seed track is resolved later.
                    }else{
                        int matched_track_number = track_number - n_gbl_track;
                        lcio_track = static_cast<EVENT::Track *>(matched_tracks->getElementAt(matched_track_number));
                        matched_track_to_index_map[lcio_track] = track_number;
                        track_gbl_ref.push_back(-99);  // Seed track needs to be resolved later.
                        track_ref.push_back(track_number); // Seed track points to itself.
                    }
                    track_particle.push_back(-99); // Pointer to particle is resolved *way* later.

                    track_d0.push_back(lcio_track->getD0());
                    track_chi2.push_back(lcio_track->getChi2());
                    track_omega.push_back(lcio_track->getOmega());
                    track_phi0.push_back(lcio_track->getPhi());
                    track_tan_lambda.push_back(lcio_track->getTanLambda());
                    track_z0.push_back(lcio_track->getZ0());
                    track_type.push_back(lcio_track->getType());

                    const EVENT::TrackState* track_state
                            = lcio_track->getTrackState(EVENT::TrackState::AtCalorimeter);
                    if(track_state){
                        const float* ecal_pos = track_state->getReferencePoint();
                        track_x_at_ecal.push_back(ecal_pos[1]);  // Note: Un-rotate from tracking coordinate system.
                        track_y_at_ecal.push_back(ecal_pos[2]);
                        track_z_at_ecal.push_back(ecal_pos[0]);
                    }else{
                        track_x_at_ecal.push_back(-99.);
                        track_y_at_ecal.push_back(-99.);
                        track_z_at_ecal.push_back(-99.);
                    }
                    EVENT::LCCollection* gbl_kink_data_rel =
                            static_cast<EVENT::LCCollection*>(lcio_event->getCollection("GBLKinkDataRelations"));

                    unique_ptr<UTIL::LCRelationNavigator> gbl_kink_data_nav =
                            make_unique<UTIL::LCRelationNavigator>(gbl_kink_data_rel);

                    // Get the list of GBLKinkData associated with the LCIO Track
                    EVENT::LCObjectVec gbl_kink_data_list
                            = gbl_kink_data_nav->getRelatedFromObjects(lcio_track);

                    IMPL::LCGenericObjectImpl *gbl_kink_datum{nullptr};

                    if(gbl_kink_data_list.size() == 1) {
                        vector<double> lambdas;
                        vector<double> kinks;
                        // Get the list GBLKinkData GenericObject associated with the LCIO Track
                        gbl_kink_datum = static_cast<IMPL::LCGenericObjectImpl *>(gbl_kink_data_list.at(0));
                        for (int ikink = 0; ikink < gbl_kink_datum->getNFloat(); ++ikink) {
                            lambdas.push_back(gbl_kink_datum->getFloatVal(ikink));
                            kinks.push_back(gbl_kink_datum->getDoubleVal(ikink));
                        }
                        track_lambda_kinks.push_back(lambdas);
                        track_phi_kinks.push_back(kinks);
                    }else{
                        track_lambda_kinks.push_back(vector<double>(0));
                        track_phi_kinks.push_back(vector<double>(0));
                    }

                    // Get the collection of LCRelations between track data variables
                    // (TrackData) and the corresponding track.
                    EVENT::LCCollection* track_data = static_cast<EVENT::LCCollection*>(
                            lcio_event->getCollection("TrackDataRelations"));

                    // Instantiate an LCRelation navigator which will allow faster access
                    // to TrackData objects
                    UTIL::LCRelationNavigator* track_data_nav
                            = new UTIL::LCRelationNavigator(track_data);

                    // Get the list of TrackData associated with the LCIO Track
                    EVENT::LCObjectVec track_data_list = track_data_nav->getRelatedFromObjects(lcio_track);

                    vector<double> iso_values(14,-99.);
                    if(track_data_list.size() ==1) {
                        // There should always be one and only one for GBL tracks, 0 for matched tracks.
                        IMPL::LCGenericObjectImpl *track_info =
                                static_cast<IMPL::LCGenericObjectImpl *>(track_data_list.at(0));

                        if (track_info->getNDouble() < 12 || track_info->getNDouble() > 14 || track_info->getNFloat() != 1
                            || track_info->getNInt() != 1) {
                            cout << "Dude, you got messed up SVT tracking data!\n";
                        }

                        for (int i_iso = 0; i_iso < track_info->getNDouble(); ++i_iso) {
                            iso_values[i_iso] = track_info->getDoubleVal(i_iso);
                        }
                        track_time.push_back(track_info->getFloatVal(0));
                        track_volume.push_back(track_info->getIntVal(0));
                    }else{
                        track_time.push_back(-99.);
                        track_volume.push_back(-1);
                    }
                    track_isolation.push_back(iso_values);
                    delete track_data_nav;

                    // Get the collection of 3D hits associated with a LCIO Track
                    EVENT::TrackerHitVec gbl_tracker_hits = lcio_track->getTrackerHits();

                    track_n_hits.push_back(gbl_tracker_hits.size());
                    vector<int> svt_hits;
                    svt_hits.reserve(gbl_tracker_hits.size());
                    for(int i_trk=0; i_trk < gbl_tracker_hits.size(); ++i_trk){
                        IMPL::TrackerHitImpl* trk_hit =
                                static_cast<IMPL::TrackerHitImpl*>(gbl_tracker_hits[i_trk]);
                        int hit_index = svt_hit_map[trk_hit];
                        svt_hits.push_back(hit_index);
                    }
                    track_svt_hits.push_back(svt_hits);

                    // TODO: It would be better to permit vector<vector<float>> and store that
                    //  instead of converting to double.
                    const vector<float> cov_matrix_f = lcio_track->getCovMatrix();
                    vector<double> cov_matrix_d(cov_matrix_f.begin(),cov_matrix_f.end());
                    track_covmatrix.push_back(cov_matrix_d);
                } // End loop over tracks.

                // Get the collection of LCRelations between seed tracks and a GBL tracks.
                EVENT::LCCollection* seed_to_gbl_relations =
                        static_cast<EVENT::LCCollection*>(lcio_event->getCollection("MatchedToGBLTrackRelations"));

                // Instantiate an LCRelation navigator which will allow faster access
                // to the seed to GBL LCRelations
                UTIL::LCRelationNavigator* seed_to_gbl_relations_nav =
                        new UTIL::LCRelationNavigator(seed_to_gbl_relations);

                for( auto const &[gbl_track, gbl_track_index] : gbl_track_to_index_map){
                    EVENT::LCObjectVec seed_to_gbl_list
                            = seed_to_gbl_relations_nav->getRelatedFromObjects(gbl_track);
                    if (seed_to_gbl_list.size() != 1) {
                        cout << "Woops, I expected only one seed track for a gbl track.\n";
                    }else{
                        EVENT::Track* seed_track = static_cast<EVENT::Track*>(seed_to_gbl_list.at(0));
                        int seed_index = matched_track_to_index_map[seed_track];
                        track_ref[gbl_track_index] = seed_index;
                        if(track_gbl_ref.size() > seed_index){
                            int debug_copy_seed_index = seed_index;
                            int debug_copy_gbl_track_index = gbl_track_index;
                            track_gbl_ref[seed_index] = gbl_track_index;
                        }
                    }
                }
            } // End if(write_tracks)

            ////////////////////////////////////////////////////////////////////////////////////////////////
            ///
            /// Particles
            ///
            /// Two types of particles: 1 or 0 tracks -> single particle,  2 tracks -> vertexes.
            ///
            /// Collections are:
            /// particle_types_single{FINAL_STATE_PARTICLE, OTHER_ELECTRONS};
            /// particle_types_double{TC_V0_CANDIDATE, UC_VC_CANDIDATE};
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////////

            // Single particles.
            for( auto type : particle_types_single){
                string collection_name = Type_to_Collection[type];
                EVENT::LCCollection* particles = lcio_event->getCollection(collection_name.c_str());
                for(int i_part=0; i_part < particles->getNumberOfElements(); ++i_part){
                    EVENT::ReconstructedParticle *lcio_part =
                            static_cast<EVENT::ReconstructedParticle*>(particles->getElementAt(i_part));
                    Fill_Part_From_LCIO(&part, lcio_part);
                    part.type.push_back(type);
                }

            }

            // Double Track Particles, i.e.  Constrained V0 Candidates.
            for(int type: particle_types_double) {
                string collection_name = Type_to_Collection[type];
                EVENT::LCCollection* particles = lcio_event->getCollection(collection_name.c_str());
                for(int i_part=0; i_part < particles->getNumberOfElements(); ++i_part){
                    EVENT::ReconstructedParticle *lcio_part =
                            static_cast<EVENT::ReconstructedParticle*>(particles->getElementAt(i_part));
                    Fill_Part_From_LCIO(&v0, lcio_part);
                    EVENT::Vertex *vertex = lcio_part->getStartVertex();
                    const float *vertex_pos = vertex->getPosition();
                    v0.vertex_x.push_back(vertex_pos[0]);
                    v0.vertex_y.push_back(vertex_pos[1]);
                    v0.vertex_z.push_back(vertex_pos[2]);
                    v0.vertex_chi2.push_back(vertex->getChi2());
                    double probability = vertex->getProbability();
                    v0.vertex_prob.push_back(probability);
                    v0.vertex_type.push_back(type);
//                    const std::string &xxx= vertex->getAlgorithmType(); // TargetConstrained/Unconstrained/...
//                    cout << "Vertex of type: " << xxx << "\n";

                    EVENT::ReconstructedParticle *sub_part = vertex->getAssociatedParticle();
                    const vector<EVENT::Track *> &tracks = lcio_part->getTracks();
                    const EVENT::ClusterVec &cluster_vec = lcio_part->getClusters();
                    const EVENT::ReconstructedParticleVec &daughters = lcio_part->getParticles();
                    if(daughters.size() != 2){
                        cout << "LcioReader:: Vertex should have 2 and only 2 daughters, but this one has:" <<
                        daughters.size() << " !\n";
                    }
                    EVENT::ReconstructedParticle *electron;
                    EVENT::ReconstructedParticle *positron;
                    if( daughters[0]->getCharge() < 0 ){ // Daughter 0 is an electron.
                        electron = daughters[0];
                        positron = daughters[1];  // NOTE: for Mollers they are both electron. Whatever.
                    }else{
                        electron = daughters[1];
                        positron = daughters[0];
                    }
                    Fill_SubPart_From_LCIO(&v0.em, electron, lcio_event);
                    Fill_SubPart_From_LCIO(&v0.ep, positron, lcio_event);
                }
            }

            if(md_output_tree){
                if( md_abort_tree_fill){
                    cout << "Bad Event -- Not filling TTree \n";
                    md_abort_tree_fill = false;
                }else {
                    md_output_tree->Fill();
                }
            }
        } // End event loop.

        lcio_reader->close();
    }
    return(0);
};

void LcioReader::Fill_SubPart_From_LCIO(Sub_Particle_t *sub,EVENT::ReconstructedParticle *daughter,
                                        EVENT::LCEvent *lcio_evnet){
    // Fill the Sub_Particle_t structure from the LCIO daughter particle.

    // Find the particle assoc with this sub particle from FINAL_STATE_PARTICLES.
    EVENT::LCCollection* particles = lcio_event->getCollection(Type_to_Collection[FINAL_STATE_PARTICLE].c_str());
    bool found_daughter_particle = false;
    for(int i_part=0; i_part < particles->getNumberOfElements(); ++i_part) {
        EVENT::ReconstructedParticle *lcio_part =
                static_cast<EVENT::ReconstructedParticle *>(particles->getElementAt(i_part));
        if(lcio_part == daughter){
            sub->part.push_back(i_part);
            sub->good_pid.push_back(daughter->getGoodnessOfPID());
            const double *momentum = daughter->getMomentum();
            double daughter_p = sqrt( momentum[0]*momentum[0]+ momentum[1]*momentum[1]+momentum[2]*momentum[2]);
            sub->p.push_back(daughter_p);

            found_daughter_particle = true;
            break;
        }
    }
    if(!found_daughter_particle){
        cout << "We did not find the daughter particle for a vertex in Run " <<
            run_number << "::" << event_number << "\n";
        md_abort_tree_fill = true;
        return;
    }

    const EVENT::TrackVec &tracks = daughter->getTracks();       // Get the tracks associated with this daughter.
    if(tracks.size() !=1 ){
        /// A sub particle should *always* have a track.
        cout << "Vertex sub particle with " << tracks.size() << " tracks! Expected 1. -- Run " <<
            run_number << "::" << event_number << "\n";
        md_abort_tree_fill = true;
        return;
    }

    EVENT::Track *track = tracks[0];
    EVENT::LCCollection* gbl_tracks = lcio_event->getCollection("GBLTracks");
    bool track_found = false;
    for(int i_trk=0; i_trk < gbl_tracks->getNumberOfElements(); ++i_trk){
        EVENT::Track *lcio_track = static_cast<EVENT::Track *>(gbl_tracks->getElementAt(i_trk));
        if(lcio_track == track){
            // We found the matching track.
            // Double check that Chi2 is the same.
            if(abs(lcio_track->getChi2()- track_chi2[i_trk])>1.0E-6){
                cout << "Not sure we got the right track: chi2_1 = " << lcio_track->getChi2() <<
                "  chi2_2 = " << track_chi2[i_trk] << "\n";
            }
            sub->track.push_back(i_trk);  // Store the index.
            sub->track_time.push_back(track_time[i_trk]);
            sub->track_nhit.push_back(track_n_hits[i_trk]);
            sub->chi2.push_back(lcio_track->getChi2());
            sub->pos_ecal_x.push_back(track_x_at_ecal[i_trk]);
            sub->pos_ecal_y.push_back(track_y_at_ecal[i_trk]);
            track_found = true;
            break;
        }
    }
    if(!track_found){
        cout << "ERROR -- Fill_SubPart_From_LCIO -- Did not find associated track!!!!\n";
        md_abort_tree_fill=true;
        return;
    }


    const EVENT::ClusterVec &clusters = daughter->getClusters();
    if(clusters.size() > 1 ){
        // A sub particle can have 1 (normal) or 0 (ECal hole) clusters.
        cout << "Vertex sub particle with " << clusters.size() << " clusters! Expected 1 or 0 -- Run "<<
        run_number << "::" << event_number << "\n";
        md_abort_tree_fill = true;
        return;
    }

    if(clusters.size() == 1) {
        IMPL::ClusterImpl *clus = static_cast<IMPL::ClusterImpl *>(clusters[0]);

        EVENT::LCCollection *lcio_clusters =
                static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalClustersCorr"));
        bool cluster_found = false;
        for (int iclus = 0; iclus < lcio_clusters->getNumberOfElements(); ++iclus) {
            IMPL::ClusterImpl *lcio_clus = static_cast<IMPL::ClusterImpl *>(lcio_clusters->getElementAt(iclus));
            if (clus == lcio_clus) {
                sub->clus.push_back(iclus);
                if( abs(ecal_cluster_energy[iclus] - lcio_clus->getEnergy())> 1.E-6 ){
                    cout << "Not sure we found the correct cluster. Clus1_E = " << ecal_cluster_energy[iclus] <<
                    " Clus2_E = " << lcio_clus->getEnergy() << "\n";
                }
                sub->clus_energy.push_back(lcio_clus->getEnergy());
                sub->clus_time.push_back(ecal_cluster_time[iclus]);
                sub->clus_pos_x.push_back(ecal_cluster_x[iclus]);
                sub->clus_pos_y.push_back(ecal_cluster_y[iclus]);
                sub->clus_ix.push_back(ecal_cluster_seed_ix[iclus]);
                sub->clus_iy.push_back(ecal_cluster_seed_iy[iclus]);
                cluster_found = true;
                break;
            }
        }
        if (!cluster_found) {
            cout << "ERROR -- Fill_SubPart_From_LCIO -- Did not find associated cluster!!!!\n";
            md_abort_tree_fill=true;
            return;
        }
    }else{ // Clusterless track. Fill cluster stuff with -99.

    }
}

void LcioReader::Fill_Part_From_LCIO(Basic_Particle_t *bp, EVENT::ReconstructedParticle *lcio_part) {
    // Fill the Basic_Particle_t structure with the information from a ReconstructedParticle.
    bp->lcio_type.push_back(lcio_part->getType());
    bp->energy.push_back(lcio_part->getEnergy());
    bp->mass.push_back(lcio_part->getMass());
    EVENT::ParticleID *lcio_part_id = lcio_part->getParticleIDUsed();
    if(lcio_part_id){
        int pdg = lcio_part_id->getPDG();
        bp->pdg.push_back(pdg);
    }else{
        bp->pdg.push_back(0);
    }

    bp->charge.push_back(lcio_part->getCharge());
    bp->goodness_of_pid.push_back(lcio_part->getGoodnessOfPID());
    const double *momentum = lcio_part->getMomentum();
    bp->px.push_back(momentum[0]);
    bp->py.push_back(momentum[1]);
    bp->pz.push_back(momentum[2]);
}

void LcioReader::End(){
    MiniDst::End();
};

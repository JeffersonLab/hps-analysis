///
/// This class reads the LCIO files and parses them for the MiniDst root files.
///
/// In principle, you could also derive a new class from this class and write analysis code to directly
/// analyze the LCIO file without writing out the MiniDST root file. To do so, your Start() function would need
/// to disable the md_output_tree by setting it to nullptr, so it does not write out...
///
#include "LcioReader.h"

LcioReader::LcioReader(const string &input_file, const int debug_level) {
   if(md_Debug > 0) md_Debug = debug_level;
   cout << "LcioReader Debug level is " << md_Debug << std::endl;
   if(!input_file.empty()) input_files.push_back(input_file);
};

LcioReader::LcioReader(const vector<string> &infile_list, const int debug_level){
   md_Debug = debug_level;
   if(md_Debug > 0) cout << "LcioReader Debug level is " << md_Debug << std::endl;
   for(auto f : infile_list){
      input_files.push_back(f);
   }

#ifdef DEBUG
   {
        for(auto file: input_files){
            std::cout << "File: " << file << std::endl;
        }
    }
#endif
};

void LcioReader::Start(){
   if( md_Debug>0 ) {
      printf("LCIO READER version " LCIOReader__Version__ "\n");
   }
   // Slightly "expensive", but it is really nice to know ahead of time if we need MCParticle in the DST.
   // Also, this lets us remove any items from the branches that are not in the LCIO file.
   if(input_files.empty()){
      cout << "No input files. Cannot initialize. \n";
      // TODO: We should throw an error here, or something.
      return;
   }

   if(kf_has_not_postscript){
      for(int i=0; i<Type_to_Collection.size(); ++i){
         if(Type_to_Collection[i].find("_KF") != string::npos){
            Type_to_Collection[i] = Type_to_Collection[i].substr(0, Type_to_Collection[i].size()-3);
         }
      }
      for(int i=0; i<Type_to_VertexCollection.size(); ++i){
         if(Type_to_VertexCollection[i].find("_KF") != string::npos){
            Type_to_VertexCollection[i] = Type_to_VertexCollection[i].substr(0, Type_to_VertexCollection[i].size()-3);
         }
      }
   }

   if(gbl_has_no_postscript){
      for(auto &name: Type_to_Collection){
         if(name.find("_GBL") != string::npos){
            name = name.substr(0, name.size()-4);
         }
      }
      for(auto &name: Type_to_VertexCollection){
         if(name.find("_GBL") != string::npos){
            name = name.substr(0, name.size()-4);
         }
      }
   }

   lcio_reader->open(input_files[0]);
   lcio_event =lcio_reader->readNextEvent();
   run_number = lcio_event->getRunNumber();
   SetupLcioDataType();
   lcio_reader->close();

   MiniDst::Start();
}

void LcioReader::Clear(){
   /// Clear the event storage.
   // Make sure you also call the "super"

   MiniDst::Clear();

   // Clear all the maps
   ecal_hit_to_index_map.clear();
   ecal_id0_to_hit_index.clear();
   ecal_cluster_to_index_map.clear();
   svt_hit_to_index_map.clear();
   svt_raw_hit_to_index_map.clear();
   kf_track_to_index_map.clear();
   gbl_track_to_index_map.clear();
   matched_track_to_index_map.clear();
   any_track_to_index_map.clear();
   any_particle_to_index_map.clear();
}

void LcioReader::SetupLcioDataType(){
   // Read the LCIO file and determine what capabilities it has.
   // Set appropriate flags accordingly.

   if (!data_type_is_known) { // Determine the data type by looking at the collections
      if (md_Debug & kDebug_L1) cout << "Setting up the LCIO data. \n";

      run_number = lcio_event->getRunNumber();
      if (run_number == 0){  // Some form of MC data, probably SLIC output (?)
         is_MC_data = true;
      }
      else if (run_number < 9000) { // The 2015 or 2016 engineering runs.
         is_2016_data = true;
         if(magnetic_field < 0.0001) magnetic_field = 0.523400;  // Field strength in Tesla for 2016 run.
         if (md_Debug & kDebug_Info) cout << "LCIO -> This is 2015/2016 data. Field set to " << magnetic_field << "T.\n";
      } else if (run_number <= 10750) { // 2019 physics run
         is_2019_data = true;
         if (md_Debug & kDebug_Info) cout << "LCIO -> This is 2019 data. \n";
      } else if (run_number > 10750) { // 2019 physics run
         is_2019_data = true;             // 2021 data behaves the same as 2019 data. (?)
         if (md_Debug & kDebug_Info) cout << "LCIO -> This is 2021 data. \n";
      }

      col_names = lcio_event->getCollectionNames();
      if (md_Debug & kDebug_L1) {
         cout << "LCIO Collections found:\n";
         for (string s: *col_names) {
            cout << s << endl;
         }
      }
// In the current reco of data, 2019 LCIO files have both the TriggerBank and the TSBank, so this cannot
// help us determine what type of data we are looking at.
//                    if( s == "TriggerBank"){
//                        is_2016_data = true;
//                        if(md_Debug & kDebug_Info) cout << "LCIO -> This is 2015/2016 data. \n";
//                    }
//                    if( s == "TSBank"){
//                        is_2019_data = true;
//                        is_2016_data = false;
//                        if(md_Debug & kDebug_Info) cout << "LCIO -> This is 2019 data. \n";
//                    }
      if (has_collection("MCParticle")) {
         if (md_Debug & kDebug_Info) cout << "LCIO -> This is Monte Carlo data. \n";
         if(use_mc_particles) is_MC_data = true;
         else{
            cout << "LCIO -> Monte Carle data, but no_mc_particle flag. MCParticles not written. \n";
         }

         // Check the scoring planes.
         for(int i=0; i< scoring_planes.size(); ++i){
            if( !has_collection(scoring_planes[i].c_str())){
               scoring_planes_active[i] = false;
               if(md_Debug & kDebug_Warning) cout << "Scoring plane " << scoring_planes[i] << " not found. Turned off.\n";
            }
         }
      }else{
         is_MC_data = false;
         use_mc_particles = false;
         use_mc_scoring = false;
         use_ecal_hits_truth = false;
      }

      if(has_collection("RFHits")){
         has_rf_hits = true;
      }else{
         has_rf_hits = false;
         if(!is_MC_data){
            if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have RF Hits. Turning of RFHit reading.\n";
         }
      }

      if (is_2016_data && is_2019_data) cout << "WOA - a file that is both 2016 and 2019 data!!!\n";
      data_type_is_known = true;

      /////////////////////////////////////////////////////////////////////////////////////////////////
      ///
      ///  Check the LCIO Data file content.
      ///  Safety switches. We check the collections names to make sure the needed data is in the file.
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////////

      if (use_ecal_hits && !has_collection("EcalCalHits")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have EcalCalHits. Turning of ECal hit reading. \n";
         use_ecal_hits = false;
      }

      if (use_ecal_hits && use_ecal_hits_truth && !has_collection("EcalHits")) {
         if( (md_Debug & kDebug_Warning) ) cout << "WARNING: The LCIO file does not have EcalHits. Ecal Hits Truth turned off.\n";
         use_ecal_hits_truth = false;
      }


      if (use_ecal_cluster && !has_collection("EcalClustersCorr")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have EcalClustersCorr. Turning of ECal corrected cluster reading. \n";
         use_ecal_cluster = false;
      }

      if (use_ecal_cluster_uncor && !has_collection("EcalClusters")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have EcalClusters. Turning of ECal cluster reading. \n";
         use_ecal_cluster_uncor = false;
      }

      if (use_hodo_hits && !has_collection("HodoCalHits")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have HodoCalHits. Turning of Hodoscope hit reading. \n";
         use_hodo_hits = false;
      }

      if (use_hodo_clusters && !has_collection("HodoGenericClusters")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have HodoGenericClusters. Turning of Hodoscope cluster reading. \n";
         use_hodo_clusters = false;
      }

      if (use_svt_raw_hits && (!has_collection("SVTRawTrackerHits") ||
                               !has_collection("SVTShapeFitParameters") ||
                               !has_collection("SVTFittedRawTrackerHits"))) {
         if(md_Debug & kDebug_Warning) {
            cout << "WARNING: The LCIO file does not have SVTRawTrackerHits or " <<
                 "SVTShapeFitParameters or SVTFittedRawTrackerHits.\n";
            cout << "         Turning of SVT raw hit writing. \n";
         }
         use_svt_raw_hits = false;
      }

      if (use_svt_hits && !(has_collection("RotatedHelicalTrackHits") ||
                            has_collection("StripClusterer_SiTrackerHitStrip1D"))) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have RotatedHelicalTrackHits. Turning of svt 3D hit writing. \n";
         use_svt_hits = false;
      } else {
         if (has_collection("StripClusterer_SiTrackerHitStrip1D"))
            svt_hit_collections.push_back("StripClusterer_SiTrackerHitStrip1D");
         if (has_collection("RotatedHelicalTrackHits"))
            svt_hit_collections.push_back("RotatedHelicalTrackHits");
      }
      if ((use_kf_tracks || use_kf_particles) && !has_collection("KalmanFullTracks")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have KalmanFullTracks. Turning of KF track writing. \n";
         use_kf_tracks = false;
         use_kf_particles = false;
      }
      if ((use_gbl_tracks || use_gbl_particles) && !has_collection("GBLTracks")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have GBLTracks. Turning of GBL track writing. \n";
         use_gbl_tracks = false;
         use_gbl_particles = false;
      }
      if (use_matched_tracks && !has_collection("MatchedTracks")) {
         if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have MatchedTracks. Turning of matched track writing. \n";
         use_matched_tracks = false;
      }
      for (auto type = particle_types_single.begin(); type< particle_types_single.end(); ++type) {
         string collection_name = Type_to_Collection[*type];
         if(md_Debug & kDebug_L1){
            cout << "Checking for " << collection_name << " in lcio file.\n";
         }
         if(!has_collection(collection_name.c_str())){
            if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have " << collection_name <<". Removing from list.\n";
            // See https://www.techiedelight.com/remove-elements-vector-inside-loop-cpp
            // We need to make sure the iterator is not invalidated by the erase.
            particle_types_single.erase(type--);
         }
      }
      for (auto type = particle_types_double.begin(); type< particle_types_double.end(); ++type) {
         string collection_name = Type_to_Collection[*type];
         if(md_Debug & kDebug_L1){
            cout << "Checking for " << collection_name << " in lcio file.\n";
         }
         if(!has_collection(collection_name.c_str())){
            if(md_Debug & kDebug_Warning) cout << "WARNING: The LCIO file does not have " << collection_name <<". Removing from list.\n";
            // See https://www.techiedelight.com/remove-elements-vector-inside-loop-cpp
            // We need to make sure the iterator is not invalidated by the erase.
            particle_types_double.erase(type--);
         }
      }
      data_type_is_known = true;
   }
}

bool LcioReader::Process(Long64_t entry){
   /// Process all the information from the current LCIO event.
   Clear();  // Clear all the vectors that contain data, so we can push_back on them again.

   /////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// Event Header
   ///
   /////////////////////////////////////////////////////////////////////////////////////////////

   run_number = lcio_event->getRunNumber();
   event_number = lcio_event->getEventNumber();
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
                (svt_burstmode_noise_good << 1) +
                (svt_latency_good << 2) +
                (svt_position_good << 3) +
                (svt_event_header_good << 4);

   // Get the LCIO GenericObject collection containing the RF times
   if(!is_MC_data && has_rf_hits) {

      EVENT::LCCollection *rf_hits
            = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("RFHits"));

      if (rf_hits->getNumberOfElements() != 1) {
         if (md_Debug & kDebug_L1)
            cout << "Issue with rf_hits for event " << event_number << " == No rf hits found. \n";
         rf_time1 = -99.;
         rf_time2 = -99.;
      } else {
         EVENT::LCGenericObject *rf_times = static_cast<EVENT::LCGenericObject *>(rf_hits->getElementAt(0));
         if (rf_times->getNDouble() != 2) {
            cout << "ERROR in LCIO - wrong number of hits. \n";
         } else {
            rf_time1 = rf_times->getDoubleVal(0);
            rf_time2 = rf_times->getDoubleVal(1);
         }
      }
   }else{
      rf_time1 = 0;
      rf_time2 = 0;
   }
   //
   // For trigger bit parsing, see EvioTool::TSBank.h
   // Note that LCIO has the type of bank in location 0, so compared to the EVIO data
   // all the words are shifted +1.
   //
   if (is_2016_data) {
      EVENT::LCCollection *lcio_triggers
            = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("TriggerBank"));
      int n_trigger_banks = lcio_triggers->getNumberOfElements();
      static int n_warning{0};
      if (n_trigger_banks != 3 && n_warning < 3) {
         // There should be 3 banks from EVIO: 0xe10c (57612) 0xe10a (57610) and 0xe10f (57615)
         cout << "Holy cow, there were inconsistent number of triggers in this event! Hmmm, could be MC.\n";
         n_warning++;
      };
      for (int i = 0; i < n_trigger_banks; ++i) {
         EVENT::LCGenericObject *lcio_trig
               = static_cast<EVENT::LCGenericObject *>(lcio_triggers->getElementAt(i));
         int nvalues = lcio_trig->getNInt();
         int test_val = lcio_trig->getIntVal(0);
         if (md_Debug & kDebug_L4)
            printf("Trigger int 0 = 0x%04x (%6d) has %2d int values.\n",
                   test_val, test_val, nvalues);
         if (test_val == 57610) { // Trigger bits are here.
            unsigned int trigger_bits = lcio_trig->getIntVal(1);
            // We need to re-encode the bits so they are consistent with 2019 data. See MiniDst.h
            TriggerBits_int_t trig_code{0};
            trig_code.bits.Single_0_Top = (trigger_bits & (1 << 24)); // Note: sets to 0 or 1.
            trig_code.bits.Single_0_Bot = trig_code.bits.Single_0_Top;
            trig_code.bits.Single_1_Top = (trigger_bits & (1 << 25));
            trig_code.bits.Single_1_Bot = trig_code.bits.Single_1_Top;
            trig_code.bits.Pair_0 = (trigger_bits & (1 << 26));
            trig_code.bits.Pair_1 = (trigger_bits & (1 << 27));
            trig_code.bits.Pulser = (trigger_bits & (1 << 29));
            trigger = trig_code.intval; // Finally, set the trigger to be stored.

            time_stamp = lcio_trig->getIntVal(4) & 0xFFFF; // Time stamp high word.
            time_stamp <<= 32;
            time_stamp += lcio_trig->getIntVal(3); // Time stamp low word.
         }
      } // for each trigger_bank
   } // is_2016_data
   if (is_2019_data && !is_MC_data) {
      EVENT::LCCollection *tsbank_data
            = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("TSBank"));
      EVENT::LCGenericObject *lcio_ts_bank
            = static_cast<EVENT::LCGenericObject *>(tsbank_data->getElementAt(0));


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
   /// Hodoscope
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////

   if (use_hodo_raw_hits){
      auto ecal_raw_hits = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("HodoscopeReadoutHits"));

      for(int ihit = 0; ihit < ecal_raw_hits->getNumberOfElements(); ++ihit){
         auto raw_hit = static_cast<EVENT::TrackerRawData *>(ecal_raw_hits->getElementAt(ihit));

         Long64_t value = Long64_t(raw_hit->getCellID0() & 0xffffffff) |
                          (Long64_t(raw_hit->getCellID1()) << 32);
         hodo_hit_field_decoder.setValue(value);
         hodo_raw_ix.push_back(hodo_hit_field_decoder["ix"]);
         hodo_raw_iy.push_back(hodo_hit_field_decoder["iy"]);
         hodo_raw_hole.push_back(hodo_hit_field_decoder["hole"]);
         hodo_raw_layer.push_back(hodo_hit_field_decoder["layer"]);
         vector<short> raw_hit_adc_short = raw_hit->getADCValues();
         hodo_raw_adc.push_back(raw_hit_adc_short);
      }
   }

   /// Parse the "HodoCalHits"
   if (use_hodo_hits){
      EVENT::LCCollection *hodo_hits =
            static_cast<EVENT::LCCollection *>(lcio_event->getCollection("HodoCalHits"));
      for(int ihit = 0; ihit < hodo_hits->getNumberOfElements(); ++ihit){
         IMPL::CalorimeterHitImpl *lcio_hit =
               static_cast<IMPL::CalorimeterHitImpl *>(hodo_hits->getElementAt(ihit));
         // hodo_hit_to_index_map[lcio_hit] = ihit;
         hodo_hit_energy.push_back(lcio_hit->getEnergy());
         hodo_hit_time.push_back(lcio_hit->getTime());

         Long64_t value = Long64_t(lcio_hit->getCellID0() & 0xffffffff) |
                          (Long64_t(lcio_hit->getCellID1()) << 32);
         hodo_hit_field_decoder.setValue(value);

         hodo_hit_index_x.push_back(hodo_hit_field_decoder["ix"]);
         hodo_hit_index_y.push_back(hodo_hit_field_decoder["iy"]);
         hodo_hit_hole.push_back(hodo_hit_field_decoder["hole"]);
         hodo_hit_layer.push_back(hodo_hit_field_decoder["layer"]);
      }
   }

   /// Parse the Hodo Genergic Object: HodoGenericClusters
   if (use_hodo_clusters){
      EVENT::LCCollection *generic =
            static_cast<EVENT::LCCollection *>(lcio_event->getCollection("HodoGenericClusters"));
      auto gclus_ix = static_cast<IMPL::LCGenericObjectImpl *>(generic->getElementAt(0));
      auto gclus_iy = static_cast<IMPL::LCGenericObjectImpl *>(generic->getElementAt(1));
      auto gclus_layer = static_cast<IMPL::LCGenericObjectImpl *>(generic->getElementAt(2));
      auto gclus_energy = static_cast<IMPL::LCGenericObjectImpl *>(generic->getElementAt(3));
      auto gclus_time = static_cast<IMPL::LCGenericObjectImpl *>(generic->getElementAt(4));

      for(int i=0; i< gclus_ix->getNInt(); ++i){  // Only one loop needed, since we know all have same size.
         hodo_cluster_ix.push_back(gclus_ix->getIntVal(i));
         hodo_cluster_iy.push_back(gclus_iy->getIntVal(i));
         hodo_cluster_layer.push_back(gclus_layer->getIntVal(i));
         hodo_cluster_energy.push_back(gclus_energy->getDoubleVal(i));
         hodo_cluster_time.push_back(gclus_time->getDoubleVal(i));
      }
   }

   ////////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// Ecal
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// Ecal Hits are stored in: EcalHits, EcalCalHits, EcalUncalHits, EcalReadoutHits
   /// In each case, a corresponding hit will have an identical cellid0 and cellid1.
   /// The "EcalTruthRelations" connects the EcalReadoutHits to the EcalHits.
   ///
   /// EcalReadoutHits contain the raw ADC information from data, for MC this is derived from the EcalHits
   /// EcalHits        contain the MC true energy and position, with a list of PDG particle ids and energies.
   /// EcalUncalHits   contain the calculated uncalibrated hits from EcalReadoutHits.
   /// EcalCalHits     conatian the calculated calibrated hits from EcalUncalHits.
   ///
   /// The order of the hits appears is not consisted throughout but the cellids track, so truth relations seem
   /// to be not needed and instead a lookup can be used.
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////

   if (use_ecal_raw_hits){
      auto ecal_raw_hits = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalReadoutHits"));

      for(int ihit = 0; ihit < ecal_raw_hits->getNumberOfElements(); ++ihit){
         auto raw_hit = static_cast<EVENT::TrackerRawData *>(ecal_raw_hits->getElementAt(ihit));

         Long64_t value = Long64_t(raw_hit->getCellID0() & 0xffffffff) |
                          (Long64_t(raw_hit->getCellID1()) << 32);
         ecal_hit_field_decoder.setValue(value);
         ecal_raw_ix.push_back(ecal_hit_field_decoder["ix"]);
         ecal_raw_iy.push_back(ecal_hit_field_decoder["iy"]);
         vector<short> raw_hit_adc_short = raw_hit->getADCValues();
         ecal_raw_adc.push_back(raw_hit_adc_short);
      }
   }

   /// Parse the "EcalCalHits"
   if (use_ecal_hits) {
      auto ecal_hits = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalCalHits"));

//            EVENT::LCCollection *ecal_thruth_relation{NULL};
//            unique_ptr<UTIL::LCRelationNavigator> ecal_truth_nav;
//            if(use_mc_particles){  // Monte Carlo data has truth information.
//               ecal_thruth_relation = lcio_event->getCollection("EcalTruthRelations");
//               ecal_truth_nav = make_unique<UTIL::LCRelationNavigator>(ecal_thruth_relation);
//            }
//
      for (int ihit = 0; ihit < ecal_hits->getNumberOfElements(); ++ihit) {
         IMPL::CalorimeterHitImpl *lcio_hit
               = static_cast<IMPL::CalorimeterHitImpl *>(ecal_hits->getElementAt(ihit));

         ecal_hit_to_index_map[lcio_hit] = ihit;
         int id0 = lcio_hit->getCellID0();
         ecal_id0_to_hit_index[id0] = ihit;
         // Gets the CellID which identifies the crystal.
         // int id0 = lcio_hit->getCellID0();
         // 0.1 ns resolution is sufficient to distinguish any 2 hits on the same crystal.
         // int id1 = static_cast<int>(10.0*lcio_hit->getTime());

         if (use_ecal_hits) {  // Only store the results if we actually want them.
            ecal_hit_energy.push_back(lcio_hit->getEnergy());
            ecal_hit_time.push_back(lcio_hit->getTime());

            Long64_t value = Long64_t(lcio_hit->getCellID0() & 0xffffffff) |
                             (Long64_t(lcio_hit->getCellID1()) << 32);
            ecal_hit_field_decoder.setValue(value);

            ecal_hit_index_x.push_back(ecal_hit_field_decoder["ix"]);
            ecal_hit_index_y.push_back(ecal_hit_field_decoder["iy"]);
            const float *pos = lcio_hit->getPosition();
            if( pos != nullptr) {
               ecal_hit_x.push_back(pos[0]);
               ecal_hit_y.push_back(pos[1]);
               ecal_hit_z.push_back(pos[2]);
            }
         }
      }
   }

   /// Parse "EcalClustersCorr" -- corrected Ecal clusters.
   if (use_ecal_cluster) {
      // Run this also if we only store particles, since the particles need some of this info.
      EVENT::LCCollection *clusters =
            static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalClustersCorr"));
      for (int i_clus = 0; i_clus < clusters->getNumberOfElements(); ++i_clus) {
         auto lcio_clus = dynamic_cast<IMPL::ClusterImpl *>(clusters->getElementAt(i_clus));
         ecal_cluster_to_index_map[lcio_clus] = i_clus;
         ecal_cluster_energy.push_back(lcio_clus->getEnergy());
         const float *position = lcio_clus->getPosition();
         ecal_cluster_x.push_back(position[0]);
         ecal_cluster_y.push_back(position[1]);
         ecal_cluster_z.push_back(position[2]);

         // Link the hits to the cluster.
         EVENT::CalorimeterHitVec clus_hits = lcio_clus->getCalorimeterHits();

         ecal_cluster_nhits.push_back(clus_hits.size());

         double seed_energy{-99.};
         double seed_time{-99};
         int seed_index{-99};
         Long64_t seed_cellid0{-99};
         int seed_ix{-99};
         int seed_iy{-99};
         vector<int> clus_hit_indexes;
         for (int j_hit = 0; j_hit < clus_hits.size(); ++j_hit) {
            IMPL::CalorimeterHitImpl *hit = static_cast<IMPL::CalorimeterHitImpl *>(clus_hits[j_hit]);
            if (hit->getEnergy() > seed_energy) {
               seed_energy = hit->getEnergy();
               seed_time = hit->getTime();
               if (use_ecal_hits) {
                  seed_index = ecal_hit_to_index_map[hit];
               }
               seed_cellid0 = Long64_t(hit->getCellID0() & 0xffffffff) |
                              (Long64_t(hit->getCellID1()) << 32);

            }
            int hit_index = -99;
            if (use_ecal_hits)hit_index = ecal_hit_to_index_map[hit];
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

   /// Parse "EcalClusters" -- uncorrected Ecal clusters.
   if (use_ecal_cluster_uncor) {
      // Run this also if we only store particles, since the particles need some of this info.
      EVENT::LCCollection *clusters =
            static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalClusters"));
      for (int i_clus = 0; i_clus < clusters->getNumberOfElements(); ++i_clus) {
         auto lcio_clus = dynamic_cast<IMPL::ClusterImpl *>(clusters->getElementAt(i_clus));
         ecal_cluster_uncor_to_index_map[lcio_clus] = i_clus;
         ecal_cluster_uncor_energy.push_back(lcio_clus->getEnergy());
         const float *position = lcio_clus->getPosition();
         ecal_cluster_uncor_x.push_back(position[0]);
         ecal_cluster_uncor_y.push_back(position[1]);
         ecal_cluster_uncor_z.push_back(position[2]);

         // Link the hits to the cluster.
         EVENT::CalorimeterHitVec clus_hits = lcio_clus->getCalorimeterHits();

         ecal_cluster_uncor_nhits.push_back(clus_hits.size());

         double seed_energy{-99.};
         double seed_time{-99};
         int seed_index{-99};
         Long64_t seed_cellid0{-99};
         int seed_ix{-99};
         int seed_iy{-99};
         vector<int> clus_hit_indexes;
         for (int j_hit = 0; j_hit < clus_hits.size(); ++j_hit) {
            IMPL::CalorimeterHitImpl *hit = static_cast<IMPL::CalorimeterHitImpl *>(clus_hits[j_hit]);
            if (hit->getEnergy() > seed_energy) {
               seed_energy = hit->getEnergy();
               seed_time = hit->getTime();
               if (use_ecal_hits) {
                  seed_index = ecal_hit_to_index_map[hit];
               }
               seed_cellid0 = Long64_t(hit->getCellID0() & 0xffffffff) |
                              (Long64_t(hit->getCellID1()) << 32);

            }
            int hit_index = -99;
            if (use_ecal_hits)hit_index = ecal_hit_to_index_map[hit];
            clus_hit_indexes.push_back(hit_index);
         }
         ecal_cluster_uncor_time.push_back(seed_time);        // The cluster time = seed hit time.
         ecal_cluster_uncor_hits.push_back(clus_hit_indexes);
         ecal_cluster_uncor_seed_index.push_back(seed_index);
         ecal_cluster_uncor_seed_energy.push_back(seed_energy);
         ecal_hit_field_decoder.setValue(seed_cellid0);
         ecal_cluster_uncor_seed_ix.push_back(ecal_hit_field_decoder["ix"]);
         ecal_cluster_uncor_seed_iy.push_back(ecal_hit_field_decoder["iy"]);
      }
   }

   ////////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// SVT
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////

   /// SVT Raw hit and fit parameters.
   // Hits are in SVTRawTrackerHits, CellID + 6 ADC values.
   // Fits are in SVTShapeFitParameters, 5 doubles: t0, t0_err, Amp, Amp_err, Chi2
   // These two are related through the SVTFittedRawTrackerHits LCIORelations.
   // We navigate to the SVTShapeFitParameters using the LCIORelations, so we do not access
   // the collection directly.

   if(use_svt_raw_hits){
      EVENT::LCCollection *raw_svt_hits = lcio_event->getCollection("SVTRawTrackerHits");
      EVENT::LCCollection *raw_svt_hit_rel = lcio_event->getCollection("SVTFittedRawTrackerHits");
      unique_ptr<UTIL::LCRelationNavigator> raw_hit_nav =
            make_unique<UTIL::LCRelationNavigator>(raw_svt_hit_rel);
      for(int i_rhit=0; i_rhit < raw_svt_hits->getNumberOfElements(); ++i_rhit){
         auto raw_hit = static_cast<EVENT::TrackerRawData *>(raw_svt_hits->getElementAt(i_rhit));
         ULong64_t value = (ULong64_t(raw_hit->getCellID0()) & 0xffffffff) |
                           (ULong64_t(raw_hit->getCellID1()) << 32);
         raw_svt_hit_decoder.setValue(value);

         EVENT::LCObjectVec raw_hit_fit_result_list = raw_hit_nav->getRelatedToObjects(raw_hit);
         vector<short> raw_hit_adc = raw_hit->getADCValues();

         if(raw_hit_fit_result_list.size() < 1){
            cout << "Error retrieving the fits for a raw hit.\n";
         }else{
            pair<int, int> store_indexes{-1,-2};
            for(int i_raw_fit=0; i_raw_fit < raw_hit_fit_result_list.size(); ++i_raw_fit) {
               auto raw_hit_fit = static_cast<IMPL::LCGenericObjectImpl *>(raw_hit_fit_result_list.at(i_raw_fit));
               svt_raw_hit_adc.push_back(raw_hit_adc);
               svt_raw_hit_fit_no.push_back(i_raw_fit);
               svt_raw_hit_layer.push_back(raw_svt_hit_decoder["layer"]);
               svt_raw_hit_module.push_back(raw_svt_hit_decoder["module"]);
               svt_raw_hit_strip.push_back(raw_svt_hit_decoder["strip"]);
               svt_raw_hit_t0.push_back(raw_hit_fit->getDoubleVal(0));
               svt_raw_hit_t0_err.push_back(raw_hit_fit->getDoubleVal(1));
               svt_raw_hit_amp.push_back(raw_hit_fit->getDoubleVal(2));
               svt_raw_hit_amp_err.push_back(raw_hit_fit->getDoubleVal(3));
               svt_raw_hit_chi2.push_back(raw_hit_fit->getDoubleVal(4));
               if(i_raw_fit == 0) store_indexes.first = svt_raw_hit_fit_no.size()-1;
               if(i_raw_fit == 1) store_indexes.second = svt_raw_hit_fit_no.size()-1;
            }
            svt_raw_hit_to_index_map[raw_hit] = store_indexes;
         }
      }
   }

   /// Parse the "RotatedHelicalTrackHits"       - These are the hits used by the GBL tracker.
   /// and "StripClusterer_SiTrackerHitStrip1D"  - These are the hits used by the KF tracker.
   if (use_svt_hits) {
      int i_svt_hit_type = -1;
      for( auto collection_name: svt_hit_collections) {
         i_svt_hit_type++;
         EVENT::LCCollection *tracker_hits = lcio_event->getCollection(collection_name); //
         for (int i_svt_hit = 0; i_svt_hit < tracker_hits->getNumberOfElements(); ++i_svt_hit) {
            auto lcio_svt_hit =
                  dynamic_cast<IMPL::TrackerHitImpl *>(tracker_hits->getElementAt(i_svt_hit));
            svt_hit_to_index_map[lcio_svt_hit] = i_svt_hit;

            svt_hit_type.push_back(i_svt_hit_type);
            svt_hit_time.push_back(lcio_svt_hit->getTime());
            const double *pos = lcio_svt_hit->getPosition();
            svt_hit_x.push_back(pos[0]);
            svt_hit_y.push_back(pos[1]);
            svt_hit_z.push_back(pos[2]);
            svt_hit_edep.push_back(lcio_svt_hit->getEDep()); // Not in 2016 data.
            const vector<float> cov_mat = lcio_svt_hit->getCovMatrix();
            svt_hit_cxx.push_back(cov_mat[0]);
            svt_hit_cxy.push_back(cov_mat[1]);
            svt_hit_cyy.push_back(cov_mat[2]);
            svt_hit_cxz.push_back(cov_mat[3]);
            svt_hit_cyz.push_back(cov_mat[4]);
            svt_hit_czz.push_back(cov_mat[5]);

            //ULong64_t value2 = (ULong64_t(lcio_svt_hit->getCellID0()) & 0xffffffff) |
            //                   (ULong64_t(lcio_svt_hit->getCellID1()) << 32);

            EVENT::LCObjectVec raw_hits = lcio_svt_hit->getRawHits();
            vector<int> raw_index;
            vector<int> raw_other;
            int layer;
            int module;
            vector<int> strip;
            for (int i_hit = 0; i_hit < raw_hits.size(); ++i_hit) {
               auto lcio_raw_hit = static_cast<EVENT::TrackerRawData *>(raw_hits.at(i_hit));
               if(use_svt_raw_hits){
                  auto hit_index_ptr = svt_raw_hit_to_index_map.find(lcio_raw_hit);
                  if( hit_index_ptr != svt_raw_hit_to_index_map.end() ) {
                     /// We try to disambiguate the two possible fits. The only handle we have is time.
                     int i_use_raw_index = hit_index_ptr->second.first; // default to the first.
                     int i_use_raw_other = hit_index_ptr->second.second; // default to the first.
                     if(hit_index_ptr->second.second >= 0){   // There is a second.
                        if( abs(svt_hit_time.back() - svt_raw_hit_t0[hit_index_ptr->second.first])
                            > abs(svt_hit_time.back() - svt_raw_hit_t0[hit_index_ptr->second.second])) {
                           i_use_raw_index = hit_index_ptr->second.second;
                           i_use_raw_other = hit_index_ptr->second.first;
                        }
                     }
                     raw_index.push_back(i_use_raw_index);
                     raw_other.push_back(i_use_raw_other);
                  }else{
                     raw_index.push_back(-1);
                     raw_other.push_back(-1);
                  }
               }
               ULong64_t value = (ULong64_t(lcio_raw_hit->getCellID0()) & 0xffffffff) |
                                 (ULong64_t(lcio_raw_hit->getCellID1()) << 32);
               raw_svt_hit_decoder.setValue(value);

               layer = raw_svt_hit_decoder["layer"];
               module = raw_svt_hit_decoder["module"];
               strip.push_back(raw_svt_hit_decoder["strip"]);
            }
            svt_hit_raw_index.push_back(raw_index);
            svt_hit_raw_other.push_back(raw_other);
            svt_hit_layer.push_back(layer);
            svt_hit_module.push_back(module);
            svt_hit_strip.push_back(strip);
         }
      }
   }

   // GBL Track collections: GBLTracks <= (TrackDataRelations)=< TrackData
   //                                  -(GBLKinkDataRelations)-> GBLKinkData
   //
   // Kalman Track coll:     KalmanFullTracks <=(KFTrackDataRelations)=< KFTrackData
   //                                         -(KFGBLStripClusterDataRelations) -> KFGBLStripClusterData
   ////////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// TRACKS
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////

   if (use_kf_tracks || use_gbl_tracks || use_matched_tracks ||
       use_kf_particles || use_gbl_particles) {
      // Also run this for particles, so we don't need to duplicate this code when filling particle struc.
      EVENT::LCCollection *kf_tracks{nullptr};
      unique_ptr<UTIL::LCRelationNavigator> track_data_kf_nav;
      track_n_kf = 0;
      if (use_kf_tracks || use_kf_particles) {
         kf_tracks = lcio_event->getCollection("KalmanFullTracks");
         track_n_kf = kf_tracks->getNumberOfElements();
         EVENT::LCCollection *track_data_kf_rel = lcio_event->getCollection("KFTrackDataRelations");
         track_data_kf_nav =  make_unique<UTIL::LCRelationNavigator>(track_data_kf_rel);
      }

      EVENT::LCCollection *gbl_tracks{nullptr};
      unique_ptr<UTIL::LCRelationNavigator> track_data_gbl_nav;
      unique_ptr<UTIL::LCRelationNavigator> gbl_kink_data_nav;
      track_n_gbl = 0;
      if (use_gbl_tracks || use_gbl_particles) {
         gbl_tracks = lcio_event->getCollection("GBLTracks");
         track_n_gbl = gbl_tracks->getNumberOfElements();
         EVENT::LCCollection *track_data_gbl_rel = static_cast<EVENT::LCCollection *>(
               lcio_event->getCollection("TrackDataRelations"));
         track_data_gbl_nav =  make_unique<UTIL::LCRelationNavigator>(track_data_gbl_rel);

         if(use_gbl_kink_data){
            EVENT::LCCollection *gbl_kink_data_rel =
                  static_cast<EVENT::LCCollection *>(lcio_event->getCollection("GBLKinkDataRelations"));
            gbl_kink_data_nav =
                  make_unique<UTIL::LCRelationNavigator>(gbl_kink_data_rel);
         }
      }

      EVENT::LCCollection *matched_tracks{nullptr};
      track_n_matched = 0;
      if (use_matched_tracks || use_gbl_particles) {
         matched_tracks = lcio_event->getCollection("MatchedTracks");
         track_n_matched = matched_tracks->getNumberOfElements();
      }

      int n_total_tracks = track_n_kf + track_n_gbl + track_n_matched;
      EVENT::Track *lcio_track{nullptr};

      // Relations for KF track to KFTrackData:
      //EVENT::LCCollection* track_data_kf_rel = static_cast<EVENT::LCCollection*>(
      //        lcio_event->getCollection("KFTrackDataRelations"));

      // Instantiate an LCRelation navigator which will allow faster acces to TrackData objects
      //shared_ptr<UTIL::LCRelationNavigator> track_data_kf_nav =
      //        make_shared<UTIL::LCRelationNavigator>(track_data_kf_rel);
      // UTIL::LCRelationNavigator *track_data_kf_nav = new UTIL::LCRelationNavigator(track_data_kf_rel);

      for (int track_number = 0; track_number < n_total_tracks; ++track_number) {
         bool track_is_kf = false;
         bool track_is_gbl = false;
         bool track_is_matched = false;
         if (track_number < track_n_kf) {
            track_is_kf = true;
            lcio_track = static_cast<EVENT::Track *>(kf_tracks->getElementAt(track_number));
            kf_track_to_index_map[lcio_track] = track_number;
            track_gbl_ref.push_back(-99);
            track_ref.push_back(-99);
         } else if (track_number < track_n_kf + track_n_gbl) {
            track_is_gbl = true;
            int gbl_track_number = track_number - track_n_kf;
            lcio_track = static_cast<EVENT::Track *>(gbl_tracks->getElementAt(gbl_track_number));
            gbl_track_to_index_map[lcio_track] = track_number;
            track_gbl_ref.push_back(track_number); // GBL Track points to itself.
            track_ref.push_back(-99); // Pointer to seed track is resolved later.
         } else {
            track_is_matched = true;
            int matched_track_number = track_number - track_n_kf - track_n_gbl;
            lcio_track = static_cast<EVENT::Track *>(matched_tracks->getElementAt(matched_track_number));
            matched_track_to_index_map[lcio_track] = track_number;
            track_gbl_ref.push_back(-99);  // Seed track needs to be resolved later.
            track_ref.push_back(track_number); // Seed track points to itself.
         }

         any_track_to_index_map[lcio_track] = track_number;
         track_particle.push_back(-99); // Pointer to particle is resolved *way* later.

         track_d0.push_back(lcio_track->getD0());
         track_chi2.push_back(lcio_track->getChi2());
         track_omega.push_back(lcio_track->getOmega());
         track_phi0.push_back(lcio_track->getPhi());
         track_tan_lambda.push_back(lcio_track->getTanLambda());
         track_z0.push_back(lcio_track->getZ0());
         track_type.push_back(lcio_track->getType());

         std::vector<EVENT::TrackState*> TrackStateVec = lcio_track->getTrackStates();

         if(use_extra_tracks) {
            const EVENT::TrackState *track_state
                  = lcio_track->getTrackState(EVENT::TrackState::AtLastHit);
            if (track_state) {
               const float *hit_pos = track_state->getReferencePoint();
               float omega = track_state->getOmega();
               float phi = track_state->getPhi();
               float tanlambda = track_state->getTanLambda();
               float z0 = track_state->getZ0();
               float d0 = track_state->getD0();
               // const float *mom_at_lasthit = lcio_track->getReferencePoint();  // This is not correct to vanilla hps-java output.
               // cout << "o,p,t,z: " << omega << ", " << phi << ", " << tanlambda << ", " << z0 << endl;
               track_x_at_lasthit.push_back(hit_pos[0]);  // Note: these are not rotated, they were put in detector
               track_y_at_lasthit.push_back(hit_pos[1]);  //       coordinates before being written to file.
               track_z_at_lasthit.push_back(hit_pos[2]);
               //track_px_at_lasthit.push_back(mom_at_lasthit[0]);
               //track_py_at_lasthit.push_back(mom_at_lasthit[1]);
               //track_pz_at_lasthit.push_back(mom_at_lasthit[2]);
               track_omega_at_lasthit.push_back(omega);
               track_phi0_at_lasthit.push_back(phi);
               track_tan_lambda_at_lasthit.push_back(tanlambda);
               track_z0_at_lasthit.push_back(z0);
               track_d0_at_lasthit.push_back(d0);
            } else {
               track_x_at_lasthit.push_back(-9999.);
               track_y_at_lasthit.push_back(-9999.);
               track_z_at_lasthit.push_back(-9999.);
               //track_px_at_lasthit.push_back(-9999.);
               //track_py_at_lasthit.push_back(-9999.);
               //track_pz_at_lasthit.push_back(-9999.);
               track_omega_at_lasthit.push_back(-9999.);
               track_phi0_at_lasthit.push_back(-9999.);
               track_tan_lambda_at_lasthit.push_back(-9999.);
               track_z0_at_lasthit.push_back(-9999.);
               track_d0_at_lasthit.push_back(-9999.);
            }
         }

         const EVENT::TrackState *track_state
              // = lcio_track->getTrackState(EVENT::TrackState::LastLocation);
              = lcio_track->getTrackState(EVENT::TrackState::AtCalorimeter);
         if (track_state) {
            const float *ecal_pos = track_state->getReferencePoint();
            track_x_at_ecal.push_back(ecal_pos[1]);  // Note: Un-rotate from tracking coordinate system.
            track_y_at_ecal.push_back(ecal_pos[2]);
            track_z_at_ecal.push_back(ecal_pos[0]);
         } else {
            track_x_at_ecal.push_back(-9999.);
            track_y_at_ecal.push_back(-9999.);
            track_z_at_ecal.push_back(-9999.);
         }

         EVENT::LCObjectVec track_data_list;
         if(track_is_kf){
            track_data_list = track_data_kf_nav->getRelatedFromObjects(lcio_track);
         }else if(track_is_gbl) {
            track_data_list = track_data_gbl_nav->getRelatedFromObjects(lcio_track);
         }

         vector<double> iso_values(14, -99.);
         if (track_data_list.size() == 1) {
            // There should always be one and only one for GBL tracks, 0 for matched tracks.
            IMPL::LCGenericObjectImpl *track_info =
                  static_cast<IMPL::LCGenericObjectImpl *>(track_data_list.at(0));

//            track_state
//                  = lcio_track->getTrackState(EVENT::TrackState::AtIP);

            double px{-999.}, py{-999.}, pz{-999.};

            // Sanity check...
            if (track_is_gbl && (track_info->getNDouble() < 12 || track_info->getNDouble() > 14 ||  /* 2016 or 2019 */
                                 track_info->getNFloat() < 4   || track_info->getNInt() < 1) ||
                track_is_kf &&  (track_info->getNFloat() < 4 || track_info->getNInt() < 1) ){
               static int n_warning{0};

               if(n_warning < 2) {  // Only show this warning twice.
                  if(md_Debug & kDebug_Warning) {
                     std::cout << "Dude! This looks like an old SLCIO file. track_info->getNFloat() = "
                               << track_info->getNFloat() << "\n";
                     std::cout << "Using magnetic field strength of " << magnetic_field << " to calculate px,py,pz "
                               << std::endl;
                  }
                  n_warning++;
               }
               // Compute the px, py, pz component from the magnetic field and the tracking parameters. See hps-java TrackUtils.java
               double omega = lcio_track->getOmega();
               if (abs(omega) < 0.0000001) {
                  omega = 0.0000001;
               }
               double Pt = abs((1. / omega) * magnetic_field * 2.99792458E-4);
               // Get px, py, pz already taking into account the tracking coordinate to detector coordinate conversion.
               pz = Pt * cos(lcio_track->getPhi());
               px = Pt * sin(lcio_track->getPhi());
               py = Pt * lcio_track->getTanLambda();
            }

            for (int i_iso = 0; i_iso < track_info->getNDouble(); ++i_iso) {
               iso_values[i_iso] = track_info->getDoubleVal(i_iso);
            }
            track_time.push_back(track_info->getFloatVal(0));
            track_volume.push_back(track_info->getIntVal(0));

            if(track_info->getNFloat() >= 4) {
               px = track_info->getFloatVal(1);
               py = track_info->getFloatVal(2);
               pz = track_info->getFloatVal(3);
            }
            track_px.push_back(px);
            track_py.push_back(py);
            track_pz.push_back(pz);


         } else {
            if(track_is_gbl) cout << "Track without TrackData for type GBL\n";
            if(track_is_kf) cout << "Track without KFTrackData for type KF\n";
            track_time.push_back(-999.);
            track_volume.push_back(-1);
            track_px.push_back(-999.);
            track_py.push_back(-999.);
            track_pz.push_back(-999.);
         }
         track_isolation.push_back(iso_values);

         EVENT::TrackerHitVec tracker_hits = lcio_track->getTrackerHits();
         track_n_hits.push_back(tracker_hits.size());

         if(use_svt_hits ) {
            // Get the collection of 3D hits associated with a LCIO Track

            vector<int> svt_hits;
            svt_hits.reserve(tracker_hits.size());
            for (int i_trk = 0; i_trk < tracker_hits.size(); ++i_trk) {
               auto trk_hit = dynamic_cast<IMPL::TrackerHitImpl *>(tracker_hits[i_trk]);
               auto svt_hit_ptr = svt_hit_to_index_map.find(trk_hit);
               if (svt_hit_ptr == svt_hit_to_index_map.end()) {
                  cout << "Track - hit association error. Hit not found in map.\n";
                  svt_hits.push_back(-99);
               } else{
                  int hit_index = svt_hit_to_index_map[trk_hit];
                  svt_hits.push_back(hit_index);
               }
            }
            track_svt_hits.push_back(svt_hits);
         }else{
            track_n_hits.push_back(-1);
         }

         // TODO: It would be better to permit vector<vector<float>> and store that
         //  instead of converting to double.
         const vector<float> cov_matrix_f = lcio_track->getCovMatrix();
         vector<double> cov_matrix_d(cov_matrix_f.begin(), cov_matrix_f.end());
         track_covmatrix.push_back(cov_matrix_d);

         /// Store the GBL Kink information, if you care to.
         if (use_gbl_kink_data && track_is_gbl) {     // This is GBL track.
            // Get the list of GBLKinkData associated with the LCIO Track

            EVENT::LCObjectVec gbl_kink_data_list = gbl_kink_data_nav->getRelatedFromObjects(lcio_track);
            IMPL::LCGenericObjectImpl *gbl_kink_datum{nullptr};
            if (gbl_kink_data_list.size() == 1) {
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
            } else {
               track_lambda_kinks.push_back(vector<double>(0));
               track_phi_kinks.push_back(vector<double>(0));
            }
         }
      } // End loop over tracks.

      if(use_gbl_tracks) {
         // Get the collection of LCRelations between seed tracks and a GBL tracks.
         EVENT::LCCollection *seed_to_gbl_relations =
               static_cast<EVENT::LCCollection *>(lcio_event->getCollection("MatchedToGBLTrackRelations"));

         // Instantiate an LCRelation navigator which will allow faster access
         // to the seed to GBL LCRelations
         unique_ptr<UTIL::LCRelationNavigator> seed_to_gbl_relations_nav =
               make_unique<UTIL::LCRelationNavigator>(UTIL::LCRelationNavigator(seed_to_gbl_relations));

         for (auto const &[gbl_track, gbl_track_index]: gbl_track_to_index_map) {
            EVENT::LCObjectVec seed_to_gbl_list
                  = seed_to_gbl_relations_nav->getRelatedFromObjects(gbl_track);
            if (seed_to_gbl_list.size() != 1) {
               cout << "Woops, I expected only one seed track for a gbl track.\n";
            } else {
               auto *seed_track = dynamic_cast<EVENT::Track *>(seed_to_gbl_list.at(0));
               int seed_index = matched_track_to_index_map[seed_track];
               track_ref[gbl_track_index] = seed_index;
               if (track_gbl_ref.size() > seed_index) {
                  int debug_copy_seed_index = seed_index;
                  int debug_copy_gbl_track_index = gbl_track_index;
                  track_gbl_ref[seed_index] = gbl_track_index;
               }
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
   /// Possible collections are:
   /// particle_types_single{FINAL_STATE_PARTICLE, OTHER_ELECTRONS, etc};
   /// particle_types_double{TC_V0_CANDIDATE, UC_VC_CANDIDATE, etc};
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////

   if(use_kf_particles || use_gbl_particles ){
      // Single particles: These are part_xxx in the output tree.
      int i_part_all =0;
      for (int type : particle_types_single) {
         string collection_name = Type_to_Collection[type];
         EVENT::LCCollection *particles = lcio_event->getCollection(collection_name);
         for (int i_part = 0; i_part < particles->getNumberOfElements(); ++i_part) {
            auto *lcio_part = dynamic_cast<EVENT::ReconstructedParticle *>(particles->getElementAt(i_part));
            Fill_Single_Particle_From_LCIO(&part, lcio_part, type);
            any_particle_to_index_map[lcio_part] = i_part_all;
            i_part_all++;
         }
      }

      // Double Track Particles, i.e.  Constrained V0 Candidates: These are v0_xxx in the output tree.
      for (int type: particle_types_double) {
         string collection_name = Type_to_VertexCollection[type];
         EVENT::LCCollection *vertexes = lcio_event->getCollection(collection_name);
         for (int i_vertex = 0; i_vertex < vertexes->getNumberOfElements(); ++i_vertex) {
            auto *lcio_vertex = dynamic_cast<EVENT::Vertex *>(vertexes->getElementAt(i_vertex));
            Fill_Vertex_From_LCIO(&v0, lcio_vertex, type);
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// MCParticles - Monte Carlo specific information.
   ///
   ///////////////////////////////////////////////////////////////////////////////////////////////

   if(use_mc_particles) {
      EVENT::LCCollection *mc_part_col = lcio_event->getCollection("MCParticle");
      map<int, int> id_to_id;
      for (int i_part = 0; i_part < mc_part_col->getNumberOfElements(); ++i_part) {
         auto mc_part = dynamic_cast<EVENT::MCParticle *>(mc_part_col->getElementAt(i_part));
         int this_id = mc_part->id();
         id_to_id[this_id] = i_part;
         mc_part_id.push_back(this_id);
         mc_part_pdg.push_back(mc_part->getPDG());
         mc_part_energy.push_back(mc_part->getEnergy());
         mc_part_mass.push_back(mc_part->getMass());
         int simulation_status = mc_part->getGeneratorStatus() | mc_part->getSimulatorStatus();
         mc_part_sim_status.push_back(simulation_status);
         mc_part_time.push_back(mc_part->getTime());
         mc_part_charge.push_back(mc_part->getCharge());
         const double *part_vertex = mc_part->getVertex();
         mc_part_x.push_back(part_vertex[0]);
         mc_part_y.push_back(part_vertex[1]);
         mc_part_z.push_back(part_vertex[2]);
         const double *part_mom = mc_part->getMomentum();
         mc_part_px.push_back(part_mom[0]);
         mc_part_py.push_back(part_mom[1]);
         mc_part_pz.push_back(part_mom[2]);
         const double *part_end = mc_part->getEndpoint();
         mc_part_end_x.push_back(part_end[0]);
         mc_part_end_y.push_back(part_end[1]);
         mc_part_end_z.push_back(part_end[2]);
//               const double *part_end_mom = mc_part->getMomentumAtEndpoint();
//               mc_part_end_px.push_back(part_end_mom[0]);
//               mc_part_end_py.push_back(part_end_mom[1]);
//               mc_part_end_pz.push_back(part_end_mom[2]);
      }

      // Now we loop again to resolve the parent and daughter ids correctly.
      for (int i_part = 0; i_part < mc_part_col->getNumberOfElements(); ++i_part) {
         auto mc_part = dynamic_cast<EVENT::MCParticle *>(mc_part_col->getElementAt(i_part));
         EVENT::MCParticleVec parent_particles = mc_part->getParents();
         vector<int> parents;
         for (auto parent: parent_particles) {
            int parent_id = parent->id();
            auto idid = id_to_id.find(parent_id);
            if (idid != id_to_id.end()) {
               int parent_id_id = idid->second;
               parents.push_back(parent_id_id);
            } else {
               cout << "MCParticle: unidentified parent.\n";
               parents.push_back(-1);
            }
         }
         mc_part_parents.push_back(parents);
         EVENT::MCParticleVec daughter_particles = mc_part->getDaughters();
         vector<int> daughters;
         for (auto daughter: daughter_particles) {
            int daughter_id = daughter->id();
            auto idid = id_to_id.find(daughter_id);
            if (idid != id_to_id.end()) {
               int daughter_id_id = idid->second;
               daughters.push_back(daughter_id_id);
            } else {
               cout << "MCParticle: unidentified daughter.\n";
               daughters.push_back(-1);
            }
         }
         mc_part_daughters.push_back(daughters);
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////
      /// MC Scoring planes.
      ////////////////////////////////////////////////////////////////////////////////////////////////
      if (use_mc_scoring) {
         for (int type = 0; type < scoring_planes.size(); ++type) {
            if(scoring_planes_active[type]) {
               EVENT::LCCollection *mc_simtrackerhit_col = lcio_event->getCollection(
                     scoring_planes[type].c_str());
               for (int i = 0; i < mc_simtrackerhit_col->getNumberOfElements(); ++i) {
                  auto mc_score = dynamic_cast<EVENT::SimTrackerHit *>(mc_simtrackerhit_col->getElementAt(
                        i));
                  mc_score_type.push_back(type);
                  auto mc_particle = mc_score->getMCParticle();
                  int mc_part_id = mc_particle->id();
                  auto idid = id_to_id.find(mc_part_id);
                  if (idid != id_to_id.end()) {
                     mc_score_part_idx.push_back(idid->second);
                  } else {
                     mc_score_part_idx.push_back(-1);
                  }
                  const float *score_part_mom = mc_score->getMomentum();
                  mc_score_px.push_back(score_part_mom[0]);
                  mc_score_py.push_back(score_part_mom[1]);
                  mc_score_pz.push_back(score_part_mom[2]);
                  const double *score_hit_pos = mc_score->getPosition();
                  mc_score_x.push_back(score_hit_pos[0]);
                  mc_score_y.push_back(score_hit_pos[1]);
                  mc_score_z.push_back(score_hit_pos[2]);
                  mc_score_time.push_back(mc_score->getTime());
                  mc_score_pdg.push_back(mc_score->getMCParticle()->getPDG());
               }
            }
         }
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///
      /// ADD ECal Truth for MC data.
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////

      vector<int> ecal_truth_to_hit_index;  // A table for each truth hit pointing to the ecal hit.
      if(use_ecal_hits && use_ecal_hits_truth) {
         // We want to add the "truth" information to hits. However, NOT EACH HIT HAS TRUTH.
         // No idea why, but just look in the LCIO file, EcalCalHits is often larger than EcalHits.
         // A cursory check shows some low energy hits from EcalCalHits do not have a corresponding item in EcalHits.

         auto ecal_hits = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalCalHits"));
         auto ecal_truth = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("EcalHits"));

         // Go through all the TRUTH hits from the EcalHits collection.
         for(int i_truth = 0; i_truth < ecal_truth->getNumberOfElements(); ++i_truth){
            IMPL::SimCalorimeterHitImpl *lcio_ecal_truth  =
                  static_cast<IMPL::SimCalorimeterHitImpl *>(ecal_truth->getElementAt(i_truth));

            // Find the corresponding hit in the EcalCalHits collection, and check is they are the same.
            int hit_idx = ecal_id0_to_hit_index[lcio_ecal_truth->getCellID0()];
            IMPL::CalorimeterHitImpl *lcio_hit
                  = static_cast<IMPL::CalorimeterHitImpl *>(ecal_hits->getElementAt(hit_idx));
            if( lcio_hit->getCellID0() != lcio_ecal_truth->getCellID0() ) {
               int id_hit = lcio_hit->getCellID0();
               int id_truth = lcio_ecal_truth->getCellID0();
               printf("======+++== The TRUTH did not match the ECAL HIT for idx = %2d,%2d %08d != %08d \n",
                      i_truth, hit_idx, id_hit, id_truth);
               ecal_truth_to_hit_index.push_back(-1);
               continue;
            }
            ecal_truth_to_hit_index.push_back(hit_idx);  // Store the truth -> ecal hit id.
            // We get the list of MC particles related to this hit. Can be more than one, can be quite a lot.
            int nmcc = lcio_ecal_truth->getNMCContributions();
            vector<int> mc_part_index_list;  // Index of the MC Particle
            vector<int> mc_part_pdg_list;    // PDG
            vector<double> mc_part_ec_list;  // Energy contribution
            int ultimate_parent_idx = -1;
            int ultimate_parent_pdg = -9999;
            double ultimate_parent_energy_contribution = -1.;
            map<int, double> map_id_to_ec_sum;

            // For each of the contributing MC Particles, collect the energy contributed to the hit,
            // and the MC particle that deposited the energy.
            for(int i_mcp=0; i_mcp<nmcc; ++i_mcp){
               // int truth_pdg = lcio_ecal_truth->getPDGCont(i_mcp);     // This is useless, always = 0.
               double truth_e = lcio_ecal_truth->getEnergyCont(i_mcp); //
               mc_part_ec_list.push_back(truth_e);
               auto mc_particle = lcio_ecal_truth->getParticleCont(i_mcp);
               int mc_part_id = mc_particle->id();
               auto idid = id_to_id.find(mc_part_id);
               if (idid != id_to_id.end()) {
                  mc_part_index_list.push_back(idid->second);
                  mc_part_pdg_list.push_back(mc_part_pdg[idid->second]);
               } else {
                  mc_part_index_list.push_back(-1);
                  mc_part_pdg_list.push_back(-9999);  // invalid PDG.
               }

               // We now want to find the ultimate parent MC particle.
               // Under normal conditions, if there are multiple hits in this crystal, the ultimate parent
               // will be the same particle for these hits.
               // If this is not the case, we take the hit that contributed most
               // energy in the crystal, and take the parent of that.
               auto parent_particle = mc_particle;
               int n_parents = parent_particle->getParents().size();

               // We walk up the parent tree until there is no more parent. This is the ultimate parent.
               while( n_parents > 0){
                  //parent_particle = mc_particle->getParents()[0];  // In our case, only one parent per MCParticle.
                  auto parents = parent_particle->getParents();
                  parent_particle = parents[0];
                  n_parents = parent_particle->getParents().size();
                  if(n_parents > 1){
                     printf("More than one parent particle????? \n");
                  }
               }

               // We add the contributions of each MC Particle to the energy of the parent in a map.
               idid = id_to_id.find(parent_particle->id());
               auto idec = map_id_to_ec_sum.find(idid->second);
               if (idec != map_id_to_ec_sum.end()){ // Already had one, so add it.
                  map_id_to_ec_sum[idid->second] += truth_e;
               }else{ // Not found, add it.
                  map_id_to_ec_sum[idid->second] = truth_e;
               }
            }

            // Determine the ultimate parent by choosing that one that contributed most to the hit.
            // Usually this is trivial.
            // Find the highest contribution, then use the info relating to that one.
            double max_energy = 0;
            for(auto &itt: map_id_to_ec_sum){
               if(max_energy < itt.second){
                  max_energy = itt.second;
                  ultimate_parent_idx = itt.first;
                  ultimate_parent_pdg = mc_part_pdg[ultimate_parent_idx];
               }
            }

            ecal_hit_mc_contrib_id.push_back(mc_part_index_list);
            ecal_hit_mc_contrib_pdg.push_back(mc_part_pdg_list);
            ecal_hit_mc_contrib_ec.push_back(mc_part_ec_list);
            ecal_hit_mc_parent_id.push_back(ultimate_parent_idx);
            ecal_hit_mc_parent_pdg.push_back(ultimate_parent_pdg);
         }
      }

      if(use_ecal_cluster){
         // Sort through the cluster hits to determine the best guess parentage of the cluster.
         for(int ic=0; ic< ecal_cluster_hits.size(); ++ic){
            double n_tot=0.;
            map<int,double> pdg_count;  // Assumes auto initialization to zero of new elements
            for(int ih=0; ih< ecal_cluster_hits[ic].size(); ++ih){
               int hit_id = ecal_cluster_hits[ic][ih];
               auto found = std::find(ecal_truth_to_hit_index.begin(),ecal_truth_to_hit_index.end(),hit_id);
               if(found == ecal_truth_to_hit_index.end()){
                  // This hit does not have any truth, so we just ignore it.
               }else{
                  int truth_id = std::distance(ecal_truth_to_hit_index.begin(),found);
                  int p_id = ecal_hit_mc_parent_id[truth_id];
                  double weight = ecal_hit_energy[truth_id];
                  pdg_count[p_id] += weight;
                  n_tot += weight;
               }
            }
            // Find the maximum item in the pdg_count map.
            if(pdg_count.size()) {
               auto mymax = std::max_element(pdg_count.begin(), pdg_count.end(),
                                             [](const std::pair<int, double> &a,
                                                const std::pair<int, double> &b) -> bool {
                                                return a.second < b.second;
                                             });
               int parent_id = mymax->first;
               ecal_cluster_mc_id.push_back(parent_id);
               ecal_cluster_mc_pdg.push_back(mc_part_pdg[parent_id]);
               ecal_cluster_mc_pdg_purity.push_back(mymax->second / n_tot);
            }else{
               // There was NO truth hit for this cluster. So fill out invalids.
               ecal_cluster_mc_id.push_back(-1);
               ecal_cluster_mc_pdg.push_back(-999);
               ecal_cluster_mc_pdg_purity.push_back(0);
            }
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// TODO: ADD SVT Truth relation for MC data.
   ///
   ///////////////////////////////////////////////////////////////////////////////////////////////

   return true;
}

long LcioReader::Run(int max_event) {


   for (const string &file: input_files) {

      if (md_Debug & kDebug_L1) cout << "Opening file : " << file << endl;


      is_2016_data = false;
      is_2019_data = false;
      data_type_is_known = false;

      lcio_reader->open(file);

      while ((lcio_event = lcio_reader->readNextEvent())) {

         if(!data_type_is_known) SetupLcioDataType();

         if (md_Debug & kDebug_Info) {
            if ((++evt_count) % Counter_Freq == 0) {
               printf("i: %'10lu   event: %'10d  run: %5d\n", evt_count, event_number, run_number);
            }
         }
         if (max_event > 0 && evt_count > max_event) break;  // End the loop, we are done.

         Process(evt_count);
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

void LcioReader::Fill_Vertex_From_LCIO(Vertex_Particle_t *vp, EVENT::Vertex *lcio_vert, int type) {
   /// Fill a vertex from the LCIO EVENT::Vertex object.
   /// Note that the Parameters has a different (ad hoc) meaning for 2016 (older hps-java) and 2019 (later hps-java) LCIO files.
   const float *vertex_pos = lcio_vert->getPosition();
   vp->type.push_back(type);
   vp->vertex_x.push_back(vertex_pos[0]);
   vp->vertex_y.push_back(vertex_pos[1]);
   vp->vertex_z.push_back(vertex_pos[2]);
   vp->vertex_chi2.push_back(lcio_vert->getChi2());
   vp->vertex_prob.push_back(lcio_vert->getProbability());
   // You can check that the algorithm does indeed correspond to the "type"
   // const std::string &xxx= lcio_vert->getAlgorithmType(); // TargetConstrained/Unconstrained/...
   // cout << "Vertex of type: " << xxx << " is type: " << type << "\n";
   std::vector<float> params = lcio_vert->getParameters();
   if (params.size() == 8 ){  // The older algorithm, mostly for 2016 data.
      // invM,p1X, p2Y, p2X, p1Z, p2Z, p1Y,invMerr
      vp->mass.push_back(params[0]);
      vp->mass_err.push_back(params[7]);
      TVector3 p1(params[1],params[6],params[4]);
      TVector3 p2(params[3],params[2],params[5]);
      TVector3 pv = p1+p2;
      vp->px.push_back(pv.Px());
      vp->py.push_back(pv.Py());
      vp->pz.push_back(pv.Pz());

   }else{  // Newer algorithm, 2019 data.
      //0 V0PzErr, 1 invMass, 2 V0Pz, 3 vXErr, 4 V0Py, 5 V0Px, 6 V0PErr, 7 V0TargProjY, 8 vZErr,
      //9 V0TargProjXErr, 10 vYErr, 11 V0TargProjYErr, 12 invMassError, 13 p1X, 14 p2Y, 15 p2X,
      // 16 V0P, 17 p1Z, 18 p1Y, 19 p2Z, 20 V0TargProjX, 21 layerCode, 22 V0PxErr, 23 V0PyErr
      vp->mass.push_back(params[1]);
      vp->mass_err.push_back(params[12]);
//        TVector3 p1(params[13],params[18],params[17]);
//        TVector3 p2(params[15],params[14],params[19]);
//        TVector3 pv = p1+p2;
      vp->px.push_back(params[5]);
      vp->py.push_back(params[4]);
      vp->pz.push_back(params[2]);
   }

   //
   // ToDo: Call Fill_Basic_Particle_From_LCIO instead, but then don't set px,py,pz in Fill_Vertex_From_LCIO!
   //
   EVENT::ReconstructedParticle *vertex_part = lcio_vert->getAssociatedParticle();
   Fill_Basic_Particle_From_LCIO(vp, vertex_part, false);

   double check_vx_energy = vertex_part->getEnergy();
//    v0.energy.push_back(check_vx_energy);
   EVENT::ParticleID *lcio_part_id = vertex_part->getParticleIDUsed();
   if(lcio_part_id){
      int pdg = lcio_part_id->getPDG();
      //       v0.pdg.push_back(pdg);
   }else{
      //     v0.pdg.push_back(0);
   }
   // v0.charge.push_back(vertex_part->getCharge());
   // v0.goodness_of_pid.push_back(vertex_part->getGoodnessOfPID());

   double check_vx_mass = vertex_part->getMass();
   const double *check_vx_mom = vertex_part->getMomentum();
   double check_vx_px = check_vx_mom[0];
   double check_vx_py = check_vx_mom[1];
   double check_vx_pz = check_vx_mom[2];

   //const vector<EVENT::Track *> &tracks = vertex_part->getTracks();
   //const EVENT::ClusterVec &cluster_vec = vertex_part->getClusters();
   const EVENT::ReconstructedParticleVec &daughters = vertex_part->getParticles();
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
   Fill_SubPart_From_LCIO(&v0.em, electron, type);
   Fill_SubPart_From_LCIO(&v0.ep, positron, type);
}

void LcioReader::Fill_SubPart_From_LCIO(Sub_Particle_t *sub,EVENT::ReconstructedParticle *daughter, int type){
   // Fill the Sub_Particle_t structure from the LCIO daughter particle.

   // Find the index to the daughter particle in the store.
   auto i_daughter_ptr = any_particle_to_index_map.find(daughter);
   int i_part = -99;
   if(i_daughter_ptr == any_particle_to_index_map.end()) { // We did not find the daughter particle.
      if( md_Debug & kDebug_Info ) {
         cout << "We did not find the daughter particle of type " << type << " for a vertex in Run " <<
              run_number << "::" << event_number << ". Adding it.\n";
      }
      // Add the particle to the particles list. We mark it by giving it the type of the vertex collection.
      Fill_Single_Particle_From_LCIO(&part, daughter, type);
      i_part = part.type.size();
      any_particle_to_index_map[daughter] = i_part;

   }else{
      i_part = i_daughter_ptr->second;
   }

   sub->part.push_back(i_part);
   sub->good_pid.push_back(daughter->getGoodnessOfPID());
   const double *momentum = daughter->getMomentum();
   double daughter_p = sqrt( momentum[0]*momentum[0]+ momentum[1]*momentum[1]+momentum[2]*momentum[2]);
   sub->p.push_back(daughter_p);

// Find the particle assoc with this sub particle from FINAL_STATE_PARTICLES.
//    EVENT::LCCollection *particles;
//    if( type <= 8 ){  // Kalman Filter track type.
//        particles = lcio_event->getCollection(
//                Type_to_Collection[FINAL_STATE_PARTICLE_KF].c_str());
//    }else {  // GBL track type
//        particles = lcio_event->getCollection(
//                Type_to_Collection[FINAL_STATE_PARTICLE_GBL].c_str());
//    }
//    bool found_daughter_particle = false;
//    for(int i_part=0; i_part < particles->getNumberOfElements(); ++i_part) {
//        EVENT::ReconstructedParticle *lcio_part =
//                static_cast<EVENT::ReconstructedParticle *>(particles->getElementAt(i_part));
//        if(lcio_part == daughter){
//            sub->part.push_back(i_part);
//            sub->good_pid.push_back(daughter->getGoodnessOfPID());
//            const double *momentum = daughter->getMomentum();
//            double daughter_p = sqrt( momentum[0]*momentum[0]+ momentum[1]*momentum[1]+momentum[2]*momentum[2]);
//            sub->p.push_back(daughter_p);
//
//            found_daughter_particle = true;
//            break;
//        }
//    }

   const EVENT::TrackVec &tracks = daughter->getTracks();       // Get the tracks associated with this daughter.
   if(tracks.size() !=1 ){
      /// A sub particle should *always* have a track.
      cout << "Vertex sub particle with " << tracks.size() << " tracks! Expected 1. -- Run " <<
           run_number << "::" << event_number << "\n";
      md_abort_tree_fill = true;
      return;
   }

   // Here we will look to match the track to the tracks that were stored.
   //
   // ToDo: Clean this up, instead of running over the tracks, use the previously found track_ptr to index maps.
   // ToDo: Since we check both GBL and Matched tracks, make sure we can distinguish the two.
   //
   // The tracks associated with a vertex should always come from the GBLTracks collection.
   // At August 2020, for the 2019 data set, this is not the case, so we also search the MatchedTracks collection.
   // With the new tracking, it may come from the KF tracking set too.
   EVENT::Track *track = tracks[0]; // There should be only one track.
   auto i_stored_track_ptr = any_track_to_index_map.find(track);

   if( i_stored_track_ptr == any_track_to_index_map.end() ) {
      cout << "Track was not found in the any_track_to_index_map. \n";
      sub->track_time.push_back(-9999);
      sub->track_nhit.push_back(-99);
      sub->pos_ecal_x.push_back(-9999);
      sub->pos_ecal_y.push_back(-9999);
      sub->track.push_back(-99);
   }else {
      int i_stored_track = i_stored_track_ptr->second;
      if (abs(track->getChi2() - track_chi2[i_stored_track]) > 1.0E-6) {
         cout << "Not sure we got the right track: chi2_1 = " << track->getChi2() <<
              "  track_chi2_2[" << i_stored_track <<"] = " << track_chi2[i_stored_track] << "\n";
      }
      sub->track_time.push_back(track_time[i_stored_track]);
      sub->track_nhit.push_back(track_n_hits[i_stored_track]);
      sub->pos_ecal_x.push_back(track_x_at_ecal[i_stored_track]);
      sub->pos_ecal_y.push_back(track_y_at_ecal[i_stored_track]);
      sub->track.push_back(i_stored_track);
   }

   // Chi2.
   sub->chi2.push_back(track->getChi2());

   const EVENT::ClusterVec &clusters = daughter->getClusters();
   if(clusters.size() > 1 ){
      // A sub particle can have 1 (normal) or 0 (ECal hole) clusters.
      cout << "Vertex sub particle with " << clusters.size() << " clusters! Expected 1 or 0 -- Run "<<
           run_number << "::" << event_number << "\n";
      md_abort_tree_fill = true;
      return;
   }

   if(clusters.size() == 1) {
      UTIL::BitField64 ecal_hit_field_decoder("system:6,layer:2,ix:-8,iy:-6");

      IMPL::ClusterImpl *clus = static_cast<IMPL::ClusterImpl *>(clusters[0]);
      auto i_clus_ptr = ecal_cluster_to_index_map.find(clus);
      if( i_clus_ptr == ecal_cluster_to_index_map.end()){
         cout << "ERROR == Cluster pointer from LCIO without associated index. CHECK THIS!! \n";
         md_abort_tree_fill=true;
         return;
      }

      int iclus = i_clus_ptr->second;
      sub->clus.push_back(iclus);
      sub->clus_energy.push_back(clus->getEnergy());

      const float *position = clus->getPosition();
      sub->clus_pos_x.push_back(position[0]);
      sub->clus_pos_y.push_back(position[1]);

      // Link the hits to the cluster to find the seed. (hit with most energy.
      double seed_energy{-99.};
      int seed_ix{-99};
      int seed_iy{-99};
      double seed_time{-99};

      EVENT::CalorimeterHitVec clus_hits = clus->getCalorimeterHits();
      for (int j_hit = 0; j_hit < clus_hits.size(); ++j_hit) {
         IMPL::CalorimeterHitImpl *hit = static_cast<IMPL::CalorimeterHitImpl *>(clus_hits[j_hit]);
         if (hit->getEnergy() > seed_energy) {
            seed_energy = hit->getEnergy();
            seed_time = hit->getTime();
            Long64_t value = Long64_t(hit->getCellID0() & 0xffffffff) |
                             (Long64_t(hit->getCellID1()) << 32);
            ecal_hit_field_decoder.setValue(value);
            seed_ix = ecal_hit_field_decoder["ix"];
            seed_iy = ecal_hit_field_decoder["iy"];
         }
      }
      sub->clus_time.push_back(seed_time);
      sub->clus_ix.push_back(seed_ix);
      sub->clus_iy.push_back(seed_iy);
   }else{ // Clusterless track. Fill cluster stuff with -99.
      sub->clus.push_back(-99);
      sub->clus_energy.push_back(-99.);
      sub->clus_time.push_back(-99.);
      sub->clus_pos_x.push_back(-99.);
      sub->clus_pos_y.push_back(-99.);
      sub->clus_ix.push_back(-99.);
      sub->clus_iy.push_back(-99.);
   }
}

void LcioReader::Fill_Single_Particle_From_LCIO(Single_Particle_t *bp, EVENT::ReconstructedParticle *lcio_part,
                                                int type) {
   Fill_Basic_Particle_From_LCIO(dynamic_cast<Basic_Particle_t *>(bp),lcio_part);
   bp->type.push_back(type);

   const vector<EVENT::Track *> &tracks = lcio_part->getTracks();
#ifdef DEBUG
   if(tracks.size()>1) {
        cout << "Single Particle with more than one associated track. Hmmm, check the LCIO!! \n";
    }
#endif
   if( tracks.size() == 1) {
      EVENT::Track *track = tracks[0];
      auto i_stored_track_ptr = any_track_to_index_map.find(track);

      if( i_stored_track_ptr == any_track_to_index_map.end() ) {
         if(gDebug>2) cout << "Track was not found in the any_track_to_index_map. \n";
         bp->track.push_back(-99);
      }else {
         bp->track.push_back(i_stored_track_ptr->second);
         track_particle[i_stored_track_ptr->second] = bp->track.size()-1;
      }
      bp->track_chi2.push_back(track->getChi2());
   }else{ // No track for this particle, then it is most likely a neutral.
      bp->track.push_back(-1);
      bp->track_chi2.push_back(-99.);
   }

   const vector<EVENT::Cluster *> &clusters = lcio_part->getClusters();
#ifdef DEBUG
   if(clusters.size()>1){
        cout << "Single particle with more than one cluster. Yikes. Check the LCIO!! \n";
    }
#endif
   if( clusters.size() == 1){
      auto cluster = dynamic_cast<IMPL::ClusterImpl*>(clusters[0]);
      if(use_ecal_cluster && ecal_cluster_to_index_map.find(cluster) == ecal_cluster_to_index_map.end() ){
         cout << "Um, sorry, but I could not find the cluster associated with the particle, though I expected one. \n";
         bp->ecal_cluster.push_back(-101);
      }else {
         int cluster_index = ecal_cluster_to_index_map[cluster];
         bp->ecal_cluster.push_back(cluster_index);
      }
   }else{  // Ecal hole track.
      bp->ecal_cluster.push_back(-1);
   }
}

void LcioReader::Fill_Basic_Particle_From_LCIO(Basic_Particle_t *bp, EVENT::ReconstructedParticle *lcio_part,
                                               bool fill_momentum) {
   // Fill the Basic_Particle_t structure with the information from a ReconstructedParticle.
   bp->lcio_type.push_back(lcio_part->getType());
   bp->energy.push_back(lcio_part->getEnergy());
   EVENT::ParticleID *lcio_part_id = lcio_part->getParticleIDUsed();
   if(lcio_part_id){
      int pdg = lcio_part_id->getPDG();
      bp->pdg.push_back(pdg);
   }else{
      bp->pdg.push_back(0);
   }

   bp->charge.push_back(lcio_part->getCharge());
   bp->goodness_of_pid.push_back(lcio_part->getGoodnessOfPID());
   if(fill_momentum) {
      // For Vertex particles, it is "better" to fill these directly from the vertex parameters.
      // Better? Not sure, it seems the information is simply duplicated between the vertex and the vertex-particle.
      bp->mass.push_back(lcio_part->getMass());
      const double *momentum = lcio_part->getMomentum();
      bp->px.push_back(momentum[0]);
      bp->py.push_back(momentum[1]);
      bp->pz.push_back(momentum[2]);
   }
}

void LcioReader::End(){
   MiniDst::End();
};

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
        is_2016_data = false;
        lcio_reader->open(file);

        while( (lcio_event = lcio_reader->readNextEvent())){

            if( !data_type_is_known){ // Determine the data type by looking at the collections
                col_names = lcio_event->getCollectionNames();
                if(md_Debug & kDebug_L1){
                    cout << "LCIO Collections in " << file << ":\n";
                    for(string s: *col_names){
                        cout << s << endl;
                    }
                }
                for(string s: *col_names){
                    if( s == "TriggerBank"){
                        is_2016_data = true;
                        if(md_Debug & kInfo) cout << "LCIO -> This is 2015/2016 data. \n";
                    }
                    if( s == "TSBank"){
                        is_2019_data = true;
                        if(md_Debug & kInfo) cout << "LCIO -> This is 2019 data. \n";
                    }
                    if( s == "MCParticles"){
                        is_MC_data = true;
                        if(md_Debug & kInfo) cout << "LCIO -> This is Monte Carlo data. \n";
                    }
                }
                if( is_2016_data && is_2019_data) cout << "WOA - a file that is both 2016 and 2019 data!!!\n";
                data_type_is_known = true;
            }


            run_number = lcio_event->getRunNumber();
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

                        unsigned int time_stamp_hi =
                        time_stamp = lcio_trig->getIntVal(3); // Time stamp high word.
                        time_stamp <<= 32;
                        time_stamp += lcio_trig->getIntVal(2); // Time stamp low word.
                    }
                } // for each trigger_bank
            } // is_2016_data
            if(is_2019_data){
                EVENT::LCCollection *tsbank_data
                        = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("TSBank"));
                EVENT::LCCollection *vtpbank_data
                        = static_cast<EVENT::LCCollection *>(lcio_event->getCollection("VTPBank"));

            }

            if(md_output_tree) md_output_tree->Fill();
        }

        lcio_reader->close();
    }
    return(0);
};

void LcioReader::End(){
    MiniDst::End();
};

//
// Created by Maurik Holtrop on 8/30/24.
//

#include "HpstrReader.h"

HpstrReader::HpstrReader(const string &input_file, const int debug_level) {
   if(md_Debug > 0) md_Debug = debug_level;
   cout << "LcioReader Debug level is " << md_Debug << std::endl;
   if(!input_file.empty()) input_files.push_back(input_file);
   HpstrReader(input_files, debug_level);
};

HpstrReader::HpstrReader(const vector<string> &infile_list, const int debug_level){
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

   fChain = new TChain("HPS_Event");
   for(auto file: input_files){
      fChain->Add(file.c_str());
   }
   Init();
};

HpstrReader::HpstrReader(TChain *chain, int debug_level) {
   md_Debug = debug_level;
   if(md_Debug > 0) cout << "LcioReader Debug level is " << md_Debug << std::endl;
   fChain = chain;
   Init();
}

void HpstrReader::Init(){
   /// Initialze the branches on the TTree.
   if(fChain == nullptr){
      cout << "No TChain. Cannot initialize. \n";
      return;
   }
   fChain->SetBranchAddress("EventHeader", &event_header, &b_event_header);
   fChain->SetBranchAddress("HodoHit", &hodo_hit, &b_hodo_hit);
   fChain->SetBranchAddress("HodoCluster", &hodo_cluster, &b_hodo_cluster);


}

void HpstrReader::Start(){
   if( md_Debug>0 ) {
      printf("Hpstr READER version " HpstrReader__Version__ "\n");
   }

   MiniDst::Start();
}

void HpstrReader::Clear(){
   /// Clear event storage.
   MiniDst::Clear();

}


bool HpstrReader::Process(Long64_t entry) {
   /// Process all the information from the current LCIO event.
   Clear();  // Clear all the vectors that contain data, so we can push_back on them again.

   /////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// Event Header
   ///
   /////////////////////////////////////////////////////////////////////////////////////////////

   run_number = event_header->getRunNumber();
   event_number = event_header->getEventNumber();
   time_stamp = event_header->getEventTime();

   TriggerBits_int_t trigbits{0};
   trigbits.bits.Pair_0 = event_header->isPair0Trigger();
   trigbits.bits.Pair_1 = event_header->isPair1Trigger();
   trigbits.bits.Single_0_Bot = event_header->isSingle0Trigger();
   trigbits.bits.Single_0_Top = event_header->isSingle0Trigger();
   trigbits.bits.Single_1_Bot = event_header->isSingle1Trigger();
   trigbits.bits.Single_1_Top = event_header->isSingle1Trigger();

   svt_status = event_header->isSvtBiasOn() +
                ( (!event_header->hasSvtBurstModeNoise()) << 1) +
                ( event_header->isSvtLatencyGood() << 2) +
                ( event_header->isSvtClosed() << 3) +
                ( (!event_header->hasSvtEventHeaderErrors()) << 4);

   trigger = trigbits.intval;
   ext_trigger = 0;
   rf_time1 = event_header->getRfTime(0);
   rf_time2 = event_header->getRfTime(1);

   ////////////////////////////////////////////////////////////////////////////////////////////////
   ///
   /// Hodoscope
   ///
   ////////////////////////////////////////////////////////////////////////////////////////////////



   return true;

}
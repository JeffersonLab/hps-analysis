//
// Created by Maurik Holtrop on 1/15/24.
//
#include "Hbooker.h"

COMMON_BLOCK_DEF(PAWC_DEF,PAWC);
extern PAWC_DEF PAWC;
COMMON_BLOCK_DEF(DATAC_DEF,DATAC);
extern DATAC_DEF DATAC;


Hbooker::Hbooker(): Hbooker(nullptr, nullptr, "minidst-ntuple.hbook"){
};

Hbooker::Hbooker(MiniDst *mini_dst, TChain *in_tree, string fileout):
         mdst(mini_dst), tree(in_tree), hb_output_file_name(fileout) {
// Initialize class
// Initialize HBOOK system
   cout << "Call HLIMIT()\n";
   HLIMIT(PAWC_SIZE);
   cout << "Class initialized\n";
}

void Hbooker::Start() {

   if( hb_Debug & MiniDst::kDebug_Info  ) {
      printf("Hbooker version " __Hbooker__Version__ "\n");
   }

   HROPEN(1,(char *)"EX",hb_output_file_name.data(),(char *)"N",record_size,istat);

   char dir[]="//EX";
   HCDIR(dir,(char *)" ");
   HBNT(hbook_id,(char *)"NT",(char *)"");

   // Add the items to the NTuple.
   HBNAME(hbook_id,(char *)"DAT",DATAC.run_number,(char *)"RNUM:I*4");
   HBNAME(hbook_id,(char *)"DAT",DATAC.event_number,(char *)"EVTNUM:I*4");
   HBNAME(hbook_id,(char *)"DAT",DATAC.trigger,(char *)"TRIG:I*4");
   HBNAME(hbook_id,(char *)"DAT",DATAC.rf_time1,(char *)"RF1:R*8");
   HBNAME(hbook_id,(char *)"DAT",DATAC.rf_time2,(char *)"RF2:R*8");

};


//   HBNAME(hbook_id,(char *)"DAT",DATAC.x,(char *)"X(10)");
//   HBNAME(hbook_id,(char *)"DAT",DATAC.y,(char *)"Y(10)");
//   HBNAME(hbook_id,(char *)"DAT",DATAC.z,(char *)"Z(10)");

//   for (int i=0;i<1000000;i++){
//      DATAC.x[0] = 100.*rand()/RAND_MAX;
//      DATAC.y[1] = 200.*rand()/RAND_MAX;
//      DATAC.z[2] = 300.*rand()/RAND_MAX;
//      HFNT(10);
//   }

void Hbooker::Run(long num_evt) {
   // Run through the events in the mdst and process each.
   if(num_evt == 0 || num_evt > tree->GetEntries())
      num_evt = tree->GetEntries();
   for(long evt=0; evt< num_evt; ++evt){
      tree->GetEntry(evt);
      Process();
   }
}

void Hbooker::Process() {
   // Process one event. It is assumed the event is already "loaded" into mdst.
   // Unfortunately, we need to **copy** all the data into the common block.
   DATAC.run_number = mdst->run_number;
   DATAC.event_number = mdst->event_number;
   DATAC.trigger = mdst->trigger;
   DATAC.rf_time1 = mdst->rf_time1;
   DATAC.rf_time2 = mdst->rf_time2;
   HFNT(hbook_id);
}


void Hbooker::End() {
   // Finalize the writing of the file.
   int icycle = 0;
   HROUT(0,icycle,(char *)" ");
   HREND((char *)"EX");
}

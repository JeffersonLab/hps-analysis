//
// Created by Maurik Holtrop on 1/15/24.
//
#include "Hbooker.h"

COMMON_BLOCK_DEF(PAWC_DEF,PAWC);
PAWC_DEF PAWC;
COMMON_BLOCK_DEF(DATAC_DEF,DATAC);
DATAC_DEF DATAC;


Hbooker::Hbooker(): Hbooker("minidst_ntuple.hbook"){
};

Hbooker::Hbooker(string fileout): hb_output_file_name(fileout) {
// Initialize class
// Initialize HBOOK system
   HLIMIT(PAWC_SIZE);
}

void Hbooker::Start() {

   if( md_Debug>0 ) {
      printf("Hbooker version " __Hbooker__Version__ "\n");
   }
   MiniDst::Start();

   HROPEN(1,(char *)"EX",hb_output_file_name.data(),(char *)"N",record_size,istat);

   char dir[]="//EX";
   HCDIR(dir,(char *)" ");
   HBNT(hbook_id,(char *)"NT",(char *)"");

   // Add the items to the NTuple.
   HBNAME(hbook_id,(char *)"DAT",run_number,(char *)"RNUM:I*4");
   HBNAME(hbook_id,(char *)"DAT",event_number,(char *)"EVTNUM:I*4");
   HBNAME(hbook_id,(char *)"DAT",trigger,(char *)"TRIG:I*4");
   HBNAME(hbook_id,(char *)"DAT",rf_time1,(char *)"RF1:R*8");
   HBNAME(hbook_id,(char *)"DAT",rf_time2,(char *)"RF2:R*8");

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

void Hbooker::End() {
   // Finalize the writing of the file.
   int icycle = 0;
   HROUT(0,icycle,(char *)" ");
   HREND((char *)"EX");
   MiniDst::End();
}

void Hbooker::Process() {
   // Process one event.
   HFNT(hbook_id);
}

void Hbooker::Clear() {
   /// Clear the event storage.
   // Make sure you also call the "super"
   MiniDst::Clear();
}
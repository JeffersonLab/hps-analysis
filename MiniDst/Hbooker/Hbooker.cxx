//
// Created by Maurik Holtrop on 1/15/24.
//
#include "Hbooker.h"

COMMON_BLOCK_DEF(PAWC_DEF,PAWC);
extern PAWC_DEF PAWC;
COMMON_BLOCK_DEF(EVNT_DEF, EVNT);
extern EVNT_DEF EVNT;
COMMON_BLOCK_DEF(HODO_DEF, HODO);
extern HODO_DEF HODO;
COMMON_BLOCK_DEF(ECAL_DEF, ECAL);
extern ECAL_DEF ECAL;
COMMON_BLOCK_DEF(PART_DEF, PART);
extern PART_DEF PART;

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
   HBNAME(hbook_id,(char *)"EVNT",EVNT.run_number,(char *)"RNUM:I*4");
   HBNAME(hbook_id,(char *)"EVNT",EVNT.event_number,(char *)"EVTNUM:I*4");
   HBNAME(hbook_id,(char *)"EVNT",EVNT.trigger,(char *)"TRIG:I*4");
   HBNAME(hbook_id,(char *)"EVNT",EVNT.rf_time1,(char *)"RF1:R*8");
   HBNAME(hbook_id,(char *)"EVNT",EVNT.rf_time2,(char *)"RF2:R*8");

   if(mdst->use_hodo_clusters){
      HBNAME(hbook_id,(char *)"HODO",HODO.n_hodo_cluster,(char *)"HODO_I[0," QQ(MAX_HODO_CLUSTER) "]:I*4");
      HBNAME(hbook_id,(char *)"HODO",HODO.hodo_cluster_energy,(char *)"HODO_E(HODO_I):R*4");
      HBNAME(hbook_id,(char *)"HODO",HODO.hodo_cluster_time,(char *)"HODO_T(HODO_I):R*4");
      HBNAME(hbook_id,(char *)"HODO",HODO.hodo_cluster_ix,(char *)"HODO_IX(HODO_I):I*4");
      HBNAME(hbook_id,(char *)"HODO",HODO.hodo_cluster_iy,(char *)"HODO_IY(HODO_I):I*4");
      HBNAME(hbook_id,(char *)"HODO",HODO.hodo_cluster_layer,(char *)"HODO_L(HODO_I):I*4");
   }

   if(mdst->use_ecal_cluster){
      HBNAME(hbook_id,(char *)"ECAL",ECAL.n_ecal_cluster,(char *)"ECAL_I[0," QQ(MAX_ECAL_CLUSTER) "]:I*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_energy,(char *)"ECAL_E(ECAL_I):R*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_time,(char *)"ECAL_T(ECAL_I):R*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_x,(char *)"ECAL_X(ECAL_I):R*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_y,(char *)"ECAL_Y(ECAL_I):R*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_z,(char *)"ECAL_Z(ECAL_I):R*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_seed_index,(char *)"ECAL_SI(ECAL_I):I*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_seed_ix,(char *)"ECAL_SIX(ECAL_I):I*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_seed_iy,(char *)"ECAL_SIY(ECAL_I):I*4");
      HBNAME(hbook_id,(char *)"ECAL",ECAL.ecal_cluster_seed_energy,(char *)"ECAL_SE(ECAL_I):R*4");
   }

   if(mdst->use_kf_particles || mdst->use_gbl_particles){
      HBNAME(hbook_id,(char *)"PART",PART.n_particles,(char *)"P_I[0," QQ(MAX_PARTICLES) "]:I*4");
      HBNAME(hbook_id,(char *)"PART",PART.type,(char *)"P_TYPE(P_I):I*4");
      HBNAME(hbook_id,(char *)"PART",PART.lcio_type,(char *)"P_LTYPE(P_I):I*4");
      HBNAME(hbook_id,(char *)"PART",PART.energy,(char *)"P_E(P_I):R*4");
      HBNAME(hbook_id,(char *)"PART",PART.pdg,(char *)"P_PDG(P_I):I*4");
      HBNAME(hbook_id,(char *)"PART",PART.charge,(char *)"P_CHARGE(P_I)[-1,1]:I*4");
      HBNAME(hbook_id,(char *)"PART",PART.goodness_of_pid,(char *)"P_GPID(P_I):R*4");
      HBNAME(hbook_id,(char *)"PART",PART.px,(char *)"P_PX(P_I):R*4");
      HBNAME(hbook_id,(char *)"PART",PART.py,(char *)"P_PY(P_I):R*4");
      HBNAME(hbook_id,(char *)"PART",PART.pz,(char *)"P_PZ(P_I):R*4");
      HBNAME(hbook_id,(char *)"PART",PART.ecal_index,(char *)"P_EI(P_I):I*4");
      HBNAME(hbook_id,(char *)"PART",PART.track_index,(char *)"P_TI(P_I):I*4");
      HBNAME(hbook_id,(char *)"PART",PART.track_chi2,(char *)"P_TCHI(P_I):R*4");
   }

};

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
   EVNT.run_number = mdst->run_number;
   EVNT.event_number = mdst->event_number;
   EVNT.trigger = mdst->trigger;
   EVNT.rf_time1 = mdst->rf_time1;
   EVNT.rf_time2 = mdst->rf_time2;

   if(mdst->use_hodo_clusters){
      int max = mdst->hodo_cluster_energy.size();
      if(max>=MAX_HODO_CLUSTER) max = MAX_HODO_CLUSTER-1;
      HODO.n_hodo_cluster = max;
      for(int i=0; i < max; ++i){
         HODO.hodo_cluster_energy[i] = (float)mdst->hodo_cluster_energy[i];
         HODO.hodo_cluster_time[i] = (float)mdst->hodo_cluster_time[i];
         HODO.hodo_cluster_ix[i] = mdst->hodo_cluster_ix[i];
         HODO.hodo_cluster_iy[i] = mdst->hodo_cluster_ix[i];
         HODO.hodo_cluster_layer[i] = mdst->hodo_cluster_layer[i];
      }
      for(int i= mdst->hodo_cluster_energy.size(); i< MAX_HODO_CLUSTER; ++i){
         HODO.hodo_cluster_energy[i] = 0;
         HODO.hodo_cluster_time[i] = 0;
         HODO.hodo_cluster_ix[i] = 0;
         HODO.hodo_cluster_iy[i] = 0;
         HODO.hodo_cluster_layer[i] = 0;
      }

   }

   if(mdst->use_ecal_cluster){
      int max = mdst->ecal_cluster_energy.size();
      if(max>=MAX_ECAL_CLUSTER) max = MAX_ECAL_CLUSTER-1;
      ECAL.n_ecal_cluster = max;
      for(int i=0; i<max; ++i) {
         ECAL.ecal_cluster_energy[i] = (float) mdst->ecal_cluster_energy[i];
         ECAL.ecal_cluster_time[i] = (float) mdst->ecal_cluster_time[i];
         ECAL.ecal_cluster_x[i] = (float) mdst->ecal_cluster_x[i];
         ECAL.ecal_cluster_y[i] = (float) mdst->ecal_cluster_y[i];
         ECAL.ecal_cluster_z[i] = (float) mdst->ecal_cluster_z[i];
         ECAL.ecal_cluster_seed_index[i] = mdst->ecal_cluster_seed_index[i];
         ECAL.ecal_cluster_seed_ix[i] = mdst->ecal_cluster_seed_ix[i];
         ECAL.ecal_cluster_seed_iy[i] = mdst->ecal_cluster_seed_iy[i];
         ECAL.ecal_cluster_seed_energy[i] = (float) mdst->ecal_cluster_seed_energy[i];
      }
      for(int i=mdst->ecal_cluster_energy.size(); i<MAX_ECAL_CLUSTER; ++i) {
         ECAL.ecal_cluster_energy[i] = 0;
         ECAL.ecal_cluster_time[i] = 0;
         ECAL.ecal_cluster_x[i] = 0;
         ECAL.ecal_cluster_y[i] = 0;
         ECAL.ecal_cluster_z[i] = 0;
         ECAL.ecal_cluster_seed_index[i] = 0;
         ECAL.ecal_cluster_seed_ix[i] = 0;
         ECAL.ecal_cluster_seed_iy[i] = 0;
         ECAL.ecal_cluster_seed_energy[i] = 0;
      }

   }

#define CHECK_NAN(x,y)  (x =(isnan(y)?-999:(float)y ))

   if(mdst->use_kf_particles || mdst->use_gbl_particles){
      int max = mdst->part.type.size();
      if(max>=MAX_PARTICLES) max = MAX_PARTICLES-1;
      PART.n_particles = max+1;
      for(int i=0; i<max; ++i) {
         PART.type[i] = mdst->part.type[i];
         PART.lcio_type[i] = mdst->part.lcio_type[i];
         CHECK_NAN(PART.energy[i],mdst->part.energy[i]);
//         if(isnan( mdst->part.energy[i])) PART.energy[i]=-999;
//         else PART.energy[i] = (float) mdst->part.energy[i];
         PART.pdg[i] = mdst->part.pdg[i];
         PART.charge[i] = mdst->part.charge[i];
         CHECK_NAN(PART.goodness_of_pid[i], mdst->part.goodness_of_pid[i]);
         CHECK_NAN(PART.px[i],mdst->part.px[i]);
         CHECK_NAN(PART.py[i],mdst->part.py[i]);
         CHECK_NAN(PART.pz[i],mdst->part.pz[i]);
         PART.ecal_index[i] = mdst->part.ecal_cluster[i];
         PART.track_index[i] = mdst->part.track[i];
         CHECK_NAN(PART.track_chi2[i],mdst->part.track_chi2[i]);
      }
   }

   HFNT(hbook_id);
}


void Hbooker::End() {
   // Finalize the writing of the file.
   int icycle = 0;
   HROUT(0,icycle,(char *)" ");
   HREND((char *)"EX");
}

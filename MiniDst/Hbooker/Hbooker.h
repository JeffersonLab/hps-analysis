//
// Created by Maurik Holtrop on 7/24/20.
//
// Class and code to write out HBook files.
// (Yes, old stuff.)
//
#ifndef HBOOKER_H
#define HBOOKER_H
#include <iostream>
#include <variant>
#include <chrono>
#include "TObject.h"
#include "MiniDst.h"

#include <stdlib.h>

#define g77Fortran 1

#include "cfortran/cfortran.h"
#include "cfortran/packlib.h"

#define __Hbooker__Version__ "1.0.0"

using namespace std;

// This is a trick to quote items.
#define Q(x) #x
#define QQ(x) Q(x)

#define PAWC_SIZE 500000

typedef struct { float PAW[PAWC_SIZE]; } PAWC_DEF;
#define PAWC COMMON_BLOCK(PAWC,pawc)
COMMON_BLOCK_DEF(PAWC_DEF,PAWC);

typedef struct {
   int run_number;
   int event_number;
   int trigger;
   double rf_time1;
   double rf_time2;
} EVNT_DEF;
#define EVNT COMMON_BLOCK(EVNT,evnt)
COMMON_BLOCK_DEF(EVNT_DEF, EVNT);

#define MAX_HODO_CLUSTER 5
typedef struct{
   int n_hodo_cluster;
   float hodo_cluster_energy[MAX_HODO_CLUSTER];
   float hodo_cluster_time[MAX_HODO_CLUSTER];
   int hodo_cluster_ix[MAX_HODO_CLUSTER];
   int hodo_cluster_iy[MAX_HODO_CLUSTER];
   int hodo_cluster_layer[MAX_HODO_CLUSTER];
} HODO_DEF;
#define HODO COMMON_BLOCK(HODO,hodo)
COMMON_BLOCK_DEF(HODO_DEF,hodo);

#define MAX_ECAL_CLUSTER 10
typedef struct{
   int n_ecal_cluster;
   float ecal_cluster_energy[MAX_ECAL_CLUSTER];
   float ecal_cluster_time[MAX_ECAL_CLUSTER];
   float ecal_cluster_x[MAX_ECAL_CLUSTER];
   float ecal_cluster_y[MAX_ECAL_CLUSTER];
   float ecal_cluster_z[MAX_ECAL_CLUSTER];
   int ecal_cluster_seed_index[MAX_ECAL_CLUSTER];
   int ecal_cluster_seed_ix[MAX_ECAL_CLUSTER];
   int ecal_cluster_seed_iy[MAX_ECAL_CLUSTER];
   float ecal_cluster_seed_energy[MAX_ECAL_CLUSTER];
   int ecal_cluster_nhits[MAX_ECAL_CLUSTER];
}ECAL_DEF;
#define ECAL COMMON_BLOCK(ECAL,ecal)
COMMON_BLOCK_DEF(ECAL_DEF, ECAL);

#define MAX_PARTICLES 10
typedef struct{
   int n_particles;
   int type[MAX_PARTICLES];
   int lcio_type[MAX_PARTICLES];
   float energy[MAX_PARTICLES];
   int pdg[MAX_PARTICLES];
   int charge[MAX_PARTICLES];
   float goodness_of_pid[MAX_PARTICLES];
   float px[MAX_PARTICLES];
   float py[MAX_PARTICLES];
   float pz[MAX_PARTICLES];
   float p[MAX_PARTICLES];
   int ecal_index[MAX_PARTICLES];
   int track_index[MAX_PARTICLES];
   float track_chi2[MAX_PARTICLES];
   int track_nhit[MAX_PARTICLES];
   float pos_ecal_x[MAX_PARTICLES];
   float pos_ecal_y[MAX_PARTICLES];
}PART_DEF;
#define PART COMMON_BLOCK(PART,part)
COMMON_BLOCK_DEF(PART_DEF, PART);

#define MAX_VERTEXES 10
typedef struct{
   int n_vertexes;
   float vertex_x[MAX_VERTEXES];
   float vertex_y[MAX_VERTEXES];
   double vertex_z[MAX_VERTEXES];
   double vertex_chi2[MAX_VERTEXES];
   int ele_index[MAX_VERTEXES];
   int pos_index[MAX_VERTEXES];
}VERT_DEF;
#define VERT COMMON_BLOCK(VERT,vert)
COMMON_BLOCK_DEF(VERT_DEF, VERT);

class Hbooker {

   MiniDst *mdst;
   TChain *tree;
   string hb_output_file_name;
   int hb_Debug = MiniDst::kDebug_Info & MiniDst::kDebug_L1;
   int hbook_id=10;
   int record_size=8191;
   int istat;

public:
   static string _version_() { return (__Hbooker__Version__); };
   Hbooker();
   explicit Hbooker(MiniDst *mini_dst, TChain *in_tree, string fileout);
   void SetMiniDst(MiniDst *mini_dst){mdst = mini_dst; }
   void SetOutputFileName(string outfile){hb_output_file_name = outfile;}
   void Start();
   void Run(long numevt);
   bool Process(Long64_t entry);
   void End();
};

#endif

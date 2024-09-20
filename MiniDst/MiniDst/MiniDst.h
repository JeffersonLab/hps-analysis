//
// Created by Maurik Holtrop on 7/24/20.
//
// This file is the header for the miniDST class.
// The miniDST class contains all the information to make the ROOT TTree objects for analyzing HPS data.
// These trees would normally be filled from the LCIO data files that are produced in the data cooking.
//
// Note that this class is needed to write the TTrees, but it is not needed at all to use the resulting
// ROOT files. This header file does contain useful information for understanding the contents of the trees,
// and some of the definitions that can be useful in your analysis.
//
#ifndef MINIDST_MINIDST_H
#define MINIDST_MINIDST_H
#include <iostream>
#include <variant>
#include <chrono>
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"

#define MiniDst_Version_ "1.2.1"
//
// The following construction defines a "Variant", a type in C++17 and later that can contain different kinds of object.
// In our case the variant contains each of the possible types that we have in the output tree.
// Note that using a variant is not always straight forward, and often requires a "visitor". Sometimes the visitor is
// simple, sometimes more complicated. See the "Clear()" function.
// **Note** that each added type in the variant will need a corresponding type in the visitor used in Clear().
//
using Multi_Value = std::variant < int *, double *, unsigned int *, ULong64_t *, std::vector<double>* ,
      std::vector<int>*, std::vector< std::vector<short> >*,  std::vector< std::vector<int> >*, std::vector< std::vector<double> >*>;
//
// Each variable that goes into the output is stored into the Multi_Branch object.
// The Multi_Branch maps the name of the TTree branch to the variable in the code.
//
using Multi_Branch = std::map<std::string, Multi_Value>;

//
//
// Same structure/union for 2019 trigger as EVIOTool::TSBank
//
using namespace std;

using TriggerBits_t  = struct{  // Trigger structure for the 2019 data set.
   bool Single_0_Top: 1; //  0   ( 150-8191) MeV (-31,31)   Low energy cluster
   bool Single_1_Top: 1; //  1   ( 300-3000) MeV (  5,31)   e+
   bool Single_2_Top: 1; //  2   ( 300-3000) MeV (  5,31)   e+ : Position dependent energy cut
   bool Single_3_Top: 1; //  3   ( 300-3000) MeV (  5,31)   e+ : HODO L1*L2  Match with cluster
   bool Single_0_Bot: 1; //  4   ( 150-8191) MeV (-31,31)   Low energy cluster
   bool Single_1_Bot: 1; //  5   ( 300-3000) MeV (  5,31)   e+
   bool Single_2_Bot: 1; //  6   ( 300-3000) MeV (  5,31)   e+ : Position dependent energy cut
   bool Single_3_Bot: 1; //  7   ( 300-3000) MeV (  5,31)   e+ : HODO L1*L2  Match with cluster
   bool Pair_0      : 1; //  8    A-prime
   bool Pair_1      : 1; //  9    Moller
   bool Pair_2      : 1; // 10    pi0
   bool Pair_3      : 1; // 11    -
   bool LED         : 1; // 12    LED
   bool Cosmic      : 1; // 13    Cosmic
   bool Hodoscope   : 1; // 14    Hodoscope
   bool Pulser      : 1; // 15    Pulser
   bool Mult_0      : 1; // 16    Multiplicity-0 2 Cluster Trigger
   bool Mult_1      : 1; // 17    Multiplicity-1 3 Cluster trigger
   bool FEE_Top     : 1; // 18    FEE Top       ( 2600-5200)
   bool FEE_Bot     : 1; // 19    FEE Bot       ( 2600-5200)
   unsigned int NA  :12; // 20-31 Not used
};

using TriggerBits_int_t = union{
   unsigned int intval;
   TriggerBits_t  bits;
};

class MiniDst : public TObject {

public:
   static string _version_(){return(MiniDst_Version_);};
   MiniDst(): md_output_file_name("minidst.root"){};
   explicit MiniDst(string output_file_name): md_output_file_name(output_file_name){};
   ~MiniDst() override = default;
   virtual void Clear();  /// Clear all the vectors and event storage.
   virtual void Start();
   virtual void DefineBranchMap(bool use_all=false); /// Setup the branch_map and branch_map_active
   virtual void SetBranchMap();    /// Use the branch_map to set the branches in the TTree that are active.
   virtual void SetBranchAddressesOnTree(TTree *); /// Connect the variables/vectors to the TTree for *reading* the tree.
   virtual long Run(int nevt);

   virtual bool Process(Long64_t entry);

   virtual void End();
   virtual void SetOutputFileName(const string& outfile){md_output_file_name=outfile;};
   virtual void SetDebugLevel(const int level){ md_Debug = level;};
   //std::map<std::string, vector<double>* > &Get_brmap(){return branch_map;};
   // #define FULL_CLEAR(xxx)  { for( auto p: xxx){delete p;}; xxx.Clear(); };

   virtual std::vector<std::string> GetBranchNames();           /// Return a list of all the POSSIBLE branch names as a vector.
   virtual std::vector<std::string> GetActiveBranchNames();     /// Return a list of all the ACTIVE branch names as a vector.
   virtual void SetBranchActive(std::string, bool active=true); /// Manipulate whether a branch is active or not.

   void Add_Scoring_Plane(const string name){
      scoring_planes.push_back(name);
   }
   void Remove_Scoring_Plane(const string name) {
      auto itr = std::find(scoring_planes.begin(), scoring_planes.end(), name);
      if (itr != scoring_planes.end()) scoring_planes.erase(itr);
   }

   template<typename T> inline void FULL_CLEAR(T& xxx){
      for(auto p: xxx){ delete p; }
      xxx.clear();
   }

   // Multi-event utilities
   int Add_Ecal_Hit(MiniDst &event_in, int i_hit);
   int Add_Ecal_Cluster(MiniDst &event_in, int i_cluster, bool also_copy_ecal_hits=false);
   int Add_Track(MiniDst &event_in, int i_track, bool also_copy_svt_hits=false);
   int Add_Particle(MiniDst &event_in, int i_particle);

public:

   enum Debug_codes {
      kDebug_Quiet   = 0x00,
      kDebug_Error   = 0x01,
      kDebug_Warning = 0x02,
      kDebug_Info    = 0x04,
      kDebug_L1      = 0x10,  // Debug levels
      kDebug_L2      = 0x20,
      kDebug_L3      = 0x30,
      kDebug_L4      = 0x40
   };

   enum ParticleType {  // Copy of Vertex and Particle types from hps-dst.
      FINAL_STATE_PARTICLE = 0,        // Kalman tracks are default.
      UC_V0_CANDIDATE      = 1,
      BSC_V0_CANDIDATE     = 2,
      TC_V0_CANDIDATE      = 3,
      UC_MOLLER_CANDIDATE  = 4,
      BSC_MOLLER_CANDIDATE = 5,
      TC_MOLLER_CANDIDATE  = 6,
      OTHER_ELECTRONS      = 7,
      UC_VC_CANDIDATE      = 8,
      FINAL_STATE_PARTICLE_KF = 0,
      UC_V0_VERTICES_KF      = 1,
      BSC_V0_VERTICES_KF     = 2,
      TC_V0_VERTICES_KF      = 3,
      UC_MOLLER_VERTICES_KF  = 4,
      BSC_MOLLER_VERTICES_KF = 5,
      TC_MOLLER_VERTICES_KF  = 6,
      OTHER_ELECTRONS_KF      = 7,
      UC_VC_VERTICES_KF      = 8,
      FINAL_STATE_PARTICLE_GBL = 0+9,
      UC_V0_VERTICES_GBL      = 1 + 9,
      BSC_V0_VERTICES_GBL     = 2 + 9,
      TC_V0_VERTICES_GBL      = 3 + 9,
      UC_MOLLER_VERTICES_GBL  = 4 + 9,
      BSC_MOLLER_VERTICES_GBL = 5 + 9,
      TC_MOLLER_VERTICES_GBL  = 6 + 9,
      OTHER_ELECTRONS_GBL      = 7+9,
      UC_VC_VERTICES_GBL      = 8 + 9
   };

   const vector<string> ParticleTypeNames{
         "Final state particle KF",
         "UC V0 vertex KF",
         "BSC V0 vertex KF",
         "TC V0 vertex KF",
         "UC Moller vertex KF",
         "BSC Moller vertex KF",
         "TC Moller vertex KF",
         "Other electron KF",
         "UC VC vertex KF",
         "Final state particle GBL",
         "UC V0 vertex GBL",
         "BSC V0 vertex GBL",
         "TC V0 vertex GBL",
         "UC Moller vertex GBL",
         "BSC Moller vertex GBL",
         "TC Moller vertex GBL",
         "Other electron GBL",
         "UC VC vertex GBL"
   };
public:

   unsigned int md_Debug{0x07};  // Start at Info+Warning+Error level.
   string  md_output_file_name;
   unsigned long Counter_Freq{1000}; // How often to print a status line.

   TFile *md_output_file{nullptr};
   TTree *md_output_tree{nullptr};
   bool md_abort_tree_fill{false};

   /// Switches that allow turning output on/off.
   bool use_hodo_raw_hits{false};
   bool use_hodo_hits{true};
   bool use_hodo_clusters{true};
   bool use_ecal_raw_hits{false};
   bool use_ecal_cluster{true};
   bool use_ecal_cluster_uncor{false};
   bool use_ecal_hits{true};
   bool use_ecal_hits_truth{true};
   bool use_svt_raw_hits{false};

   bool use_svt_hits{false};
   bool use_kf_tracks{true};
   bool use_gbl_tracks{false};
   bool use_matched_tracks{false};
   bool use_extra_tracks{false};
   bool use_gbl_kink_data{false};
   bool use_mc_particles{false};
   bool use_mc_scoring{false};
   bool use_kf_particles{true};
   bool use_gbl_particles{false};

   /// These determine what types of particles will be written out.
   // No vertex
   vector<int> particle_types_single{FINAL_STATE_PARTICLE_KF, OTHER_ELECTRONS_KF,
                                     FINAL_STATE_PARTICLE_GBL, OTHER_ELECTRONS_GBL};

   // Yes vertex.
   vector<int> particle_types_double{UC_V0_VERTICES_KF, BSC_V0_VERTICES_KF, TC_V0_VERTICES_KF,
                                     UC_VC_VERTICES_KF,
                                     UC_MOLLER_VERTICES_KF, BSC_MOLLER_VERTICES_KF, TC_MOLLER_VERTICES_KF,
                                     UC_V0_VERTICES_GBL, BSC_V0_VERTICES_GBL, TC_V0_VERTICES_GBL,
                                     UC_VC_VERTICES_GBL,
                                     UC_MOLLER_VERTICES_GBL, BSC_MOLLER_VERTICES_GBL, TC_MOLLER_VERTICES_GBL };

   /// Map that contains the name and address of each branch in the TTree.
   Multi_Branch branch_map;
   std::map<std::string, bool> branch_map_active;

   /// Simple helper function that does two things, fill the branch map and add the name to the active list.
   /// Note that you may want non-active branches. These will exist and get cleared, but are not written out.
   void branch_map_try_emplace(std::string name, Multi_Value address, bool is_active=true) {
      branch_map.try_emplace(name, address);
      branch_map_active.try_emplace(name, is_active);
   }

   /// All the items that we could possibly write out are here so that direct access
   /// to the variables is possible throughout the code, and to each class that derives from the MiniDst.
   /// Note:
   /// - For each variable you want in the output you need 3 steps (the last one for each input type!)
   /// -- 1: Create a variable below for easy access.
   /// -- 2: Add this variable and the name in the branch map. Please make these names consistent!
   /// -- 3: Fill the variable in each of the input file parsers:
   /// ---  3a: LcioReader.cxx
   /// ---  3b: Dst2016.cxx
   ///

   // Event header information
   int run_number{0};
   int event_number{0};
   ULong64_t time_stamp{0};
   unsigned int svt_status{0};  // Only useful for some 2016 data.
   unsigned int trigger{0};     // Packed trigger bits. For 2019 - prescaled bits
   unsigned int ext_trigger{0}; // Ext trigger bits. For 2019 - un-prescaled bits, N/A for 2016
   double rf_time1{0};
   double rf_time2{0};

   // Hodo RAW Hits
   vector<int>    hodo_raw_ix;
   vector<int>    hodo_raw_iy;
   vector<int>    hodo_raw_hole;
   vector<int>    hodo_raw_layer;
   vector<vector<short> > hodo_raw_adc;

   // Hodo Hits
   vector<double> hodo_hit_energy;
   vector<double> hodo_hit_time;
   vector<int> hodo_hit_index_x;
   vector<int> hodo_hit_index_y;
   vector<int> hodo_hit_hole;
   vector<int> hodo_hit_layer;

   // Hodo Clusters
   vector<double> hodo_cluster_energy;
   vector<double> hodo_cluster_time;
   vector<int>    hodo_cluster_ix;
   vector<int>    hodo_cluster_iy;
   vector<int>    hodo_cluster_layer;

   // Ecal RAW Hits
   vector<int>    ecal_raw_ix;
   vector<int>    ecal_raw_iy;
   vector<vector<short> > ecal_raw_adc;

   // Ecal Hits
   vector<double> ecal_hit_energy;
   vector<double> ecal_hit_time;
   vector<int>    ecal_hit_index_x;
   vector<int>    ecal_hit_index_y;
   // It *should* be possible to derive these from ix and iy, but that would be detector dependent.
   vector<double> ecal_hit_x;
   vector<double> ecal_hit_y;
   vector<double> ecal_hit_z;
   // For MC events, we want to store the truth with the hits.
   vector<vector<int>> ecal_hit_mc_contrib_id;   // For each hit contributor, the MCParticle indexes.
   vector<vector<int>> ecal_hit_mc_contrib_pdg;  // For each hit contributor, the MCParticle pdg nums.
   vector<vector<double>> ecal_hit_mc_contrib_ec; // Energy contributed to hit by this particle;
   vector<int> ecal_hit_mc_parent_id;          // For the ultimate parent, the MCParticle index.
   vector<int> ecal_hit_mc_parent_pdg;         // For the ultimate parent, the MCParticle pdg.

   // Ecal Clusters:
   vector<double> ecal_cluster_energy;
   vector<double> ecal_cluster_time;    // Cluster time = seed hit time.
   vector<double> ecal_cluster_x;
   vector<double> ecal_cluster_y;
   vector<double> ecal_cluster_z;
   vector<int> ecal_cluster_seed_index;
   vector<int> ecal_cluster_seed_ix;
   vector<int> ecal_cluster_seed_iy;
   vector<double> ecal_cluster_seed_energy;
   vector< vector<int> > ecal_cluster_hits;
   vector<int> ecal_cluster_nhits; // Not strictly needed, but handy. (could use ecal_cluster_hits[i].size())
   vector<int> ecal_cluster_mc_id;
   vector<int> ecal_cluster_mc_pdg;
   vector<double> ecal_cluster_mc_pdg_purity;

   // Ecal Clusters Ucor:
   vector<double> ecal_cluster_uncor_energy;
   vector<double> ecal_cluster_uncor_time;    // Cluster time = seed hit time.
   vector<double> ecal_cluster_uncor_x;
   vector<double> ecal_cluster_uncor_y;
   vector<double> ecal_cluster_uncor_z;
   vector<int> ecal_cluster_uncor_seed_index;
   vector<int> ecal_cluster_uncor_seed_ix;
   vector<int> ecal_cluster_uncor_seed_iy;
   vector<double> ecal_cluster_uncor_seed_energy;
   vector< vector<int> > ecal_cluster_uncor_hits;
   vector<int> ecal_cluster_uncor_nhits; // Not strictly needed, but handy. (could use ecal_cluster_uncor_hits[i].size())



   // RAW SVT Hits  == Not available in the 2016 DST, but could be useful for some LCIO based studies.
   //
   vector<int> svt_raw_hit_layer;
   vector<int> svt_raw_hit_module;
   vector<int> svt_raw_hit_strip;
   vector<vector<short>> svt_raw_hit_adc;
   vector<int> svt_raw_hit_fit_no;
   vector<double> svt_raw_hit_t0;
   vector<double> svt_raw_hit_t0_err;
   vector<double> svt_raw_hit_amp;
   vector<double> svt_raw_hit_amp_err;
   vector<double> svt_raw_hit_chi2;

   // SVT Hits
   vector<int>    svt_hit_type; /// 0 -> StripClusterer_SiTrackerHitStrip1D ; 1 -> RotatedHelicalTrackHits
   vector<double> svt_hit_x;
   vector<double> svt_hit_y;
   vector<double> svt_hit_z;
   vector<double> svt_hit_cxx;
   vector<double> svt_hit_cxy;
   vector<double> svt_hit_cxz;
   vector<double> svt_hit_cyy;
   vector<double> svt_hit_cyz;
   vector<double> svt_hit_czz;
   vector<double> svt_hit_time;
   vector<double> svt_hit_edep;  // Energy deposit. Not is 2016 DST.
   vector<vector<int>> svt_hit_raw_index;
   vector<vector<int>> svt_hit_raw_other;
   vector<int> svt_hit_layer;
   vector<int> svt_hit_module;
   vector<vector<int>> svt_hit_strip;

   // SVT Tracks
   int track_n_gbl{0}; /// Number of GBL tracks. The rest will be matched tracks or Kalman tracks.
   int track_n_kf{0};  /// Number of Kalman tracks.
   int track_n_matched{0};
   vector<int>  track_n_hits; /** The number of 3D hits associated with this track. */
   vector<int>  track_volume; /** The volume to which this track belongs to. */
   vector<int>  track_type;   /** The track type from LCIO. 0=matched, 1=KF, 32+=GBL */
   vector<double> track_d0;   /** The distance of closest approach to the reference point. */
   vector<double> track_phi0;     // The azimuthal angle of the momentum at the position of closest approach to the reference point.
   vector<double> track_omega; // The track curvature. The curvature is positive (negative) if the particle has a positive (negative) charge.
   vector<double> track_tan_lambda; // The slope of the track in the SY plane where S is the arc length of the helix in the xz plane.
   vector<double> track_z0;      // The y position of the track at the distance of closest approach in the xz plane.
   vector<double> track_chi2; /** The chi^2 of the track fit. */
   vector<double> track_time; /** The time of the track.  This is currently the average time of all hits composing the track. **/
   vector<double> track_px;
   vector<double> track_py;
   vector<double> track_pz;
   vector<double> track_x_at_lasthit; /** The x position of the extrapolated track at last hit */
   vector<double> track_y_at_lasthit; /** The y position of the extrapolated track at last hit */
   vector<double> track_z_at_lasthit; /** The z position of the extrapolated track at last hit */
   //vector<double> track_px_at_lasthit; /** 3 momentum at last hit obtained from the reference point. **/
   //vector<double> track_py_at_lasthit; /** This requires a hack on hps-java **/
   //vector<double> track_pz_at_lasthit;

   vector<double> track_omega_at_lasthit;
   vector<double> track_tan_lambda_at_lasthit;
   vector<double> track_phi0_at_lasthit;
   vector<double> track_d0_at_lasthit;
   vector<double> track_z0_at_lasthit;
   vector<double> track_x_at_ecal; /** The x position of the extrapolated track at the Ecal face. */
   vector<double> track_y_at_ecal; /** The y position of the extrapolated track at the Ecal face. */
   vector<double> track_z_at_ecal; /** The z position of the extrapolated track at the Ecal face. */
   vector<vector<double> > track_isolation; /** Array used to store the isolation variables for each of the sensor layers. */
   vector<vector<double> > track_covmatrix;   /** The 1/2 Covariant Matrix. This is the lower 1/2. **/
   vector<vector<double> > track_lambda_kinks; /** the lambda kinks for each of the sensor layers.  Empty for non-GBL tracks **/
   vector<vector<double> > track_phi_kinks; /** the phi kinks for each of the sensor layers. Empty for non-GBL tracks **/
   vector<int>  track_particle; /** Reference to the reconstructed particle associated with this track. */
   vector<int>  track_gbl_ref;
   vector<int>  track_ref;
   vector<vector<int> > track_svt_hits;/** Reference to the 3D hits associated with this track. */

   struct Basic_Particle_t {
      vector<int>    type;
      vector<int>    lcio_type;
      vector<double> energy;
      vector<double> mass;
      vector<int>    pdg;
      vector<int>    charge;
      vector<double> goodness_of_pid;
      vector<double> px;
      vector<double> py;
      vector<double> pz;
      int Add(Basic_Particle_t &base,int id){
         /// Add particle base@id and return index of added particle.
         type.push_back(base.type[id]);
         lcio_type.push_back(base.lcio_type[id]);
         energy.push_back(base.energy[id]);
         mass.push_back(base.mass[id]);
         pdg.push_back(base.pdg[id]);
         charge.push_back(base.charge[id]);
         goodness_of_pid.push_back(base.goodness_of_pid[id]);
         px.push_back(base.px[id]);
         py.push_back(base.py[id]);
         pz.push_back(base.pz[id]);
         return int(base.type.size() -1);
      }
      void clear(void){
         type.clear();
         lcio_type.clear();
         energy.clear();
         mass.clear();
         pdg.clear();
         charge.clear();
         goodness_of_pid.clear();
         px.clear();
         py.clear();
         pz.clear();
      }
   };

   struct Single_Particle_t : Basic_Particle_t{
      vector<int>    ecal_cluster; // At most one ecal cluster per particle.
      vector<int>    track;       // At most one track per particle.
      // Useful extra items obtained from track:
      vector<double> track_chi2;
      int Add(Single_Particle_t &single, int id, int new_cluster_id, int new_track_id){
         /// Add particle single@id and return index of added particle.
         Basic_Particle_t::Add(single, id);
         ecal_cluster.push_back(new_cluster_id);
         track.push_back(new_track_id);
         track_chi2.push_back(single.track_chi2[id]);
         return int(track.size() -1);
      }
      void clear(){
         Basic_Particle_t::clear();
         ecal_cluster.clear();
         track.clear();
         track_chi2.clear();
      }
   };

   Single_Particle_t part;

   struct Sub_Particle_t{
      /// The SubParticle structure is really just for convenience.
      /// All the information here could be found by looking up the particle,
      /// and from that particle the cluster, and then retrieve the information needed.
      /// But, I am lazy, and want easy access to these numbers without all the lookups.
      vector<int>    part;  /// index to particle
      vector<int>    track; /// index to track
      vector<int>    track_nhit; /// Number of hits on track.
      vector<double> p;          /// Track momentum
      vector<double> chi2;      /// Track chi-squared.
      vector<double> good_pid;  /// goodness of pid from particle.
      vector<double> track_time; /// track time
      vector<double> pos_ecal_x; /// track position at ecal
      vector<double> pos_ecal_y;
      vector<int>    clus;       /// index to cluster
      vector<double> clus_energy; /// cluster energy
      vector<double> clus_time;  /// cluster time
      vector<int>    clus_ix;   /// Cluster seed ix.
      vector<int>    clus_iy;   /// Cluster seed iy.
      vector<double> clus_pos_x; /// Cluster x position.
      vector<double> clus_pos_y;

      void clear(void){
         part.clear();
         track.clear();
         track_nhit.clear();
         p.clear();
         chi2.clear();
         good_pid.clear();
         track_time.clear();
         pos_ecal_x.clear();
         pos_ecal_y.clear();
         clus.clear();
         clus_energy.clear();
         clus_time.clear();
         clus_ix.clear();
         clus_iy.clear();
         clus_pos_x.clear();
         clus_pos_y.clear();
      }
   };

   struct Vertex_Particle_t: Basic_Particle_t{
      /// Since this is for HPS data, the subs are always an e+ and an e-.
      vector<double> vertex_x;
      vector<double> vertex_y;
      vector<double> vertex_z;
      vector<double> vertex_chi2;
      vector<double> vertex_prob;

      // Mass uncertainty from fit:
      vector<double> mass_err;

      Sub_Particle_t em;  // Electron
      Sub_Particle_t ep;  // Positron

      // ToDo: Add the covariance matrix.
      // vector<double> covar;

      void clear(void){
         Basic_Particle_t::clear();
         vertex_x.clear();
         vertex_y.clear();
         vertex_z.clear();
         vertex_chi2.clear();
         vertex_prob.clear();
         mass_err.clear();
         em.clear();
         ep.clear();
      }
   };

   Vertex_Particle_t v0;

   // MCParticles
   vector<double> mc_part_energy;
   vector<int>    mc_part_id;
   vector<int>    mc_part_pdg;
   vector<int>    mc_part_sim_status; /** Generator Status **/
   vector<double> mc_part_time;      /** The global creation time. */
   vector<double> mc_part_x;      /** The X vertex. */
   vector<double> mc_part_y;      /** The Y vertex. */
   vector<double> mc_part_z;      /** The Z vertex. */
   vector<double> mc_part_px;
   vector<double> mc_part_py;
   vector<double> mc_part_pz;
   vector<double> mc_part_end_x;  /** The X end point. */
   vector<double> mc_part_end_y;  /** The Y end point. */
   vector<double> mc_part_end_z;  /** The Z end point. */
//   vector<double> mc_part_end_px;  /** The X end point. */
//   vector<double> mc_part_end_py;  /** The Y end point. */
//   vector<double> mc_part_end_pz;  /** The Z end point. */
   vector<double> mc_part_mass;
   vector<double> mc_part_charge;
   vector< vector<int> > mc_part_daughters;
   vector< vector<int> > mc_part_parents;

   // Data from the scoring planes which are in SimTrackerHit collection type.
   // Store the LCIO collection names for the scoring planes to examine.
   vector<string> scoring_planes{"TrackerHits", "TrackerHitsECal", "HodoscopeHits"};
   vector<bool> scoring_planes_active{true, true, true};
   vector<int>  mc_score_type;   // Which scoring plane.
   vector<int>  mc_score_part_idx; // Index to the MCParticle that made the hit.
   vector<double> mc_score_x;    // The x position.
   vector<double> mc_score_y;    // The y position.
   vector<double> mc_score_z;    // The z position.
   vector<double> mc_score_px;   // Momentum x
   vector<double> mc_score_py;   // Momentum y
   vector<double> mc_score_pz;   // Momentum z
   vector<double> mc_score_time; // Time of crossing
   vector<int>    mc_score_pdg;     // PDG id of particle.

};

#endif //MINIDST_MINIDST_H

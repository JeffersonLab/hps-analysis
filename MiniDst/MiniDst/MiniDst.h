//
// Created by Maurik Holtrop on 7/24/20.
//

#ifndef MINIDST_MINIDST_H
#define MINIDST_MINIDST_H
#include <iostream>
#include <variant>
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

// In theory, the Multi_Value variant could also contain "int" and "double" types, however, doing so
// will complicate the "visit" functions for operations like clear(), while really providing little benefit.
// So only vector types are in the variant.
using Multi_Value = std::variant < std::vector<double>* , std::vector<int>*,
        std::vector< std::vector<int> >*, std::vector< std::vector<double> >*>;
using Multi_Branch = std::map<std::string, Multi_Value >;

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
    bool Pair_0      : 1; //  8    A'
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
        FINAL_STATE_PARTICLE = 0,
        UC_V0_CANDIDATE      = 1,
        BSC_V0_CANDIDATE     = 2,
        TC_V0_CANDIDATE      = 3,
        UC_MOLLER_CANDIDATE  = 4,
        BSC_MOLLER_CANDIDATE = 5,
        TC_MOLLER_CANDIDATE  = 6,
        OTHER_ELECTRONS      = 7,
        UC_VC_CANDIDATE      = 8
    };

public:

    Int_t   md_Debug{0x07};  // Start at Info+Warning+Error level.
    string  md_output_file_name;
    TFile *md_output_file{nullptr};
    TTree *md_output_tree{nullptr};
    bool md_abort_tree_fill{false};

    /// Switches that allow turning output on/off.
    bool write_ecal_cluster{true};
    bool write_ecal_hits{true};
    bool write_svt_hits{true};
    bool write_tracks{true};
    bool write_only_gbl_tracks{false};
    bool write_mc_particles{false};

    /// These determine what types of particles will be written out.
    vector<int> particle_types_single{FINAL_STATE_PARTICLE, OTHER_ELECTRONS}; // No vertex
    vector<int> particle_types_double{TC_V0_CANDIDATE, UC_VC_CANDIDATE}; // Yes vertex.

    // Vector type helpers.
    Multi_Branch branch_map;

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
    unsigned int ext_trigger{0};     // Ext trigger bits. For 2019 - un-prescaled bits, N/A for 2016
    double rf_time1;
    double rf_time2;

    // Ecal Hits
    vector<double> ecal_hit_energy;
    vector<double> ecal_hit_time;
    vector<int>    ecal_hit_index_x;
    vector<int>    ecal_hit_index_y;

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

    // RAW SVT Hits  == Not available in the 2016 DST, but could be useful for some LCIO based studies.
    // ToDo: Add these later?
    //
    // svt_raw_hit_layer
    // svt_raw_hit_strip
    // svt_raw_hit_charge
    // svt_raw_hit_time
    //

    // SVT Hits
    vector<int> svt_hit_layer;
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

    // SVT Tracks
    vector<int>  track_n_hits; /** The number of 3D hits associated with this track. */
    vector<int>  track_volume; /** The volume to which this track belongs to. */
    vector<int>  track_type;   /** The track type. */
    vector<double> track_d0;   /** The distance of closest approach to the reference point. */
    vector<double> track_phi0;     // The azimuthal angle of the momentum at the position of closest approach to the reference point.
    vector<double> track_omega; // The track curvature. The curvature is positive (negative) if the particle has a positive (negative) charge.
    vector<double> track_tan_lambda; // The slope of the track in the SY plane where S is the arc length of the helix in the xz plane.
    vector<double> track_z0;      // The y position of the track at the distance of closest approach in the xz plane.
    vector<double> track_chi2; /** The chi^2 of the track fit. */
    vector<double> track_time; /** The time of the track.  This is currently the average time of all hits composing the track. **/
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
//        vector<double> corr_px;
//        vector<double> corr_py;
//        vector<double> corr_pz;
    };

    struct Single_Particle_t : Basic_Particle_t{
        vector<int>    track;       // At most one track per particle.
        vector<int>    ecal_cluster; // At most one ecal cluster per particle.
    };

    Single_Particle_t part;

    struct Sub_Particle_t{
        /// The SubParticle structure is really just for convenience.
        /// All the information here could be found by looking up the particle,
        /// and from that particle the cluster, and then retreive the information needed.
        /// But, I am lazy, and want easy access to these numbers without all the lookups.
        vector<int>    part;  /// index to particle
        vector<int>    track; /// index to track
        vector<int>    track_nhit; // Number of hits on track.
        vector<double> p;          /// Track momentum
        vector<double> chi2;      /// Track chi-squared.
        vector<double> good_pid;  /// goodness of pid from particle.
        vector<double> track_time; /// track time
        vector<double> pos_ecal_x; /// track position at ecal
        vector<double> pos_ecal_y;
        vector<int>    clus;       /// index to cluster
        vector<double> clus_energy; /// cluster energy
        vector<double> clus_time;  /// cluster time
        vector<int>    clus_ix;   // Cluster seed ix.
        vector<int>    clus_iy;   // Cluster seed iy.
        vector<double> clus_pos_x; /// Cluster x position.
        vector<double> clus_pos_y;
    };

    struct Vertex_Particle_t: Basic_Particle_t{
        vector<double> n_daughter; // Not needed, should alsways be = 2.
        /// Since this is for HPS data, the subs are always an e+ and an e-.
        vector<double> vertex_x;
        vector<double> vertex_y;
        vector<double> vertex_z;
        vector<double> vertex_chi2;
        vector<double> vertex_prob;
        vector<int>    vertex_type;

        Sub_Particle_t em;  // Electron
        Sub_Particle_t ep;  // Positron
    };

    Vertex_Particle_t v0;

    // MCParticles
    vector<double> mc_part_energy;
    vector<int>    mc_part_pdg_id;
    vector<int>    mc_part_gen_status; /** Generator Status **/
    vector<double> mc_part_time;      /** The global creation time. */
    vector<double> mc_part_x;      /** The X vertex. */
    vector<double> mc_part_y;      /** The Y vertex. */
    vector<double> mc_part_z;      /** The Z vertex. */
    vector<double> mc_part_end_x;  /** The X end point. */
    vector<double> mc_part_end_y;  /** The Y end point. */
    vector<double> mc_part_end_z;  /** The Z end point. */
    vector<double> mc_part_px;
    vector<double> mc_part_py;
    vector<double> mc_part_pz;
    vector<double> mc_part_mass;
    vector<double> mc_part_charge;
    vector< vector<int> > mc_part_daughters;
    vector< vector<int> > mc_part_parents;

public:
    MiniDst(): md_output_file_name("minidst.root"){};
    MiniDst(string_view output_file_name): md_output_file_name(output_file_name){};
    ~MiniDst(){};
    virtual void Start();
    virtual long Run(int nevt=0);
    virtual void End();
    virtual void SetOutputFileName(const string& outfile){md_output_file_name=outfile;};
    virtual void SetDebugLevel(const int level){ md_Debug = level;};
    void clear();  // Clear all the vectors.
    //std::map<std::string, vector<double>* > &Get_brmap(){return branch_map;};
    // #define FULL_CLEAR(xxx)  { for( auto p: xxx){delete p;}; xxx.clear(); };

    template<typename T> inline void FULL_CLEAR(T& xxx){
        for(auto p: xxx){ delete p; };
        xxx.clear();
    }

};


#endif //MINIDST_MINIDST_H
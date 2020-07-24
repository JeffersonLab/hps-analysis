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

using namespace std;

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

    /// All the items that we could possibly write out are here so that direct access
    /// to the variables is possible throughout the code.

    int trigger{0};
    int run_number{0};
    int event_number{0};

    // Switches that allow turning output on/off.
    bool write_ecal_cluster{true};
    bool write_ecal_hits{true};
    bool write_svt_hits{true};
    bool write_tracks{true};
    bool write_only_gbl_tracks{false};
    bool write_mc_particles{false};

    // These determine what types of particles will be written out.
    vector<int> particle_types_single{FINAL_STATE_PARTICLE, OTHER_ELECTRONS}; // No vertex
    vector<int> particle_types_double{TC_V0_CANDIDATE, UC_VC_CANDIDATE}; // Yes vertex.

    // Vector type helpers.
    Multi_Branch branch_map;

    // Ecal Hits
    vector<double> ecal_hit_energy;
    vector<double> ecal_hit_time;
    vector<int>    ecal_hit_index_x;
    vector<int>    ecal_hit_index_y;

    // Ecal Clusters:
    vector<double> ecal_cluster_energy;
    vector<double> ecal_cluster_time;
    vector<double> ecal_cluster_x;
    vector<double> ecal_cluster_y;
    vector<double> ecal_cluster_z;
    vector<int> ecal_cluster_seed_idx;
    vector< vector<int> > ecal_cluster_hits;
    vector<int> ecal_cluster_nhits; // Not strictly needed, but handy. (could use ecal_cluster_hits[i].size())

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

    // Final State Particles.
    vector<int>    part_type;
    vector<double> part_energy;
    vector<int>    part_pdg;
    vector<int>    part_charge;
    vector<double> part_goodness_of_pid;
    vector<double> part_px;
    vector<double> part_py;
    vector<double> part_pz;
    vector<double> part_corr_px;
    vector<double> part_corr_py;
    vector<double> part_corr_pz;
    vector<double> part_vertex_x;
    vector<double> part_vertex_y;
    vector<double> part_vertex_z;
    vector<double> part_vertex_chi2;
    vector<int>    part_track;       // At most one track per particle.
    vector<int>    part_ecal_cluster; // At most one ecal cluster per particle.


    // Target Constrained V0 Particles;
    vector<int>    v0_type;
    vector<double> v0_energy;
    vector<double> v0_mass;
    vector<int>    v0_pdg;
    vector<int>    v0_charge;
    vector<double> v0_goodness_of_pid;
    vector<double> v0_px;
    vector<double> v0_py;
    vector<double> v0_pz;
    vector<double> v0_corr_px;
    vector<double> v0_corr_py;
    vector<double> v0_corr_pz;
    vector<double> v0_vertex_x;
    vector<double> v0_vertex_y;
    vector<double> v0_vertex_z;
    vector<double> v0_vertex_chi2;
    vector<double> v0_n_daughter;
//    vector< vector< int> > v0_parts;
//    vector< vector< int> > v0_tracks;
//    vector< vector< int> > v0_ecal_clusters;

///
/// The following are for convenience, essentially duplicating information
/// from the particle, track or cluster collections.
///
/// Since this is HPS we are interested in the e+ and e- vertexes.
    vector<int>    v0_em_part;
    vector<int>    v0_ep_part;
    vector<int>    v0_em_track;
    vector<int>    v0_ep_track;
    vector<int>    v0_em_track_nhit;
    vector<int>    v0_ep_track_nhit;
    vector<double> v0_em_p;
    vector<double> v0_ep_p;
    vector<double> v0_em_chi2;
    vector<double> v0_ep_chi2;
    vector<double> v0_em_good_pid;
    vector<double> v0_ep_good_pid;
    vector<double> v0_em_track_time;
    vector<double> v0_ep_track_time;
    vector<double> v0_em_clus_time;
    vector<double> v0_ep_clus_time;
    vector<double> v0_em_pos_ecal_x;
    vector<double> v0_em_pos_ecal_y;
    vector<double> v0_ep_pos_ecal_x;
    vector<double> v0_ep_pos_ecal_y;

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
    MiniDst(){
        cout << "MiniDstLib starting up. \n";
    };
    MiniDst(string_view output_file_name): md_output_file_name(output_file_name){};
    ~MiniDst(){};

    void Setup();

};


#endif //MINIDST_MINIDST_H

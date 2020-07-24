//
// Created by Maurik Holtrop on 7/24/20.
//

#include "MiniDst.h"

void MiniDst::Setup() {
    // Create the output TTree
    if( !md_output_file){
        if(md_Debug & kDebug_Info) std::cout << "Opening output file: " << md_output_file_name << std::endl;
        md_output_file = new TFile(md_output_file_name.data(),"RECREATE");
    }
    else std::cout << "The output file should not be open already!\n\n\n\n";

    md_output_tree = new TTree("MiniDST","HPS mini-DST");

    md_output_tree->Branch("trigger",&trigger);
    md_output_tree->Branch("run_number",&run_number);
    md_output_tree->Branch("event_number", &event_number);

    // You do not HAVE to create the handles to the vectors in "this":
    //    branch_map.try_emplace("hidden_double", new std::vector<double>());
    // It is just a bit more cumbersome to use, plus costs you a map lookup every time you use a variable from the tree.
    //    std::get<std::vector<double>* >(branch_map["hidden_double"])->push_back(run_number);

    // Thus this is the preferred way to add variables to the output tree.
    // Doing this, rather than direct output_tree->Branch(...), allow for automatic clear() and other operations.
    // This way, the output tree is also automatically sorted alphabetically by name :-)

    if( write_ecal_hits) {
        branch_map.try_emplace("ecal_hit_index_x", &ecal_hit_index_x);
        branch_map.try_emplace("ecal_hit_index_y", &ecal_hit_index_y);
        branch_map.try_emplace("ecal_hit_energy", &ecal_hit_energy);
        branch_map.try_emplace("ecal_hit_time", &ecal_hit_time);
    }

    if( write_ecal_cluster) {
        branch_map.try_emplace("ecal_cluster_energy", &ecal_cluster_energy);
        branch_map.try_emplace("ecal_cluster_time", &ecal_cluster_time);
        branch_map.try_emplace("ecal_cluster_x", &ecal_cluster_x);
        branch_map.try_emplace("ecal_cluster_y", &ecal_cluster_y);
        branch_map.try_emplace("ecal_cluster_z", &ecal_cluster_z);
        branch_map.try_emplace("ecal_cluster_seed_idx", &ecal_cluster_seed_idx);
        branch_map.try_emplace("ecal_cluster_hits", &ecal_cluster_hits);
        branch_map.try_emplace("ecal_cluster_nhits", &ecal_cluster_nhits);
    }

    if( write_svt_hits ){
        branch_map.try_emplace("svt_hit_layer", &svt_hit_layer);
        branch_map.try_emplace("svt_hit_x", &svt_hit_x);
        branch_map.try_emplace("svt_hit_y", &svt_hit_y);
        branch_map.try_emplace("svt_hit_z", &svt_hit_z);
        branch_map.try_emplace("svt_hit_cxx", &svt_hit_cxx);
        branch_map.try_emplace("svt_hit_cxy", &svt_hit_cxy);
        branch_map.try_emplace("svt_hit_cxz", &svt_hit_cxz);
        branch_map.try_emplace("svt_hit_cyy", &svt_hit_cyy);
        branch_map.try_emplace("svt_hit_cyz", &svt_hit_cyz);
        branch_map.try_emplace("svt_hit_czz", &svt_hit_czz);
        branch_map.try_emplace("svt_hit_time", &svt_hit_time);
    }

    if( write_tracks){
        branch_map.try_emplace("track_n_hits", & track_n_hits); /** The number of 3D hits associated with this track. */
        branch_map.try_emplace("track_volume", & track_volume); /** The volume to which this track belongs to. */
        branch_map.try_emplace("track_type", & track_type);   /** The track type. */
        branch_map.try_emplace("track_d0", &track_d0);   /** The distance of closest approach to the reference point. */
        branch_map.try_emplace("track_phi0", &track_phi0);     //* The azimuthal angle of the momentum at the position of closest approach to the reference point.
        branch_map.try_emplace("track_omega", &track_omega); //* The track curvature. The curvature is positive (negative) if the particle has a positive (negative) charge.
        branch_map.try_emplace("track_tan_lambda", &track_tan_lambda); // The slope of the track in the SY plane where S is the arc length of the helix in the xz plane.
        branch_map.try_emplace("track_z0", &track_z0);      // The y position of the track at the distance of closest approach in the xz plane.
        branch_map.try_emplace("track_chi2", &track_chi2); /** The chi^2 of the track fit. */
        branch_map.try_emplace("track_time", &track_time); /** The time of the track.  This is currently the average time of all hits composing the track. **/
        branch_map.try_emplace("track_x_at_ecal", &track_x_at_ecal); /** The x position of the extrapolated track at the Ecal face. */
        branch_map.try_emplace("track_y_at_ecal", &track_y_at_ecal); /** The y position of the extrapolated track at the Ecal face. */
        branch_map.try_emplace("track_z_at_ecal", &track_z_at_ecal); /** The z position of the extrapolated track at the Ecal face. */
        branch_map.try_emplace("track_isolation", &track_isolation); /** Array used to store the isolation variables for each of the sensor layers. */
        branch_map.try_emplace("track_covmatrix", &track_covmatrix);   /** The 1/2 Covariant Matrix. This is the lower 1/2. **/
        branch_map.try_emplace("track_lambda_kinks", &track_lambda_kinks);
        branch_map.try_emplace("track_phi_kinks", &track_phi_kinks);
        branch_map.try_emplace("track_particle", & track_particle); /** Reference to the reconstructed particle associated with this track. */
        branch_map.try_emplace("track_svt_hits",&track_svt_hits);/** Reference to the 3D hits associated with this track. */
        branch_map.try_emplace("track_gbl_ref", & track_gbl_ref);
        branch_map.try_emplace("track_ref", & track_ref);
    }

    branch_map.try_emplace("part_type", &part_type);
    branch_map.try_emplace("part_energy", &part_energy);
    branch_map.try_emplace("part_pdg", &part_pdg);
    branch_map.try_emplace("part_charge", &part_charge);
    branch_map.try_emplace("part_goodness_of_pid", &part_goodness_of_pid);
    branch_map.try_emplace("part_px", &part_px);
    branch_map.try_emplace("part_py", &part_py);
    branch_map.try_emplace("part_pz", &part_pz);
    branch_map.try_emplace("part_corr_px", &part_corr_px);
    branch_map.try_emplace("part_corr_py", &part_corr_py);
    branch_map.try_emplace("part_corr_pz", &part_corr_pz);
    branch_map.try_emplace("part_vertex_x", &part_vertex_x);
    branch_map.try_emplace("part_vertex_y", &part_vertex_y);
    branch_map.try_emplace("part_vertex_z", &part_vertex_z);
    branch_map.try_emplace("part_vertex_chi2", &part_vertex_chi2); // Always zero.
    branch_map.try_emplace("part_track", &part_track);
    branch_map.try_emplace("part_ecal_cluster", &part_ecal_cluster);

    branch_map.try_emplace("v0_type", &v0_type);          // Always 3, but should be able to check.
    branch_map.try_emplace("v0_energy", &v0_energy);
    branch_map.try_emplace("v0_mass", &v0_mass);
    branch_map.try_emplace("v0_px", &v0_px);
    branch_map.try_emplace("v0_py", &v0_py);
    branch_map.try_emplace("v0_pz", &v0_pz);
    branch_map.try_emplace("v0_corr_px", &v0_corr_px);
    branch_map.try_emplace("v0_corr_py", &v0_corr_py);
    branch_map.try_emplace("v0_corr_pz", &v0_corr_pz);
    branch_map.try_emplace("v0_vertex_x", &v0_vertex_x);
    branch_map.try_emplace("v0_vertex_y", &v0_vertex_y);
    branch_map.try_emplace("v0_vertex_z", &v0_vertex_z);
    branch_map.try_emplace("v0_vertex_chi2", &v0_vertex_chi2);
    branch_map.try_emplace("v0_n_daughter", &v0_n_daughter);
//    branch_map.try_emplace("v0_parts", &v0_parts);
//    branch_map.try_emplace("v0_tracks", &v0_tracks);
//    branch_map.try_emplace("v0_ecal_clusters", &v0_ecal_clusters);

    branch_map.try_emplace("v0_em_part", &v0_em_part);
    branch_map.try_emplace("v0_ep_part", &v0_ep_part);
    branch_map.try_emplace("v0_ep_track", &v0_ep_track);
    branch_map.try_emplace("v0_em_track", &v0_em_track);
    branch_map.try_emplace("v0_em_track_nhit", &v0_em_track_nhit);
    branch_map.try_emplace("v0_ep_track_nhit", &v0_ep_track_nhit);

    branch_map.try_emplace("v0_em_p", &v0_em_p);
    branch_map.try_emplace("v0_ep_p", &v0_ep_p);

    branch_map.try_emplace("v0_ep_chi2", &v0_ep_chi2);
    branch_map.try_emplace("v0_em_chi2", &v0_em_chi2);
    branch_map.try_emplace("v0_ep_good_pid", &v0_ep_good_pid);
    branch_map.try_emplace("v0_em_good_pid", &v0_em_good_pid);
    branch_map.try_emplace("v0_em_track_time", &v0_em_track_time);
    branch_map.try_emplace("v0_ep_track_time", &v0_ep_track_time);
    branch_map.try_emplace("v0_em_clus_time", &v0_em_clus_time);
    branch_map.try_emplace("v0_ep_clus_time", &v0_ep_clus_time);
    branch_map.try_emplace("v0_em_pos_ecal_x", &v0_em_pos_ecal_x);
    branch_map.try_emplace("v0_em_pos_ecal_y", &v0_em_pos_ecal_y);
    branch_map.try_emplace("v0_ep_pos_ecal_x", &v0_ep_pos_ecal_x);
    branch_map.try_emplace("v0_ep_pos_ecal_y", &v0_ep_pos_ecal_y);

    branch_map.try_emplace("v0_pdg", &v0_pdg);            // Not usefull, always zero ???
    branch_map.try_emplace("v0_charge", &v0_charge);      // Not usefull, always zero ???
    branch_map.try_emplace("v0_goodness_of_pid", &v0_goodness_of_pid); // Not usefull, always zero ???


    if(write_mc_particles) {
        // MCParticles
        branch_map.try_emplace("mc_part_energy", &mc_part_energy);
        branch_map.try_emplace("mc_part_pdg_id", &mc_part_pdg_id);
        branch_map.try_emplace("mc_part_gen_status", &mc_part_gen_status); /** Generator Status **/
        branch_map.try_emplace("mc_part_time", &mc_part_time);      /** The global creation time. */
        branch_map.try_emplace("mc_part_x", &mc_part_x);      /** The X vertex. */
        branch_map.try_emplace("mc_part_y", &mc_part_y);      /** The Y vertex. */
        branch_map.try_emplace("mc_part_z", &mc_part_z);      /** The Z vertex. */
        branch_map.try_emplace("mc_part_end_x", &mc_part_end_x);  /** The X end point. */
        branch_map.try_emplace("mc_part_end_y", &mc_part_end_y);  /** The Y end point. */
        branch_map.try_emplace("mc_part_end_z", &mc_part_end_z);  /** The Z end point. */
        branch_map.try_emplace("mc_part_px", &mc_part_px);
        branch_map.try_emplace("mc_part_py", &mc_part_py);
        branch_map.try_emplace("mc_part_pz", &mc_part_pz);
        branch_map.try_emplace("mc_part_mass", &mc_part_mass);
        branch_map.try_emplace("mc_part_charge", &mc_part_charge);
        branch_map.try_emplace("mc_part_daughters",&mc_part_daughters);
        branch_map.try_emplace("mc_part_parents", &mc_part_parents);
    }

    // Now set all the branches in the tree:
    // for( auto const& [nam, bran] : branch_map ){ // auto does not work with a Lambda, since no type declaration at all!
    for(std::map<std::string,Multi_Value>::iterator it = branch_map.begin(); it != branch_map.end(); ++it ){
        if(md_Debug & kDebug_L1) cout << "Creating branch " << it->first << endl;
        //
        // Explicit specification:
        //if(bran.index() == 0){ output_tree->Branch(nam.c_str(), std::get< vector<double>*>(bran));}
        //if(bran.index() == 1){ output_tree->Branch(nam.c_str(), std::get< vector<int>*>(bran));}
        std::visit([this,&it](auto&& arg){md_output_tree->Branch(it->first.c_str(), arg);}, it->second);
    }
}
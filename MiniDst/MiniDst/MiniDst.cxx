//
// Created by Maurik Holtrop on 7/24/20.
//

#include "MiniDst.h"

void MiniDst::Start() {
   DefineBranchMap();
   SetBranchMap();    // Opens output file and sets the branches.
   WriteStateToFile();
}

void MiniDst::WriteStateToFile() {
   /// Write the relevant state information into the ROOT file so you can later
   /// know how the file was written.
   /// Must be called *after* SetBranchMap so the output file is open.
   /// Implement this in the derived classes for it to actually do anything!
}

std::vector<std::string> MiniDst::GetBranchNames() {
   std::vector<std::string> names;
   for(auto & it : branch_map){
      names.push_back(it.first);
   }
   return names;
}

std::vector<std::string> MiniDst::GetActiveBranchNames() {
   std::vector<std::string> names;
   for(auto & it : branch_map_active){
      if(it.second) names.push_back(it.first);
   }
   return names;
}


void MiniDst::SetBranchActive(std::string name, bool active) {
   if( branch_map_active.find(name) != branch_map_active.end() ) {
      branch_map_active[name] = active;
   }else{
      cout << "ERROR: " << name << " is not in the branch map.";
   }
}

void MiniDst::DefineBranchMap(bool use_all) {
   /// Setup the branch_map and branch_map_active.
   /// Default is use_all = false. Set this to true to place everything on the branch.

   // Notes on branch_map:
   //    You do not HAVE to create the handles to the vectors in "this":
   //         branch_map_try_emplace("hidden_double", new std::vector<double>());
   //    It is just a bit more cumbersome to use, plus costs you a map lookup every time you
   //    use a variable from the tree.
   //    std::get<std::vector<double>* >(branch_map["hidden_double"])->push_back(run_number);
   //
   // Thus the preferred way is to add the variables to the output tree with a handle.
   // Doing this, rather than direct output_tree->Branch(...), allow for automatic Clear() and other operations.
   // (See note above on implementing a visitor.)
   // This way, the output tree is also automatically sorted alphabetically by name :-)
   //
   // The helper function branch_map_try_emplace does a branch_map.try_emplace() and branch_map_active.try_emplace()
   // just to shorten the code.

   branch_map_try_emplace("run_number", &run_number);
   branch_map_try_emplace("event_number", &event_number);
   branch_map_try_emplace("time_stamp", &time_stamp);
   branch_map_try_emplace("svt_status", &svt_status);
   branch_map_try_emplace("trigger", &trigger);
   branch_map_try_emplace("ext_trigger", &ext_trigger);
   branch_map_try_emplace("rf_time1", &rf_time1);
   branch_map_try_emplace("rf_time2", &rf_time2);
   branch_map_try_emplace("track_n_gbl", &track_n_gbl);

   branch_map_try_emplace("hodo_raw_ix", &hodo_raw_ix, use_hodo_raw_hits | use_all );
   branch_map_try_emplace("hodo_raw_iy", &hodo_raw_iy, use_hodo_raw_hits | use_all );
   branch_map_try_emplace("hodo_raw_hole", &hodo_raw_hole, use_hodo_raw_hits | use_all );
   branch_map_try_emplace("hodo_raw_layer", &hodo_raw_layer, use_hodo_raw_hits | use_all );
   branch_map_try_emplace("hodo_raw_adc", &hodo_raw_adc, use_hodo_raw_hits | use_all );

   branch_map_try_emplace("hodo_hit_energy", &hodo_hit_energy, use_hodo_hits | use_all );
   branch_map_try_emplace("hodo_hit_time", &hodo_hit_time, use_hodo_hits | use_all );
   branch_map_try_emplace("hodo_hit_index_x", &hodo_hit_index_x, use_hodo_hits | use_all );
   branch_map_try_emplace("hodo_hit_index_y", &hodo_hit_index_y, use_hodo_hits | use_all );
   branch_map_try_emplace("hodo_hit_hole", &hodo_hit_hole, use_hodo_hits | use_all );
   branch_map_try_emplace("hodo_hit_layer", &hodo_hit_layer, use_hodo_hits | use_all );

   branch_map_try_emplace("hodo_cluster_energy", &hodo_cluster_energy, use_hodo_clusters | use_all );
   branch_map_try_emplace("hodo_cluster_time", &hodo_cluster_time, use_hodo_clusters | use_all );
   branch_map_try_emplace("hodo_cluster_ix", &hodo_cluster_ix, use_hodo_clusters | use_all );
   branch_map_try_emplace("hodo_cluster_iy", &hodo_cluster_iy, use_hodo_clusters | use_all );
   branch_map_try_emplace("hodo_cluster_layer", &hodo_cluster_layer, use_hodo_clusters | use_all );

   branch_map_try_emplace("ecal_raw_ix", &ecal_raw_ix, use_ecal_raw_hits | use_all );
   branch_map_try_emplace("ecal_raw_iy", &ecal_raw_iy, use_ecal_raw_hits | use_all );
   branch_map_try_emplace("ecal_raw_adc", &ecal_raw_adc, use_ecal_raw_hits | use_all );

   branch_map_try_emplace("ecal_hit_energy", &ecal_hit_energy, use_ecal_hits | use_all );
   branch_map_try_emplace("ecal_hit_time", &ecal_hit_time, use_ecal_hits | use_all );
   branch_map_try_emplace("ecal_hit_index_x", &ecal_hit_index_x, use_ecal_hits | use_all );
   branch_map_try_emplace("ecal_hit_index_y", &ecal_hit_index_y, use_ecal_hits | use_all );
   branch_map_try_emplace("ecal_hit_x", &ecal_hit_x, use_ecal_hits | use_all );
   branch_map_try_emplace("ecal_hit_y", &ecal_hit_y, use_ecal_hits | use_all );
   branch_map_try_emplace("ecal_hit_z", &ecal_hit_z, use_ecal_hits | use_all );

   branch_map_try_emplace("ecal_hit_mc_contrib_id", &ecal_hit_mc_contrib_id, (use_ecal_hits & use_mc_particles) | use_all);
   branch_map_try_emplace("ecal_hit_mc_contrib_pdg",  &ecal_hit_mc_contrib_pdg, (use_ecal_hits & use_mc_particles) | use_all);
   branch_map_try_emplace("ecal_hit_mc_contrib_ec", &ecal_hit_mc_contrib_ec, (use_ecal_hits & use_mc_particles) | use_all);
   branch_map_try_emplace("ecal_hit_mc_parent_id",  &ecal_hit_mc_parent_id, (use_ecal_hits & use_mc_particles) | use_all);
   branch_map_try_emplace("ecal_hit_mc_parent_pdg", &ecal_hit_mc_parent_pdg, (use_ecal_hits & use_mc_particles) | use_all);

   branch_map_try_emplace("ecal_cluster_energy", &ecal_cluster_energy, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_time", &ecal_cluster_time, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_x", &ecal_cluster_x, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_y", &ecal_cluster_y, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_z", &ecal_cluster_z, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_seed_index", &ecal_cluster_seed_index, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_seed_ix", &ecal_cluster_seed_ix, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_seed_iy", &ecal_cluster_seed_iy, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_seed_energy", &ecal_cluster_seed_energy, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_hits", &ecal_cluster_hits, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_nhits", &ecal_cluster_nhits, use_ecal_cluster | use_all );
   branch_map_try_emplace("ecal_cluster_mc_id", &ecal_cluster_mc_id, (use_ecal_cluster & use_mc_particles) | use_all);
   branch_map_try_emplace("ecal_cluster_mc_pdg", &ecal_cluster_mc_pdg, (use_ecal_cluster & use_mc_particles) | use_all);
   branch_map_try_emplace("ecal_cluster_mc_pdg_purity", &ecal_cluster_mc_pdg_purity, (use_ecal_cluster & use_mc_particles) | use_all);

   branch_map_try_emplace("ecal_cluster_uncor_energy", &ecal_cluster_uncor_energy, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_time", &ecal_cluster_uncor_time, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_x", &ecal_cluster_uncor_x, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_y", &ecal_cluster_uncor_y, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_z", &ecal_cluster_uncor_z, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_seed_index", &ecal_cluster_uncor_seed_index, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_seed_ix", &ecal_cluster_uncor_seed_ix, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_seed_iy", &ecal_cluster_uncor_seed_iy, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_seed_energy", &ecal_cluster_uncor_seed_energy, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_hits", &ecal_cluster_uncor_hits, use_ecal_cluster_uncor | use_all );
   branch_map_try_emplace("ecal_cluster_uncor_nhits", &ecal_cluster_uncor_nhits, use_ecal_cluster_uncor | use_all );


   branch_map_try_emplace("svt_raw_hit_layer", &svt_raw_hit_layer, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_module", &svt_raw_hit_module, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_strip", &svt_raw_hit_strip, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_adc", &svt_raw_hit_adc, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_t0", &svt_raw_hit_t0, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_t0_err", &svt_raw_hit_t0_err, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_amp", &svt_raw_hit_amp, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_amp_err", &svt_raw_hit_amp_err, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_chi2", &svt_raw_hit_chi2, use_svt_raw_hits | use_all );
   branch_map_try_emplace("svt_raw_hit_fit_no", &svt_raw_hit_fit_no, use_svt_raw_hits | use_all );

   branch_map_try_emplace("svt_hit_type", &svt_hit_type, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_x", &svt_hit_x, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_y", &svt_hit_y, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_z", &svt_hit_z, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_cxx", &svt_hit_cxx, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_cxy", &svt_hit_cxy, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_cxz", &svt_hit_cxz, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_cyy", &svt_hit_cyy, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_cyz", &svt_hit_cyz, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_czz", &svt_hit_czz, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_time", &svt_hit_time, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_edep", &svt_hit_edep, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_raw_index", &svt_hit_raw_index, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_raw_other", &svt_hit_raw_other, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_layer", &svt_hit_layer, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_module", &svt_hit_module, use_svt_hits | use_all );
   branch_map_try_emplace("svt_hit_strip", &svt_hit_strip, use_svt_hits | use_all );


   bool write_any_tracks = use_kf_tracks || use_gbl_tracks || use_matched_tracks || use_all;
   branch_map_try_emplace("track_n_kf", &track_n_kf, write_any_tracks);
   branch_map_try_emplace("track_n_gbl", &track_n_gbl, write_any_tracks);
   branch_map_try_emplace("track_n_matched", &track_n_matched, write_any_tracks);
   branch_map_try_emplace("track_n_hits", &track_n_hits, write_any_tracks);
   branch_map_try_emplace("track_volume", &track_volume, write_any_tracks);
   branch_map_try_emplace("track_type", &track_type, write_any_tracks);
   branch_map_try_emplace("track_d0", &track_d0, write_any_tracks);
   branch_map_try_emplace("track_phi0",&track_phi0, write_any_tracks);
   branch_map_try_emplace("track_omega",&track_omega, write_any_tracks);
   branch_map_try_emplace("track_tan_lambda",&track_tan_lambda, write_any_tracks);
   branch_map_try_emplace("track_z0",&track_z0, write_any_tracks);
   branch_map_try_emplace("track_chi2", &track_chi2, write_any_tracks);
   branch_map_try_emplace("track_time",&track_time, write_any_tracks);
   branch_map_try_emplace("track_px",&track_px, write_any_tracks);
   branch_map_try_emplace("track_py",&track_py, write_any_tracks);
   branch_map_try_emplace("track_pz",&track_pz, write_any_tracks);
   branch_map_try_emplace("track_x_at_target",&track_x_at_target, write_any_tracks);
   branch_map_try_emplace("track_y_at_target",&track_y_at_target, write_any_tracks);
   branch_map_try_emplace("track_z_at_target",&track_z_at_target, write_any_tracks);
#
#ifdef DEBUG
   branch_map_try_emplace("track_px_old",&track_px_old, write_any_tracks);
   branch_map_try_emplace("track_py_old",&track_py_old, write_any_tracks);
   branch_map_try_emplace("track_pz_old",&track_pz_old, write_any_tracks);
   branch_map_try_emplace("track_bfield_at_target",&track_bfield_at_target, write_any_tracks);
   branch_map_try_emplace("track_bfield_alt",&track_bfield_alt, write_any_tracks);
   branch_map_try_emplace("track_omega_old",&track_omega_old, write_any_tracks);
   branch_map_try_emplace("track_bfield_at_lasthit",&track_bfield_at_lasthit,write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_bfield_at_ecal",&track_bfield_at_ecal,write_any_tracks && use_extra_tracks);
#endif
   branch_map_try_emplace("track_x_at_lasthit",&track_x_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_y_at_lasthit",&track_y_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_z_at_lasthit",&track_z_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_px_at_lasthit",&track_px_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_py_at_lasthit",&track_py_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_pz_at_lasthit",&track_pz_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_omega_at_lasthit",&track_omega_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_tan_lambda_at_lasthit",&track_tan_lambda_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_phi0_at_lasthit",&track_phi0_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_d0_at_lasthit",&track_d0_at_lasthit, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_z0_at_lasthit",&track_z0_at_lasthit, write_any_tracks && use_extra_tracks);

   branch_map_try_emplace("track_x_at_ecal",&track_x_at_ecal, write_any_tracks);
   branch_map_try_emplace("track_y_at_ecal",&track_y_at_ecal, write_any_tracks);
   branch_map_try_emplace("track_z_at_ecal",&track_z_at_ecal, write_any_tracks);
   branch_map_try_emplace("track_px_at_ecal",&track_px_at_ecal, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_py_at_ecal",&track_py_at_ecal, write_any_tracks && use_extra_tracks);
   branch_map_try_emplace("track_pz_at_ecal",&track_pz_at_ecal, write_any_tracks && use_extra_tracks);

   branch_map_try_emplace("track_isolation",&track_isolation, write_any_tracks);
   branch_map_try_emplace("track_covmatrix",&track_covmatrix, write_any_tracks);
   branch_map_try_emplace("track_lambda_kinks", &track_lambda_kinks, write_any_tracks);
   branch_map_try_emplace("track_phi_kinks", &track_phi_kinks, write_any_tracks);
   branch_map_try_emplace("track_particle",&track_particle, write_any_tracks);
   branch_map_try_emplace("track_gbl_ref", &track_gbl_ref, write_any_tracks);
   branch_map_try_emplace("track_ref", &track_ref, write_any_tracks);
   branch_map_try_emplace("track_svt_hits",&track_svt_hits, write_any_tracks);
   branch_map_try_emplace("track_hit_pattern",&track_hit_pattern, write_any_tracks);

   bool write_particles = use_kf_particles || use_gbl_particles;
   branch_map_try_emplace("part_type", &part.type, write_particles);
   branch_map_try_emplace("part_lcio_type", &part.lcio_type, write_particles);
   branch_map_try_emplace("part_energy", &part.energy, write_particles);
   branch_map_try_emplace("part_mass", &part.mass, write_particles);
   branch_map_try_emplace("part_pdg", &part.pdg, write_particles);
   branch_map_try_emplace("part_charge", &part.charge, write_particles);
   branch_map_try_emplace("part_goodness_of_pid", &part.goodness_of_pid, write_particles);
   branch_map_try_emplace("part_px", &part.px, write_particles);
   branch_map_try_emplace("part_py", &part.py, write_particles);
   branch_map_try_emplace("part_pz", &part.pz, write_particles);
   branch_map_try_emplace("part_track", &part.track, write_particles);
   branch_map_try_emplace("part_track_chi2", &part.track_chi2, write_particles);
   branch_map_try_emplace("part_ecal_cluster", &part.ecal_cluster, write_particles);

   branch_map_try_emplace("v0_type", &v0.type, write_particles);
   branch_map_try_emplace("v0_lcio_type", &v0.lcio_type, write_particles);
   branch_map_try_emplace("v0_energy", &v0.energy, write_particles);
   branch_map_try_emplace("v0_mass", &v0.mass, write_particles);
   branch_map_try_emplace("v0_mass_err", &v0.mass_err, write_particles);

   branch_map_try_emplace("v0_pdg", &v0.pdg, write_particles);
   branch_map_try_emplace("v0_charge", &v0.charge, write_particles);
   branch_map_try_emplace("v0_goodness_of_pid", &v0.goodness_of_pid, write_particles);

   branch_map_try_emplace("v0_p", &v0.p, write_particles);
   branch_map_try_emplace("v0_p_err", &v0.p_err, write_particles);

   branch_map_try_emplace("v0_px", &v0.px, write_particles);
   branch_map_try_emplace("v0_py", &v0.py, write_particles);
   branch_map_try_emplace("v0_pz", &v0.pz, write_particles);
   branch_map_try_emplace("v0_px_err", &v0.px_err, write_particles);
   branch_map_try_emplace("v0_py_err", &v0.py_err, write_particles);
   branch_map_try_emplace("v0_pz_err", &v0.pz_err, write_particles);

   branch_map_try_emplace("v0_vertex_x", &v0.vertex_x, write_particles);
   branch_map_try_emplace("v0_vertex_y", &v0.vertex_y, write_particles);
   branch_map_try_emplace("v0_vertex_z", &v0.vertex_z, write_particles);
   branch_map_try_emplace("v0_vertex_x_err", &v0.vertex_x, write_particles);
   branch_map_try_emplace("v0_vertex_y_err", &v0.vertex_y, write_particles);
   branch_map_try_emplace("v0_vertex_z_err", &v0.vertex_z, write_particles);

   branch_map_try_emplace("v0_vertex_chi2", &v0.vertex_chi2, write_particles);
   branch_map_try_emplace("v0_vertex_prob", &v0.vertex_prob, write_particles);

   branch_map_try_emplace("v0_em_part", &v0.em.part, write_particles);
   branch_map_try_emplace("v0_em_track", &v0.em.track, write_particles);
   branch_map_try_emplace("v0_em_track_nhit", &v0.em.track_nhit, write_particles);
   branch_map_try_emplace("v0_em_p", &v0.em.p, write_particles);
   branch_map_try_emplace("v0_em_chi2", &v0.em.chi2, write_particles);
   branch_map_try_emplace("v0_em_good_pid", &v0.em.good_pid, write_particles);
   branch_map_try_emplace("v0_em_track_time", &v0.em.track_time, write_particles);
   branch_map_try_emplace("v0_em_pos_ecal_x", &v0.em.pos_ecal_x, write_particles);
   branch_map_try_emplace("v0_em_pos_ecal_y", &v0.em.pos_ecal_y, write_particles);
   branch_map_try_emplace("v0_em_clus", &v0.em.clus, write_particles);
   branch_map_try_emplace("v0_em_clus_energy", &v0.em.clus_energy, write_particles);
   branch_map_try_emplace("v0_em_clus_time", &v0.em.clus_time, write_particles);
   branch_map_try_emplace("v0_em_clus_ix", &v0.em.clus_ix, write_particles);
   branch_map_try_emplace("v0_em_clus_iy", &v0.em.clus_iy, write_particles);
   branch_map_try_emplace("v0_em_clus_pos_x", &v0.em.clus_pos_x, write_particles);
   branch_map_try_emplace("v0_em_clus_pos_y", &v0.em.clus_pos_y, write_particles);
   branch_map_try_emplace("v0_em_px", &v0.em_px_refit, write_particles);
   branch_map_try_emplace("v0_em_py", &v0.em_py_refit, write_particles);
   branch_map_try_emplace("v0_em_pz", &v0.em_pz_refit, write_particles);

   branch_map_try_emplace("v0_ep_part", &v0.ep.part, write_particles);
   branch_map_try_emplace("v0_ep_track", &v0.ep.track, write_particles);
   branch_map_try_emplace("v0_ep_track_nhit", &v0.ep.track_nhit, write_particles);
   branch_map_try_emplace("v0_ep_p", &v0.ep.p, write_particles);
   branch_map_try_emplace("v0_ep_chi2", &v0.ep.chi2, write_particles);
   branch_map_try_emplace("v0_ep_good_pid", &v0.ep.good_pid, write_particles);
   branch_map_try_emplace("v0_ep_track_time", &v0.ep.track_time, write_particles);
   branch_map_try_emplace("v0_ep_pos_ecal_x", &v0.ep.pos_ecal_x, write_particles);
   branch_map_try_emplace("v0_ep_pos_ecal_y", &v0.ep.pos_ecal_y, write_particles);
   branch_map_try_emplace("v0_ep_clus", &v0.ep.clus, write_particles);
   branch_map_try_emplace("v0_ep_clus_energy", &v0.ep.clus_energy, write_particles);
   branch_map_try_emplace("v0_ep_clus_time", &v0.ep.clus_time, write_particles);
   branch_map_try_emplace("v0_ep_clus_ix", &v0.ep.clus_ix, write_particles);
   branch_map_try_emplace("v0_ep_clus_iy", &v0.ep.clus_iy, write_particles);
   branch_map_try_emplace("v0_ep_clus_pos_x", &v0.ep.clus_pos_x, write_particles);
   branch_map_try_emplace("v0_ep_clus_pos_y", &v0.ep.clus_pos_y, write_particles);

   branch_map_try_emplace("v0_ep_px", &v0.ep_px_refit, write_particles);
   branch_map_try_emplace("v0_ep_py", &v0.ep_py_refit, write_particles);
   branch_map_try_emplace("v0_ep_pz", &v0.ep_pz_refit, write_particles);


   // MCParticles
   branch_map_try_emplace("mc_part_energy", &mc_part_energy, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_pdg", &mc_part_pdg, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_id", &mc_part_id, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_sim_status", &mc_part_sim_status, use_mc_particles | use_all ); /** Generator Status **/
   branch_map_try_emplace("mc_part_time", &mc_part_time, use_mc_particles | use_all );      /** The global creation time. */
   branch_map_try_emplace("mc_part_x", &mc_part_x, use_mc_particles | use_all );      /** The X vertex. */
   branch_map_try_emplace("mc_part_y", &mc_part_y, use_mc_particles | use_all );      /** The Y vertex. */
   branch_map_try_emplace("mc_part_z", &mc_part_z, use_mc_particles | use_all );      /** The Z vertex. */
   branch_map_try_emplace("mc_part_px", &mc_part_px, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_py", &mc_part_py, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_pz", &mc_part_pz, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_end_x", &mc_part_end_x, use_mc_particles | use_all );  /** The X end point. */
   branch_map_try_emplace("mc_part_end_y", &mc_part_end_y, use_mc_particles | use_all );  /** The Y end point. */
   branch_map_try_emplace("mc_part_end_z", &mc_part_end_z, use_mc_particles | use_all );  /** The Z end point. */
//   branch_map_try_emplace("mc_part_end_px", &mc_part_end_px, use_mc_particles | use_all );  /** The X end point. */
//   branch_map_try_emplace("mc_part_end_py", &mc_part_end_py, use_mc_particles | use_all );  /** The Y end point. */
//   branch_map_try_emplace("mc_part_end_pz", &mc_part_end_pz, use_mc_particles | use_all );  /** The Z end point. */

   branch_map_try_emplace("mc_part_mass", &mc_part_mass, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_charge", &mc_part_charge, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_daughters", &mc_part_daughters, use_mc_particles | use_all );
   branch_map_try_emplace("mc_part_parents", &mc_part_parents, use_mc_particles | use_all );

   branch_map_try_emplace("mc_score_type", &mc_score_type, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_part_idx", &mc_score_part_idx, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_x", &mc_score_x, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_y", &mc_score_y, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_z", &mc_score_z, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_px", &mc_score_px, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_py", &mc_score_py, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_pz", &mc_score_pz, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_time", &mc_score_time, use_mc_scoring | use_all );
   branch_map_try_emplace("mc_score_pdg", &mc_score_pdg, use_mc_scoring | use_all );
}

void MiniDst::SetBranchMap() {
   /// Open the output file, create the output TTree and set the branches in the TTree (if active)

   // Create the output TTree
   if (!md_output_file) {
      if (md_Debug & kDebug_Info) std::cout << "Opening output file: " << md_output_file_name << std::endl;
      md_output_file = new TFile(md_output_file_name.data(), "RECREATE","MiniDST");
   } else std::cout << "The output file should not be open already!\n\n\n\n";

   md_output_tree = new TTree("MiniDST", "HPS mini-DST");

   // Now set all the branches in the tree:
   // for( auto const& [nam, bran] : branch_map ){ // <-- Cannot do this, "auto" does not work with a Lambda, since no type declaration at all!
   //
   for(std::map<std::string,Multi_Value>::iterator it = branch_map.begin(); it != branch_map.end(); ++it ){  // Iterate over the map.

      if(branch_map_active[it->first]) {
         //
         // Simple visitor works here since each argument has exactly the same operation.
         //
         std::visit([this, &it](auto &&arg) { md_output_tree->Branch(it->first.c_str(), arg); }, it->second);
         if(md_Debug & kDebug_L1) cout << "Branch " << it->first << " is active\n";
      }
   }
}

// The following templates are needed to get an overloaded visitor.
// helper type for the visitor
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
// explicit deduction guide (not needed as of C++20)
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

void MiniDst::SetBranchAddressesOnTree(TTree *read_tree) {
   /// Connect the variables/vectors to the TTree for *reading* the tree.

   // Get the list of branches on the tree, so we can check if we can set the addresses.
   TObjArray *branch_list = read_tree->GetListOfBranches();

   for(std::map<std::string,Multi_Value>::iterator it = branch_map.begin(); it != branch_map.end(); ++it ){  // Iterate over the map.

      if(md_Debug & kDebug_L1) cout << "Looking at setting the branch for " << it->first << " which is " << branch_map_active[it->first] << endl;

      if(branch_map_active[it->first]){ // Only the active branches.
         if( branch_list->FindObject(it->first.c_str()) != nullptr){  // Check that the tree being read has the variable
            // We need an overloaded visitor:
            // For primitive, we set the branch to the address of the primitive:
            //             int run_number;
            //             t->SetBranchAddress("run_number", &run_number);
            // For STL containers, we set the branch to the address of the pointer to the vector<>:
            //             vector<int> *part_pdg;
            //             t->SetBranchAddress("part_pdg", &part_pdg);
            std::visit( overloaded{
                  [&read_tree, &it, this](int *arg){
                     if(md_Debug & kDebug_L1) cout << "Setting " << it->first.c_str() << " int arg_ptr: " << arg << endl;
                     read_tree->SetBranchAddress(it->first.c_str(), arg);
                  },
                  [&read_tree, &it, this](unsigned int *arg){
                     if(md_Debug & kDebug_L1) cout << "Setting " << it->first.c_str() << " uns int arg_ptr: " << arg << endl;
                     read_tree->SetBranchAddress(it->first.c_str(), arg);},
                  [&read_tree, &it, this](double *arg){
                     if(md_Debug & kDebug_L1) cout << "Setting " << it->first.c_str() << " double  arg_ptr: " << arg << endl;
                     read_tree->SetBranchAddress(it->first.c_str(), arg);},
                  [&read_tree, &it, this](ULong64_t *arg){
                     if(md_Debug & kDebug_L1) cout << "Setting " << it->first.c_str() << " ulong   arg_ptr: " << arg << endl;
                     read_tree->SetBranchAddress(it->first.c_str(), arg);},
                  [&read_tree, &it, this](auto &&arg){
                     if(md_Debug & kDebug_L1) cout << "Setting " << it->first.c_str() << " vector  arg_ptr: " << arg << " with address: " << &arg << endl;
                     read_tree->SetBranchAddress(it->first.c_str(), &arg);}    // All other's are vectors
            }, it->second);
         }
      }
   }
}

void MiniDst::Clear(){
   // Clear the event storage vectors.
   // We use the map to Clear all the vectors. Note that this means you *must* have each vector in the branch_map,
   // otherwise you will create an ever growing vector.

   for( auto const& [nam, bran] : branch_map ){
      // std::visit([](auto &&arg){arg->clear();},bran);  Using simple lambda if the variant only contains vectors.
      //
      // Instead, we use an overloaded lambda. See: https://en.cppreference.com/w/cpp/utility/variant/visit
      // Also: https://arne-mertz.de/2018/05/overload-build-a-variant-visitor-on-the-fly/
      // We need one line for each type that is contained in our variant.
      //
      std::visit(overloaded{
            [](int *arg)    { (*arg) = 0; },
            [](unsigned int *arg)    { (*arg) = 0; },
            [](double *arg)    { (*arg) = 0; },
            [](ULong64_t *arg)    { (*arg) = 0; },
            [](vector<int> *arg)    { arg->clear(); },
            [](vector<double> *arg) { arg->clear(); },
            [](vector< vector<short> > *arg) { arg->clear(); },
            [](vector< vector<int> > *arg) { arg->clear(); },
            [](vector< vector<double> >*arg) {arg->clear(); },
      },bran);
   }
   v0.clear();
   part.clear();
}

long MiniDst::Run(int nevt){
   std::cout << "MiniDst::Run() called. =====================<<<<<<<<<<<<<======\n";
   return(0);
}

bool MiniDst::Process(Long64_t entry){
   /// Process a single event.
   /// Intended to be overridden by sub-classes that have an event loop.

   // For vanilla MiniDst there is nothing to process.
  return true;
}

void MiniDst::End() {
   if (md_output_file->IsOpen()) {
      if (!md_output_file->IsWritable()) {
         cout << "End -- file is open but not writable!!??\n";
      } else {
         md_output_file->Write();
         md_output_file->Close();
      }
   }
}

int MiniDst::Add_Ecal_Hit(MiniDst &event_in, int i_hit) {
   /// Copy an ECal hit from event_in at i_hit and add it to the current MiniDst event.
   ///
   ecal_hit_energy.push_back(event_in.ecal_hit_energy[i_hit]);
   ecal_hit_x.push_back(event_in.ecal_hit_x[i_hit]);
   ecal_hit_y.push_back(event_in.ecal_hit_y[i_hit]);
   ecal_hit_z.push_back(event_in.ecal_hit_z[i_hit]);
   ecal_hit_time.push_back(event_in.ecal_hit_time[i_hit]);
   ecal_hit_index_x.push_back(event_in.ecal_hit_index_x[i_hit]);
   ecal_hit_index_y.push_back(event_in.ecal_hit_index_y[i_hit]);
   return ecal_hit_energy.size()-1;
}

int MiniDst::Add_Ecal_Cluster(MiniDst &event_in, int i_cluster, bool also_copy_ecal_hits) {
   /// Copy Ecal Cluster at i_cluster from event_in and add it to the current event.
   /// This will also add the hits from ecal_hit_ associated with the cluster.
   /// event_in = event to copy from
   /// i_cluster = cluster number in that event.
   ///
   ecal_cluster_energy.push_back(event_in.ecal_cluster_energy[i_cluster]);
   ecal_cluster_time.push_back(event_in.ecal_cluster_time[i_cluster]);
   ecal_cluster_x.push_back(event_in.ecal_cluster_x[i_cluster]);
   ecal_cluster_y.push_back(event_in.ecal_cluster_y[i_cluster]);
   ecal_cluster_z.push_back(event_in.ecal_cluster_z[i_cluster]);
   ecal_cluster_seed_ix.push_back(event_in.ecal_cluster_seed_ix[i_cluster]);
   ecal_cluster_seed_iy.push_back(event_in.ecal_cluster_seed_iy[i_cluster]);
   ecal_cluster_seed_energy.push_back(event_in.ecal_cluster_seed_energy[i_cluster]);

   if(also_copy_ecal_hits) {
      vector<int> copied_cluster_hits;
      bool seed_hit_copied = false;
      for (int i_hit = 0; i_hit < event_in.ecal_cluster_hits[i_cluster].size(); ++i_hit) {
         int new_hit = Add_Ecal_Hit(event_in, event_in.ecal_cluster_hits[i_cluster][i_hit]);
         if (event_in.ecal_cluster_hits[i_cluster][i_hit] == event_in.ecal_cluster_seed_index[i_cluster]) {
            ecal_cluster_seed_index.push_back(new_hit);
            seed_hit_copied = true;
         }
         copied_cluster_hits.push_back(new_hit);
      }
      if(!seed_hit_copied){
         cout << "Seed hit not copied for cluster " << i_cluster << endl;
         ecal_cluster_seed_index.push_back(-1);
      }
      ecal_cluster_hits.push_back(copied_cluster_hits);
      ecal_cluster_nhits.push_back(copied_cluster_hits.size());
   }else{
      ecal_cluster_seed_index.push_back(-1);
      ecal_cluster_hits.push_back({-1});
      ecal_cluster_nhits.push_back(event_in.ecal_cluster_nhits[i_cluster]);
   }

   if(event_in.ecal_cluster_mc_id.size() > 0) {
      ecal_cluster_mc_id.push_back(event_in.ecal_cluster_mc_id[i_cluster]);
      ecal_cluster_mc_pdg.push_back(event_in.ecal_cluster_mc_pdg[i_cluster]);
      ecal_cluster_mc_pdg_purity.push_back(event_in.ecal_cluster_mc_pdg_purity[i_cluster]);
   }

   return ecal_cluster_energy.size()-1;
}

int MiniDst::Add_Track(MiniDst &event_in, int i_track, bool also_copy_svt_hits) {
   /// Copy a track from event_in at i_track and add it to the current event.
   ///
   if(also_copy_svt_hits){
      cout << "Woops, implement the SVT hit copy first!\n";
   }

   if( event_in.track_type[i_track]==0){
      track_n_matched++;
   }else if(event_in.track_type[i_track]<32){
      track_n_kf++;
   }else{
      track_n_gbl++;
   }
   track_type.push_back(event_in.track_type[i_track]);
   track_n_hits.push_back(event_in.track_n_hits[i_track]);
   track_volume.push_back(event_in.track_volume[i_track]);
   track_d0.push_back(event_in.track_d0[i_track]);
   track_phi0.push_back(event_in.track_phi0[i_track]);
   track_omega.push_back(event_in.track_omega[i_track]);
   track_tan_lambda.push_back(event_in.track_tan_lambda[i_track]);
   track_z0.push_back(event_in.track_z0[i_track]);
   track_chi2.push_back(event_in.track_chi2[i_track]);
   track_time.push_back(event_in.track_time[i_track]);
   track_px.push_back(event_in.track_px[i_track]);
   track_py.push_back(event_in.track_py[i_track]);
   track_pz.push_back(event_in.track_pz[i_track]);
   if(event_in.track_x_at_lasthit.size() > i_track) { // The lasthit information is available
      track_x_at_lasthit.push_back(event_in.track_x_at_lasthit[i_track]);
      track_y_at_lasthit.push_back(event_in.track_y_at_lasthit[i_track]);
      track_z_at_lasthit.push_back(event_in.track_z_at_lasthit[i_track]);
      track_omega_at_lasthit.push_back(event_in.track_omega_at_lasthit[i_track]);
      track_tan_lambda_at_lasthit.push_back(event_in.track_tan_lambda_at_lasthit[i_track]);
      track_phi0_at_lasthit.push_back(event_in.track_phi0_at_lasthit[i_track]);
      track_d0_at_lasthit.push_back(event_in.track_d0_at_lasthit[i_track]);
      track_z0_at_lasthit.push_back(event_in.track_z0_at_lasthit[i_track]);
   }
   track_x_at_ecal.push_back(event_in.track_x_at_ecal[i_track]);
   track_y_at_ecal.push_back(event_in.track_y_at_ecal[i_track]);
   track_z_at_ecal.push_back(event_in.track_z_at_ecal[i_track]);
   track_isolation.push_back(event_in.track_isolation[i_track]);  // Places pointer to existing vector.
   track_covmatrix.push_back(event_in.track_covmatrix[i_track]);

   if(event_in.track_lambda_kinks.size() > i_track) {   // KF tracks don't have kinks.
      track_lambda_kinks.push_back(event_in.track_lambda_kinks[i_track]);
      track_phi_kinks.push_back(event_in.track_phi_kinks[i_track]);
   }
   track_particle.push_back(-1);
   track_gbl_ref.push_back(-1);
   track_ref.push_back(-1);
   track_svt_hits.push_back({-1});
   track_hit_pattern.push_back({-2});

   return track_type.size()-1;
}

int MiniDst::Add_Particle(MiniDst &event_in, int i_particle) {
   /// Copy a particle from event_in at i_particle and add it to the current event.
   ///
   int i_ecal = -1;
   if(event_in.part.ecal_cluster[i_particle] >= 0) {
      i_ecal = Add_Ecal_Cluster(event_in, event_in.part.ecal_cluster[i_particle], true);
   }
   int i_track = -1;
   if(event_in.part.track[i_particle] >= 0) {
      i_track = Add_Track(event_in, event_in.part.track[i_particle], false);
   }
   return part.Add(event_in.part, i_particle, i_ecal, i_track);
}

void MiniDst::Print(Option_t *option) const {
   // Pretty print the contents of the current event in the MiniDst class.
   cout << "======== Event: " << event_number << " Run: " << run_number << "  ========\n";
   cout << " Particles: \n";
   for(int i=0; i< part.pdg.size(); ++i){
      cout << "   Particle " << i << ": type: " << part.type[i] << " pdg: " << part.pdg[i] << " energy: " << part.energy[i]
           << " px: " << part.px[i] << " py: " << part.py[i] << " pz: " << part.pz[i]
           << " track: " << part.track[i] << " ecal_cluster: " << part.ecal_cluster[i] << endl;
      if(part.ecal_cluster[i] >= 0) {
         cout << "    |- Ecal Cluster: energy: " << ecal_cluster_energy[part.ecal_cluster[i]]
              << " time: " << ecal_cluster_time[part.ecal_cluster[i]]
              << " x: " << ecal_cluster_x[part.ecal_cluster[i]]
              << " y: " << ecal_cluster_y[part.ecal_cluster[i]]
              << " z: " << ecal_cluster_z[part.ecal_cluster[i]]
              << " seed_ix: " << ecal_cluster_seed_ix[part.ecal_cluster[i]]
              << " seed_iy: " << ecal_cluster_seed_iy[part.ecal_cluster[i]]
              << " seed_energy: " << ecal_cluster_seed_energy[part.ecal_cluster[i]] << endl;
      }
      if(part.track[i] >= 0) {
         cout << "    |- Track: type: " << track_type[part.track[i]]
              << " n_hits: " << track_n_hits[part.track[i]]
              << " d0: " << track_d0[part.track[i]]
              << " phi0: " << track_phi0[part.track[i]]
              << " omega: " << track_omega[part.track[i]]
              << " tan_lambda: " << track_tan_lambda[part.track[i]]
              << " z0: " << track_z0[part.track[i]]
              << " chi2: " << track_chi2[part.track[i]]
              << " time: " << track_time[part.track[i]]
              << " px: " << track_px[part.track[i]]
              << " py: " << track_py[part.track[i]]
              << " pz: " << track_pz[part.track[i]] << endl;
         cout << "     |- At ECal: " << " x: " << track_x_at_ecal[part.track[i]]
              << " y: " << track_y_at_ecal[part.track[i]]
              << " z: " << track_z_at_ecal[part.track[i]] << endl;
      }
   }
}
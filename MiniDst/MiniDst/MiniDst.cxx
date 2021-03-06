//
// Created by Maurik Holtrop on 7/24/20.
//

#include "MiniDst.h"

void MiniDst::Start() {
    DefineBranchMap();
    SetBranchMap();
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

void MiniDst::DefineBranchMap() {
    /// Setup the branch_map and branch_map_active

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
    branch_map_try_emplace("track_n_gbl", &track_n_gbl); /// Number of GBL tracks. The rest are matched tracks.

    if (write_ecal_hits) {
        branch_map_try_emplace("ecal_hit_energy", &ecal_hit_energy);
        branch_map_try_emplace("ecal_hit_time", &ecal_hit_time);
        branch_map_try_emplace("ecal_hit_index_x", &ecal_hit_index_x);
        branch_map_try_emplace("ecal_hit_index_y", &ecal_hit_index_y);

    }

    if (write_ecal_cluster) {
        branch_map_try_emplace("ecal_cluster_energy", &ecal_cluster_energy);
        branch_map_try_emplace("ecal_cluster_time", &ecal_cluster_time);
        branch_map_try_emplace("ecal_cluster_x", &ecal_cluster_x);
        branch_map_try_emplace("ecal_cluster_y", &ecal_cluster_y);
        branch_map_try_emplace("ecal_cluster_z", &ecal_cluster_z);
        branch_map_try_emplace("ecal_cluster_seed_index", &ecal_cluster_seed_index);
        branch_map_try_emplace("ecal_cluster_seed_ix", &ecal_cluster_seed_ix);
        branch_map_try_emplace("ecal_cluster_seed_iy", &ecal_cluster_seed_iy);
        branch_map_try_emplace("ecal_cluster_seed_energy", &ecal_cluster_seed_energy);
        branch_map_try_emplace("ecal_cluster_hits", &ecal_cluster_hits);
        branch_map_try_emplace("ecal_cluster_nhits", &ecal_cluster_nhits);
    }

    if (write_svt_hits) {
        branch_map_try_emplace("svt_hit_layer", &svt_hit_layer);
        branch_map_try_emplace("svt_hit_x", &svt_hit_x);
        branch_map_try_emplace("svt_hit_y", &svt_hit_y);
        branch_map_try_emplace("svt_hit_z", &svt_hit_z);
        branch_map_try_emplace("svt_hit_cxx", &svt_hit_cxx);
        branch_map_try_emplace("svt_hit_cxy", &svt_hit_cxy);
        branch_map_try_emplace("svt_hit_cxz", &svt_hit_cxz);
        branch_map_try_emplace("svt_hit_cyy", &svt_hit_cyy);
        branch_map_try_emplace("svt_hit_cyz", &svt_hit_cyz);
        branch_map_try_emplace("svt_hit_czz", &svt_hit_czz);
        branch_map_try_emplace("svt_hit_time", &svt_hit_time);
        branch_map_try_emplace("svt_hit_edep", &svt_hit_edep);
    }

    if (write_tracks) {
        branch_map_try_emplace("track_n_hits", &track_n_hits); /// The number of 3D hits associated with this track.
        branch_map_try_emplace("track_volume", &track_volume); /// The volume to which this track belongs to.
        branch_map_try_emplace("track_type", &track_type);   /// The track type.

        branch_map_try_emplace("track_d0", &track_d0);   /// The distance of closest approach to the reference point.
        branch_map_try_emplace("track_phi0",
                               &track_phi0);     /// The azimuthal angle of the momentum at the position of closest approach to the reference point.
        branch_map_try_emplace("track_omega",
                               &track_omega); /// The track curvature. The curvature is positive (negative) if the particle has a positive (negative) charge.
        branch_map_try_emplace("track_tan_lambda",
                               &track_tan_lambda); /// The slope of the track in the SY plane where S is the arc length of the helix in the xz plane.
        branch_map_try_emplace("track_z0",
                               &track_z0);      /// The y position of the track at the distance of closest approach in the xz plane.
        branch_map_try_emplace("track_chi2", &track_chi2); /// The chi^2 of the track fit.
        branch_map_try_emplace("track_time",
                               &track_time); /// The time of the track.  This is currently the average time of all hits composing the track.
        branch_map_try_emplace("track_x_at_ecal",
                               &track_x_at_ecal); /// The x position of the extrapolated track at the Ecal face.
        branch_map_try_emplace("track_y_at_ecal",
                               &track_y_at_ecal); /// The y position of the extrapolated track at the Ecal face.
        branch_map_try_emplace("track_z_at_ecal",
                               &track_z_at_ecal); /// The z position of the extrapolated track at the Ecal face.
        branch_map_try_emplace("track_isolation",
                               &track_isolation); /// Array used to store the isolation variables for each of the sensor layers.
        branch_map_try_emplace("track_covmatrix",
                               &track_covmatrix);   /// The 1/2 Covariant Matrix. This is the lower 1/2.
        branch_map_try_emplace("track_lambda_kinks", &track_lambda_kinks);
        branch_map_try_emplace("track_phi_kinks", &track_phi_kinks);
        branch_map_try_emplace("track_particle",
                               &track_particle); /// Reference to the reconstructed particle associated with this track.
        branch_map_try_emplace("track_gbl_ref", &track_gbl_ref);
        branch_map_try_emplace("track_ref", &track_ref);
        branch_map_try_emplace("track_svt_hits",
                               &track_svt_hits);/// Reference to the 3D hits associated with this track.
    }

    branch_map_try_emplace("part_type", &part.type);
    branch_map_try_emplace("part_lcio_type", &part.lcio_type);
    branch_map_try_emplace("part_energy", &part.energy);
    branch_map_try_emplace("part_mass", &part.mass);
    branch_map_try_emplace("part_pdg", &part.pdg);
    branch_map_try_emplace("part_charge", &part.charge);
    branch_map_try_emplace("part_goodness_of_pid", &part.goodness_of_pid);
    branch_map_try_emplace("part_px", &part.px);
    branch_map_try_emplace("part_py", &part.py);
    branch_map_try_emplace("part_pz", &part.pz);
    branch_map_try_emplace("part_track", &part.track);
    branch_map_try_emplace("part_track_chi2", &part.track_chi2);
    branch_map_try_emplace("part_ecal_cluster", &part.ecal_cluster);

    branch_map_try_emplace("v0_type", &v0.type);
    branch_map_try_emplace("v0_lcio_type", &v0.lcio_type);
    branch_map_try_emplace("v0_energy", &v0.energy);
    branch_map_try_emplace("v0_mass", &v0.mass);
    branch_map_try_emplace("v0_pdg", &v0.pdg);
    branch_map_try_emplace("v0_charge", &v0.charge);
    branch_map_try_emplace("v0_goodness_of_pid", &v0.goodness_of_pid);
    branch_map_try_emplace("v0_px", &v0.px);
    branch_map_try_emplace("v0_py", &v0.py);
    branch_map_try_emplace("v0_pz", &v0.pz);
    branch_map_try_emplace("v0_n_daughter", &v0.n_daughter);

    branch_map_try_emplace("v0_vertex_x", &v0.vertex_x);
    branch_map_try_emplace("v0_vertex_y", &v0.vertex_y);
    branch_map_try_emplace("v0_vertex_z", &v0.vertex_z);
    branch_map_try_emplace("v0_vertex_chi2", &v0.vertex_chi2);
    branch_map_try_emplace("v0_vertex_prob", &v0.vertex_prob);

    branch_map_try_emplace("v0_mass_err", &v0.mass_err);


    branch_map_try_emplace("v0_em_part", &v0.em.part);
    branch_map_try_emplace("v0_em_track", &v0.em.track);
    branch_map_try_emplace("v0_em_track_nhit", &v0.em.track_nhit);
    branch_map_try_emplace("v0_em_p", &v0.em.p);
    branch_map_try_emplace("v0_em_chi2", &v0.em.chi2);
    branch_map_try_emplace("v0_em_good_pid", &v0.em.good_pid);
    branch_map_try_emplace("v0_em_track_time", &v0.em.track_time);
    branch_map_try_emplace("v0_em_pos_ecal_x", &v0.em.pos_ecal_x);
    branch_map_try_emplace("v0_em_pos_ecal_y", &v0.em.pos_ecal_y);
    branch_map_try_emplace("v0_em_clus", &v0.em.clus);
    branch_map_try_emplace("v0_em_clus_energy", &v0.em.clus_energy);
    branch_map_try_emplace("v0_em_clus_time", &v0.em.clus_time);
    branch_map_try_emplace("v0_em_clus_ix", &v0.em.clus_ix);
    branch_map_try_emplace("v0_em_clus_iy", &v0.em.clus_iy);
    branch_map_try_emplace("v0_em_clus_pos_x", &v0.em.clus_pos_x);
    branch_map_try_emplace("v0_em_clus_pos_y", &v0.em.clus_pos_y);

    branch_map_try_emplace("v0_ep_part", &v0.ep.part);
    branch_map_try_emplace("v0_ep_track", &v0.ep.track);
    branch_map_try_emplace("v0_ep_track_nhit", &v0.ep.track_nhit);
    branch_map_try_emplace("v0_ep_p", &v0.ep.p);
    branch_map_try_emplace("v0_ep_chi2", &v0.ep.chi2);
    branch_map_try_emplace("v0_ep_good_pid", &v0.ep.good_pid);
    branch_map_try_emplace("v0_ep_track_time", &v0.ep.track_time);
    branch_map_try_emplace("v0_ep_pos_ecal_x", &v0.ep.pos_ecal_x);
    branch_map_try_emplace("v0_ep_pos_ecal_y", &v0.ep.pos_ecal_y);
    branch_map_try_emplace("v0_ep_clus", &v0.ep.clus);
    branch_map_try_emplace("v0_ep_clus_energy", &v0.ep.clus_energy);
    branch_map_try_emplace("v0_ep_clus_time", &v0.ep.clus_time);
    branch_map_try_emplace("v0_ep_clus_ix", &v0.ep.clus_ix);
    branch_map_try_emplace("v0_ep_clus_iy", &v0.ep.clus_iy);
    branch_map_try_emplace("v0_ep_clus_pos_x", &v0.ep.clus_pos_x);
    branch_map_try_emplace("v0_ep_clus_pos_y", &v0.ep.clus_pos_y);


    if (write_mc_particles) {
        // MCParticles
        branch_map_try_emplace("mc_part_energy", &mc_part_energy);
        branch_map_try_emplace("mc_part_pdg_id", &mc_part_pdg_id);
        branch_map_try_emplace("mc_part_gen_status", &mc_part_gen_status); /** Generator Status **/
        branch_map_try_emplace("mc_part_time", &mc_part_time);      /** The global creation time. */
        branch_map_try_emplace("mc_part_x", &mc_part_x);      /** The X vertex. */
        branch_map_try_emplace("mc_part_y", &mc_part_y);      /** The Y vertex. */
        branch_map_try_emplace("mc_part_z", &mc_part_z);      /** The Z vertex. */
        branch_map_try_emplace("mc_part_end_x", &mc_part_end_x);  /** The X end point. */
        branch_map_try_emplace("mc_part_end_y", &mc_part_end_y);  /** The Y end point. */
        branch_map_try_emplace("mc_part_end_z", &mc_part_end_z);  /** The Z end point. */
        branch_map_try_emplace("mc_part_px", &mc_part_px);
        branch_map_try_emplace("mc_part_py", &mc_part_py);
        branch_map_try_emplace("mc_part_pz", &mc_part_pz);
        branch_map_try_emplace("mc_part_mass", &mc_part_mass);
        branch_map_try_emplace("mc_part_charge", &mc_part_charge);
        branch_map_try_emplace("mc_part_daughters", &mc_part_daughters);
        branch_map_try_emplace("mc_part_parents", &mc_part_parents);
    }
}

void MiniDst::SetBranchMap() {
    /// Open the output file, create the output TTree and set the branches in the TTree (if active)

    // Create the output TTree
    if (!md_output_file) {
        if (md_Debug & kDebug_Info) std::cout << "Opening output file: " << md_output_file_name << std::endl;
        md_output_file = new TFile(md_output_file_name.data(), "RECREATE");
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
                [](vector< vector<int> > *arg) { arg->clear(); },
                [](vector< vector<double> >*arg) {arg->clear(); }
        },bran);
    }
}


long MiniDst::Run(int nevt){
    std::cout << "MiniDst::Run() called. =====================<<<<<<<<<<<<<======\n";
    return(0);
}

void MiniDst::End(){
    md_output_file->Write();
    md_output_file->Close();
}
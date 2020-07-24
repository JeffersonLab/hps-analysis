//
// Created by Maurik Holtrop on 7/6/20.
//

#include "Dst2016.h"
#include "TVector3.h"

ClassImp(Dst2016)

Dst2016::Dst2016(TTree *tree, string out_file_name): BaseAna(tree), MiniDst(out_file_name) {
    // SetOutputFileName("MiniDst2016.root");
    fCounter_Freq = 10000;
}

void Dst2016::Start() {
    SlaveBegin();
    MiniDst::Start();
}

long Dst2016::Run(int nevt){
    return(BaseAna::Run(nevt));
}

void Dst2016::End(){
    BaseAna::End();
    MiniDst::End();
}

void Dst2016::SlaveBegin(TTree *tree) {

    // The SlaveBegin() function is called after the Begin() function when processing a chain.
    // When running with PROOF SlaveBegin() is called on each slave server, but Begin is NOT called on each slave,
    // so SlaveBegin() needs to setup the essentials, i.e. Histograms etc.
    // The tree argument is deprecated (on PROOF 0 is passed).

    if (md_Debug & MiniDst::kDebug_L1)cout << "Dst2016::Begin(): \n";
    BaseAna::SlaveBegin(tree);

}

void Dst2016::clear(){
    // Clear the event. We can now just use the map.
    for( auto const& [nam, bran] : branch_map ){
        std::visit([](auto &&arg){arg->clear();},bran);
    }
}

Bool_t  Dst2016::Process(Long64_t entry) {

    clear();
    int stat =GetEntry(entry);
    if(  stat <= 0 ){
        if(md_Debug & MiniDst::kDebug_Error){
            cout << "GetEntry("<< entry << ") returned with status "<<  stat << endl;
            printf("i: %9ld  event: %9d\n", evt_count, event->getEventNumber());
        }
        Abort("Bad event");
        return false;
    }
    if( (evt_count++ % fCounter_Freq ) == 0) {
        printf("i: %9ld  event: %9d  fChain evt: %9lld  output evt: %9lld\n", evt_count, event->getEventNumber(),
                fChain->GetReadEntry(), md_output_tree->GetEntries());
    }

    run_number = GetRunNumber();
    event_number = GetEventNumber();
    trigger = (IsSingle0Trigger() << 3) + (IsSingle1Trigger()<<2) + (IsPair0Trigger() << 1) + (IsPair1Trigger());

//    if( md_Debug & MiniDst::kDebug_L2){
//        cout << "S0: " << IsSingle0Trigger() << " S1: " << IsSingle1Trigger() << "  P0:" << IsPair0Trigger() <<
//        " P1: " << IsPair1Trigger() << endl;
//        cout << "Trigger = " << trigger << endl;
//    }

    // ECAL hits.
    if(write_ecal_hits){
        for(int i=0; i< GetNumberOfEcalHits(); ++i){
            EcalHit *hit = GetEcalHit(i);
            ecal_hit_energy.push_back(hit->energy);
            ecal_hit_time.push_back(hit->hit_time);
            ecal_hit_index_x.push_back(hit->index_x);
            ecal_hit_index_y.push_back(hit->index_y);
        }
    }

    // ECAL Clusters.
    if(write_ecal_cluster) {
        for (int i = 0; i < GetNumberOfEcalClusters(); ++i) {
            EcalCluster *clus = GetEcalCluster(i);
            ecal_cluster_energy.push_back(clus->energy);
            ecal_cluster_time.push_back(clus->cluster_time);
            ecal_cluster_x.push_back(clus->x);
            ecal_cluster_y.push_back(clus->y);
            ecal_cluster_z.push_back(clus->z);
            ecal_cluster_nhits.push_back(clus->n_ecal_hits);

            vector<int> matched_hits;
            auto clus_seed = static_cast<EcalHit *>(clus->seed_hit.GetObject());
            for(int k=0; k < GetNumberOfEcalHits(); ++k){   // For each of the ECAL hits, we know k
                EcalHit *hit = GetEcalHit(k);               // corresponds to the location in our index.
                if(hit == clus_seed ){
                    ecal_cluster_seed_idx.push_back(k);     // Matches the seed hit, so store that index.
                }

                for(int j=0; j < clus->n_ecal_hits; ++j){   // Get each of the hits in the cluster.
                    EcalHit *test_hit = static_cast<EcalHit *>(clus->ecal_hits->At(j));

                    if( test_hit == hit ){                  // So k matches a cluster hit. Store k.
                        matched_hits.push_back(k);
                        break;                              // Already found a match, no need to keep looking.
                    }
                }
            }
            if(matched_hits.size() != clus->n_ecal_hits){   // Make sure we got all of them!
                cout << "WARNING - Did not get all hits for the cluster, need =" << clus->n_ecal_hits <<
                     " got = " << matched_hits.size() << endl;
            }
            ecal_cluster_hits.push_back(matched_hits);      // Store the vector of indexes to hits for this cluster.
        }
    }

    if(write_svt_hits){
        for(int i =0; i< GetNumberOfSvtHits(); ++i){
            SvtHit *svthit = GetSvtHit(i);
            svt_hit_x.push_back(svthit->x);
            svt_hit_y.push_back(svthit->y);
            svt_hit_z.push_back(svthit->z);
            svt_hit_cxx.push_back(svthit->cxx);
            svt_hit_cxy.push_back(svthit->cxy);
            svt_hit_cxz.push_back(svthit->cxz);
            svt_hit_cyy.push_back(svthit->cyy);
            svt_hit_cyz.push_back(svthit->cyz);
            svt_hit_czz.push_back(svthit->czz);
            svt_hit_time.push_back(svthit->time);
            svt_hit_layer.push_back(svthit->layer);

        }
    }

    if(write_tracks){
#ifdef DEBUG
        // Does EVERY GBL track have a valid SEED?
        for(int i=0;i < GetNumberOfGblTracks(); ++i) {
            GblTrack *track = GetGblTrack(i);
            SvtTrack *svt_trk = track->getSeedTrack();
            if( svt_trk == nullptr){
                cout << "Entry " << entry << " has GBL track " << i << " without seed. \n";
            }else{
                GblTrack *gbl_test = svt_trk->getGblTrack();
                if( gbl_test == nullptr){
                    cout << "Entry " << entry << " has GBL track " << i << " with seed that does not point back to GBL. \n";
                } else{
                    if( gbl_test != track){
                        cout << "Entry " << entry << " has GBL track " << i << " with seed, but seed.getGBLTrack does not point back to GBL. \n";
                    }
                }
            }

        }
#endif
        for(int i=0;i < GetNumberOfGblTracks(); ++i) {
            GblTrack *track = GetGblTrack(i);
            vector<double> iso(track->isolation, track->isolation + 14);
            track_isolation.push_back(iso);
            track_n_hits.push_back(track->n_hits);
            track_volume.push_back(track->track_volume);
            track_type.push_back(track->type);
            track_d0.push_back(track->d0);
            track_phi0.push_back(track->phi0);
            track_omega.push_back(track->omega);
            track_tan_lambda.push_back(track->tan_lambda);
            track_z0.push_back(track->z0);
            track_chi2.push_back(track->chi_squared);
            vector<double> covmat(track->covmatrix, track->covmatrix + 15);
            track_covmatrix.push_back(covmat);
            track_time.push_back(track->track_time);
             track_x_at_ecal.push_back(track->x_at_ecal);
            track_y_at_ecal.push_back(track->y_at_ecal);
            track_z_at_ecal.push_back(track->z_at_ecal);

            // Match the svt_hits to the hits on track
            if (write_svt_hits) {
                for (int j = 0; j < track->svt_hits->GetEntries(); ++j) {
                    SvtHit *svt_hit = (SvtHit *) track->svt_hits->At(j);
                    std::vector<int> svt_hit_store;
                    for (int k = 0; k < GetNumberOfSvtHits(); ++k) {
                        SvtHit *test_hit = GetSvtHit(k);
                        if (svt_hit == test_hit) {
                            svt_hit_store.push_back(k);
                            break;
                        }
                    }
                    track_svt_hits.push_back(svt_hit_store);
                }
            }
            // Link of track_particle  --> In particle store.
            // Link of track_gbl_track
            std::vector<double> lambdas(track->lambda_kinks, track->lambda_kinks + 14);
            std::vector<double> phis(track->phi_kinks, track->phi_kinks + 14);
            track_lambda_kinks.push_back(lambdas);
            track_phi_kinks.push_back(phis);
            track_gbl_ref.push_back(i);      // Reference yourself.
            track_ref.push_back(-1);         // Set to invalid. Fill this later when/if processing seed tracks.
        }
    }

    if(write_tracks && !write_only_gbl_tracks){
        for(int i=0;i < GetNumberOfTracks(); ++i) {
            SvtTrack *track = GetTrack(i);
            vector<double> iso(track->isolation, track->isolation + 14);
            track_isolation.push_back(iso);
            track_n_hits.push_back(track->n_hits);
            track_volume.push_back(track->track_volume);
            track_type.push_back(track->type);
            track_d0.push_back(track->d0);
            track_phi0.push_back(track->phi0);
            track_omega.push_back(track->omega);
            track_tan_lambda.push_back(track->tan_lambda);
            track_z0.push_back(track->z0);
            track_chi2.push_back(track->chi_squared);
            vector<double> covmat(track->covmatrix, track->covmatrix + 15);
            track_covmatrix.push_back(covmat);
            track_time.push_back(track->track_time);
            track_x_at_ecal.push_back(track->x_at_ecal);
            track_y_at_ecal.push_back(track->y_at_ecal);
            track_z_at_ecal.push_back(track->z_at_ecal);

            // Match the svt_hits to the hits on track
            if (write_svt_hits) {
                for (int j = 0; j < track->svt_hits->GetEntries(); ++j) {
                    SvtHit *svt_hit = (SvtHit *) track->svt_hits->At(j);
                    std::vector<int> svt_hit_store;
                    for (int k = 0; k < GetNumberOfSvtHits(); ++k) {
                        SvtHit *test_hit = GetSvtHit(k);
                        if (svt_hit == test_hit) {
                            svt_hit_store.push_back(k);
                            break;
                        }
                    }
                    track_svt_hits.push_back(svt_hit_store);
                }
            }
            track_ref.push_back(track_ref.size());  // References self.
            track_gbl_ref.push_back(-1);            // Fill these when running over the GBL tracks again.
        }

        // Loop over GBL tracks to link GBL -> Seed track.
        for(int i=0;i < GetNumberOfGblTracks(); ++i) {
            GblTrack *track = GetGblTrack(i);
            SvtTrack *svt_trk = track->getSeedTrack();
            if( svt_trk == nullptr){
                cout << "Entry " << entry << " has GBL track " << i << " without seed. \n";
            }else{
#ifdef DEBUG
                GblTrack *gbl_test = svt_trk->getGblTrack();
                if( gbl_test == nullptr){
                    cout << "Entry " << entry << " has GBL track " << i << " with seed that does not point back to GBL. \n";
                } else{
                    if( gbl_test != track){
                        cout << "Entry " << entry << " has GBL track " << i << " with seed, but seed.getGBLTrack does not point back to GBL. \n";
                    }
                }
#endif
                bool found = false;
                for(int j=0; j< GetNumberOfTracks(); ++j){
                    SvtTrack *test_track = GetTrack(j);
                    if( test_track == svt_trk){ // We found the seed belong to the GBL track.
                        track_ref[i] = j;
                        track_gbl_ref[j] = i;
                        found = true;
                        break;
                    }
                }
#ifdef DEBUG
                if(!found){
                    cout << "Problem finding the SEED track in the SvtTrack collection. \n";
                }
#endif
            }
        }
    }

    // Single Track Particles.
    for( auto type: particle_types_single) {
        for (int i = 0; i < GetNumberOfParticles(HpsParticle::ParticleType(type)); ++i) {
            HpsParticle *part = GetParticle(HpsParticle::ParticleType(type), i);
            part_type.push_back( int(type));
            part_energy.push_back(part->energy);
            part_pdg.push_back(part->pdg);
            part_charge.push_back(part->charge);
            part_goodness_of_pid.push_back(part->getGoodnessOfPID());
            part_px.push_back(part->px);
            part_py.push_back(part->py);
            part_pz.push_back(part->pz);
            part_corr_px.push_back(part->px_corr);
            part_corr_py.push_back(part->py_corr);
            part_corr_pz.push_back(part->pz_corr);
            part_vertex_x.push_back(part->vtx_x);
            part_vertex_y.push_back(part->vtx_y);
            part_vertex_z.push_back(part->vtx_z);
            part_vertex_chi2.push_back(part->vtx_fit_chi2);
            // Find and linkup the track associated with the particle. == Usually just one?
#ifdef DEBUG
            if(part->svt_tracks->GetEntries() > 1){
                std::cout<< "Particle with more than one track: " << part->svt_tracks->GetEntries() << endl;
            }
#endif
            int part_track_idx = -10;
            if(part->svt_tracks->GetEntries() == 1) {
                GblTrack *gbl_track = (GblTrack *)part->svt_tracks->At(0);

                for (int k = 0; k < GetNumberOfGblTracks(); k++) {
                    GblTrack *test_track = GetGblTrack(k);
                    if (test_track == gbl_track) {
                        part_track_idx = k;
                        break;
                    }
                }
            }
            part_track.push_back(part_track_idx); // No track for this particle = neutral.


            // Find and linkup the Ecal clusters. Only one!
#ifdef DEBUG
            if(part->ecal_clusters->GetEntries() > 1){
                std::cout<< "Particle with more than one ecal cluster! " << part->ecal_clusters->GetEntries() << endl;
            }
#endif
            if(part->ecal_clusters->GetEntries() ==1){
                EcalCluster *p_clus = (EcalCluster *) part->ecal_clusters->At(0);
                for(int k = 0; k < GetNumberOfEcalClusters(); ++k){
                    EcalCluster *test_clus = GetEcalCluster(k);
                    if( test_clus == p_clus){
                        part_ecal_cluster.push_back(k);
                        break;
                    }
                }
            }else{
                part_ecal_cluster.push_back(-1);
            }
        }
    }

    // Double Track Particles, i.e.  Constrained V0 Candidates.
    for(int type: particle_types_double) {
        for (int i = 0; i < GetNumberOfParticles(HpsParticle::ParticleType(type)); ++i) {
            HpsParticle *part = GetParticle(HpsParticle::ParticleType(type), i);
            v0_type.push_back( int(type));
            v0_energy.push_back(part->energy);
            v0_mass.push_back(part->getMass());
            v0_pdg.push_back(part->pdg);
            v0_charge.push_back(part->charge);
            v0_goodness_of_pid.push_back(part->getGoodnessOfPID());
            v0_px.push_back(part->px);
            v0_py.push_back(part->py);
            v0_pz.push_back(part->pz);
            v0_corr_px.push_back(part->px_corr);
            v0_corr_py.push_back(part->py_corr);
            v0_corr_pz.push_back(part->pz_corr);
            v0_vertex_x.push_back(part->vtx_x);
            v0_vertex_y.push_back(part->vtx_y);
            v0_vertex_z.push_back(part->vtx_z);
            v0_vertex_chi2.push_back(part->vtx_fit_chi2);
            v0_n_daughter.push_back(part->n_daughters);

            // We only look for and store 2 daughter particles.
#ifdef DEBUG
            if (part->n_daughters != 2)
                std::cout << "Weird, but I expected 2 and only 2 daughters, but got " << part->n_daughters << std::endl;
#endif
            /// The following is for convenience in analysis.
            /// Here we identify the electron and positron that make the V0 particle, and fill
            /// some useful quantities for each:  chi2 of track, goodness_of_pid, track time, cluster time.
            //
            // Todo: This code can be prettified: the em and ep code duplicates and could done with a struct instead.
            if (part->n_daughters > 1) {
                auto *d0 = dynamic_cast<HpsParticle *>(part->getParticles()->At(0)); // getParticles returns TRefArray.
                auto *d1 = dynamic_cast<HpsParticle *>(part->getParticles()->At(1));
                int found_idx_em = -99;
                int found_idx_ep = -99;
                int found_ep_track = -99;
                int found_em_track = -99;
                int found_em_track_nhit=-99;
                int found_ep_track_nhit=-99;

                double found_em_mom = -999.;
                double found_ep_mom = -999.;
                double found_em_chi2 = -999.;
                double found_ep_chi2 = -999.;
                double found_em_good_pid= -999.;
                double found_ep_good_pid= -999.;
                double found_em_track_time = -999.;
                double found_ep_track_time = -999.;
                double found_em_clus_time = -999.;
                double found_ep_clus_time = -999.;
                double found_em_pos_ecal_x = -999.;
                double found_em_pos_ecal_y = -999.;
                double found_ep_pos_ecal_x = -999.;
                double found_ep_pos_ecal_y = -999.;

                for (int idx = 0; idx < GetNumberOfParticles(HpsParticle::FINAL_STATE_PARTICLE); ++idx) {
                    HpsParticle *part = GetParticle(HpsParticle::FINAL_STATE_PARTICLE, idx);
                    if (d0 == part || d1 == part) {
                        if( part->getPDG() == 11) {
                            found_idx_em = idx;
                            TRefArray *track_refs = part->getTracks();
                            if(track_refs->GetEntries() == 0){  // no track?
                                found_em_chi2 = -1;
                            }else {                             // Should be only one track.
                                GblTrack *track = (GblTrack *) track_refs->At(0);
                                found_em_track = part_track[idx];
#ifdef DEBUG
                                GblTrack *test_track = GetGblTrack(found_em_track);
                                if( test_track != track){
                                    cout << "We have a electron particle -> track pointer mismatch. \n";
                                }
#endif
                                found_em_chi2 = track->getChi2();
                                found_em_track_nhit = track->getSvtHits()->GetEntries();
                                found_em_track_time = track->getTrackTime();
                                vector<double> mom = track->getMomentum();
                                found_em_mom = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
                                vector<double> pos_ecal = track->getPositionAtEcal();
                                found_em_pos_ecal_x = pos_ecal[0];
                                found_em_pos_ecal_y = pos_ecal[1];

                            }

                            TRefArray *cluster_refs = part->getClusters();
                            if( cluster_refs->GetEntries() == 0){
                                cout << "Weird electron without a cluster-ref!\n";
                            }else{
#ifdef DEBUG
                                if(cluster_refs->GetEntries()>1){
                                    cout << "Weird electron with more than one cluster!! \n";
                                }
#endif
                                EcalCluster *clust = (EcalCluster *)cluster_refs->At(0) ;
                                found_em_clus_time = clust->getClusterTime();
                            }
                            found_em_good_pid = part->getGoodnessOfPID();

                        }else if( part->getPDG() == -11){
                            found_idx_ep = idx;
                            TRefArray *track_refs = part->getTracks();
                            if(track_refs->GetEntries() == 0){
                                found_ep_chi2 = -1;
                            }else {
                                GblTrack *track = (GblTrack *) track_refs->At(0);
                                found_ep_track = part_track[idx];
#ifdef DEBUG
                                GblTrack *test_track = GetGblTrack(found_ep_track);
                                if( test_track != track){
                                    cout << "We have a positron particle -> track pointer mismatch. \n";
                                }
#endif
                                found_ep_chi2 = track->getChi2();
                                found_ep_track_nhit = track->getSvtHits()->GetEntries();
                                found_ep_track_time = track->getTrackTime();
                                vector<double> mom = track->getMomentum();
                                found_ep_mom = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
                                vector<double> pos_ecal = track->getPositionAtEcal();
                                found_ep_pos_ecal_x = pos_ecal[0];
                                found_ep_pos_ecal_y = pos_ecal[1];
                            }
                            TRefArray *cluster_refs = part->getClusters();
                            if( cluster_refs->GetEntries() == 0){
                                cout << "Weird positron without a cluster-ref!\n";
                            }else{
#ifdef DEBUG
                                if(cluster_refs->GetEntries()>1){
                                    cout << "Weird positron with more than one cluster!! \n";
                                }
#endif
                                EcalCluster *clust = (EcalCluster *)cluster_refs->At(0) ;
                                found_ep_clus_time = clust->getClusterTime();
                            }

                            found_ep_good_pid = part->getGoodnessOfPID();
                        }else{
                            std::cout << "Problem with particle, neither e- nor e+ \n";
                        }
                    }
                    if (found_idx_em > 0 && found_idx_ep > 0) break;
                }
#ifdef DEBUG
                if( found_idx_em == found_idx_ep || found_idx_em <0 || found_idx_ep <0 ){
                    std::cout << "Problem with particle. The e+ or e- was not found, or e+ == e- \n";
                }
                if( part_charge[found_idx_em]>0 || part_charge[found_idx_ep]<0 ){
                    std::cout << "Problem with particle: e- has +charge or e+ has -charge.\n";
                }
                if( track_chi2[part_track[found_idx_em]] != found_em_chi2){
                    std::cout << "We got the electron chi2 wrong: " << track_chi2[part_track[found_idx_em]] << " != " << found_em_chi2 << endl;
                }
                if( track_chi2[part_track[found_idx_ep]] != found_ep_chi2){
                    std::cout << "We got the positron chi2 wrong: " << track_chi2[part_track[found_idx_ep]] << " != " << found_ep_chi2 << endl;
                }
#endif
                v0_em_part.push_back(found_idx_em);
                v0_ep_part.push_back(found_idx_ep);
                v0_em_track.push_back(found_em_track);
                v0_ep_track.push_back(found_ep_track);
                v0_em_track_nhit.push_back(found_em_track_nhit);
                v0_ep_track_nhit.push_back(found_em_track_nhit);
                v0_em_p.push_back(found_em_mom);
                v0_ep_p.push_back(found_ep_mom);

                v0_em_chi2.push_back(found_em_chi2);
                v0_ep_chi2.push_back(found_ep_chi2);
                v0_em_good_pid.push_back(found_em_good_pid);
                v0_ep_good_pid.push_back(found_ep_good_pid);
                v0_em_track_time.push_back(found_em_track_time);
                v0_ep_track_time.push_back(found_ep_track_time);
                v0_em_clus_time.push_back(found_em_clus_time);
                v0_ep_clus_time.push_back(found_ep_clus_time);
                v0_em_pos_ecal_x.push_back(found_em_pos_ecal_x);
                v0_em_pos_ecal_y.push_back(found_em_pos_ecal_y);
                v0_ep_pos_ecal_x.push_back(found_ep_pos_ecal_x);
                v0_ep_pos_ecal_y.push_back(found_ep_pos_ecal_y);
            }

            // Find and linkup the track associated with the particle.
//            vector<int> track_list;
//            for(int j=0; j< part->svt_tracks->GetEntries(); ++j) {
//                SvtTrack *p_track = (SvtTrack *)part->svt_tracks->At(j);
//                for (int k = 0; k < GetNumberOfTracks(); k++) {
//                    SvtTrack *test_track = GetTrack(k);
//                    if (test_track == p_track) {
//                        track_list.push_back(k);
//                        break;
//                    }
//                }
//            }
//            v0_tracks.push_back(track_list);
            // Find and linkup the Ecal clusters. == Usually just one?
//            vector<int> ecal_list;
//            for(int j=0; j< part->ecal_clusters->GetEntries(); ++j){
//                EcalCluster *p_clus = (EcalCluster *) part->ecal_clusters->At(j);
//                for(int k = 0; k < GetNumberOfEcalClusters(); ++k){
//                    EcalCluster *test_clus = GetEcalCluster(k);
//                    if( test_clus == p_clus){
//                        ecal_list.push_back(k);
//                        break;
//                    }
//                }
//            }
//            v0_ecal_clusters.push_back(ecal_list);
        }
    }

    if(write_mc_particles) {
        for (int i = 0; i < GetNumberOfMCParticles(); ++i) {
            MCParticle *mc_part = GetMCParticle(i);
            mc_part_energy.push_back(mc_part->energy_);
            mc_part_pdg_id.push_back(mc_part->pdg_id_);
            mc_part_gen_status.push_back(mc_part->gen_status_);
            mc_part_time.push_back(mc_part->time_);
            mc_part_x.push_back(mc_part->x_);
            mc_part_y.push_back(mc_part->y_);
            mc_part_z.push_back(mc_part->z_);
            mc_part_end_x.push_back(mc_part->end_x_);
            mc_part_end_y.push_back(mc_part->end_y_);
            mc_part_end_z.push_back(mc_part->end_z_);
            mc_part_px.push_back(mc_part->px_);
            mc_part_py.push_back(mc_part->py_);
            mc_part_pz.push_back(mc_part->pz_);
            mc_part_mass.push_back(mc_part->mass_);
            mc_part_charge.push_back(mc_part->charge_);
            vector<int> daughters;
            vector<int> parents;
            for (int j = 0; j < GetNumberOfMCParticles(); ++j) {
                MCParticle *test_part = GetMCParticle(j);
                for (int k = 0; k < mc_part->daughters_->GetEntries(); ++k) {
                    MCParticle *daught = (MCParticle *) mc_part->daughters_->At(k);
                    if (daught == test_part) {
                        daughters.push_back(j);
                        break;
                    }
                }
                for (int k = 0; k < mc_part->parents_->GetEntries(); ++k) {
                    MCParticle *paren = (MCParticle *) mc_part->parents_->At(k);
                    if (paren == test_part) {
                        parents.push_back(j);
                        break;
                    }
                }
            }
            mc_part_daughters.push_back(daughters);
            mc_part_parents.push_back(parents);
        }
    }
    md_output_tree->Fill();
    return true;
};

void Dst2016::Terminate() {
    if(md_Debug & MiniDst::kDebug_L1) cout << "Dst2016::Terminate() \n";

    TList *list=GetOutputList();
    if(md_Debug & MiniDst::kDebug_L1) cout << "The list has "<< list->GetEntries() << " entries \n";

    // Write the contends to a file/
    if(md_Debug & MiniDst::kDebug_L1) cout << "Writing output file: " << md_output_file_name << endl;
    if( !md_output_file) md_output_file = new TFile(md_output_file_name.data(),"RECREATE");

    if(!md_output_tree){
        md_output_tree = dynamic_cast<TTree *>(list->FindObject("MiniDst"));
    }

    WriteList(list);

    md_output_file->Write();
    md_output_file->Close();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Graveyard of dead code.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  From the loop over the SvtTracks, finding the GBL Track. Now instead, we loop GblTracks.
//            int gbl_track_index = -1;
//            if (track->gbl_track.IsValid()) {
//                // Link track to GBL track.
//                GblTrack *gbl_track_ref = track->getGblTrack();
//                SvtTrack *gbl_track_ref_seed = gbl_track_ref->getSeedTrack();
//                // printf("For entry %4d with i=%2d  Ngbl = %2d Nsvt = %2d \n",entry,i, GetNumberOfGblTracks(), GetNumberOfTracks());
//                if( gbl_track_ref_seed != NULL) { // The root files are not written properly with a call to BrachRef().
//                    if (track != gbl_track_ref_seed) {
//                        cout << "Track and seed pointers are inconsistent. \n";
//                    } else if (gbl_track_ref_seed->getD0() != track->getD0()) {
//                        cout << "Inconsistence in GBL track seeed: delta d0="
//                             << gbl_track_ref_seed->getD0() - track->getD0() << endl;
//                    }
//                }
//                for(int j=0; j< GetNumberOfGblTracks(); ++j){
//                    GblTrack *test_gbl = GetGblTrack(j);
//                    SvtTrack *seed = test_gbl->getSeedTrack();
//                    if( seed == nullptr){
//                        cout << "ERROR -- We have a null SEED pointer!!!" << endl;
//                    }
//                    if( seed == track){ // This is the track that is the seed to the GBL Track. Store the pointers.
//                        if( gbl_track_ref_seed != NULL && gbl_track_ref_seed != seed ){
//                            cout << "The gbl_track_ref_seed was pointing incorrectly! ============ \n";
//                        }
//                        track_ref[j] = track_gbl_ref.size();
//                        gbl_track_index = j;
//                        break;
//                    }
//                }
//            }else{ // Seed that is marked to not have a valid GBL track so we skip it.
//                // printf("For entry %4d with i=%2d  No GBL. \n",entry,i, GetNumberOfGblTracks(), GetNumberOfTracks());
//                gbl_track_index = -2;
//            }
//            track_gbl_ref.push_back(gbl_track_index);

//
// Created by Maurik Holtrop on 8/10/20.
//

#include "GammaMixer.h"

void GammaMixer::Start() {
    /// Setup the mixer.
    DefineBranchMap();

    auto br_names = GetBranchNames();

    // Set the branches not starting with ecal_ or part_ to false.
    std::for_each(br_names.begin(), br_names.end(), [this](string n){
        // If none of the strings are found, set the branch to false.
        if( n.compare(0,5,"part_")!=0 && n.compare(0,5,"ecal_")!=0 ) {
            if(md_Debug & kDebug_L1) cout << "Setting branch " << n << " to false\n";
            this->SetBranchActive(n, false);
        }
    });
    // Add a few back in (easier than making exceptions for them in lambda.
    SetBranchActive("run_number", true);
    SetBranchActive("event_number", true);
    SetBranchActive("trigger", true);
    SetBranchActive("rf_time1", true);
    SetBranchActive("rf_time2", true);
    SetBranchActive("time_stamp", true);

    auto active_branch_names = GetActiveBranchNames(); // These are the standard active branches without the extra.
    if( md_Debug & kDebug_L1) {
        cout << "Active branches: [";
        std::for_each(active_branch_names.begin(), active_branch_names.end(), [](string n){
            cout << n << ",";
        });
        cout << "]\n";
    }
    // Add one more branch for event_number2
    branch_map_try_emplace("event_number2", &event_number2);
    SetBranchMap();  // Connect the branches to the standard output tree.

    event1.DefineBranchMap();
    event2.DefineBranchMap();
    std::for_each(br_names.begin(),br_names.end(),[this,&active_branch_names](string name){
        // if we do not find the name in the active branch name list, set that branch to false.
        if( std::find(active_branch_names.begin(),active_branch_names.end(),name) == active_branch_names.end()) {
            if( md_Debug & kDebug_L1) {cout << "GM: Branch " << name << " set to false\n"; }
            this->event1.SetBranchActive(name, false);
            this->event2.SetBranchActive(name, false);
        }
    });
}



long GammaMixer::Run(int max_event) {
    /// Run through max_event input events and mix them.

    for( auto filename: input_file_list){
        TFile infile1(filename.c_str());
        tree1 = infile1.Get<TTree>("MiniDST");
        event1.SetBranchAddressesOnTree(tree1);  // Connect event1 to the tree1 branches.
        TFile infile2(filename.c_str());    // Open file *again* to get an independent TTree.
        tree2 = infile2.Get<TTree>("MiniDST");
        event2.SetBranchAddressesOnTree(tree2);  // Connect event2 to the tree2 branches.

        long n_evt_out = 0;
        long n_evt_in1 = 0;
        long n_evt_in2 = 0;
        long max_entries = (long)tree1->GetEntries();
        for(long i_evt1 =0; i_evt1 < max_entries; ++i_evt1){
            tree1->GetEntry(i_evt1);
            n_evt_in1++;
            // Check if the event1 has a good photon pair.
            int e1_gamma_1 =0;
            int e1_gamma_2 = 1;
            bool gamma1_ok;
            int n_evt_for_current_evt1 = 0;
            while(Find_Good_Photon_Pair(event1, e1_gamma_1, e1_gamma_2 ) &&
                    n_evt_for_current_evt1 < mix_multiplyer){
                   for(long i_evt2 = i_evt1+1; i_evt2 < max_entries && n_evt_for_current_evt1 < mix_multiplyer; ++i_evt2){
                       tree2->GetEntry(i_evt2);
                       n_evt_in2++;
                       if( (n_evt_in1 + n_evt_in2)%Counter_Freq == 0){
                           printf("i1: %'10ld  i2: %'10ld  event1: %'8d  event2: %'8d  run: %5d\n",
                                  i_evt1, i_evt2, event1.event_number, event2.event_number, event1.run_number);
                       }
                       int e2_gamma_1 =0;
                       int e2_gamma_2 = 1;
                       bool gamma2_ok;
                       while(Find_Good_Photon_Pair(event2, e2_gamma_1, e2_gamma_2)){
                           ///////////////////////////////////// Event1 - Event2 Photon Pair matching conditions /////////////
                           ///
                           /// Currently none.
                           ///
                           //////////////////////////////////////////////////////////////////////////////////////////////////
                           Write_Mixed_Photon_Events(event1, e1_gamma_1, event2, e2_gamma_2);
                           n_evt_out++;
                           if(max_event>0 && n_evt_out >= max_event) return(n_evt_out);
                           e2_gamma_1++;
                           e2_gamma_2 = e2_gamma_1 +1;
                           n_evt_for_current_evt1++;
                       }

                   }
                   e1_gamma_1++;
                   e1_gamma_2 = e1_gamma_1 +1;
            }
        }
    }

    return(0);
}

bool GammaMixer::Find_Good_Photon_Pair(MiniDst &event, int &found1, int &found2) {
    /// Starting from found1 and found2, search for a pair of good photons.
    /// If found, return true and place the part indexes in found1 and found2.
    if( found1 >= found2 ) return false;

    for(int i_gamma1=found1; i_gamma1< event.part.pdg.size(); ++i_gamma1){
        if( event.part.pdg[i_gamma1] != 22) continue; // Not a photon.
        int i_clus1 = event.part.ecal_cluster[i_gamma1];
        if( i_clus1 < 0 || i_clus1 >= event.ecal_cluster_time.size() ) continue; // Bad cluster.

        for(int i_gamma2=found2; i_gamma2 < event.part.pdg.size(); ++i_gamma2){
            if( i_gamma2 <= i_gamma1 ) continue;
            if( event.part.pdg[i_gamma2] != 22) continue; // Not a photon.
            int i_clus2 = event.part.ecal_cluster[i_gamma2];
            if( i_clus2 < 0 || i_clus2 >= event.ecal_cluster_time.size() ) continue; // Bad cluster.

            //////////////////////////////// GOOD PHOTON PAIR CONDITIONS /////////////////////////////////
            ///
            ///  Delta-t < 2. ns.
            ///  Esum < Ebeam*1.1;
            ///
            /////////////////////////////////////////////////////////////////////////////////////////////
            double clus_time1 = event.ecal_cluster_time[i_clus1];
            double clus_time2 = event.ecal_cluster_time[i_clus2];
            if( abs(clus_time1 - clus_time2) > 2.) continue;  // Too far apart in time.

            double clus_energy1 = event.ecal_cluster_energy[i_clus1];
            double clus_energy2 = event.ecal_cluster_energy[i_clus2];
            if( (clus_energy1 + clus_energy2) > ecal_cluster_sum_max) continue;

            // Found a pair, so return.
            found1 = i_gamma1;
            found2 = i_gamma2;
            return(true);
        }
    }


    return false;
}

bool GammaMixer::Write_Mixed_Photon_Events(MiniDst &event1, int e1_gamma,
                                           MiniDst &event2, int e2_gamma) {
/// Write out an event to the output TTree with two gamma's in it, one from event1 and the other from event2.
/// Also populate the run_number, event_numbers etc.

// ToDo: Consider a more elegant event copy scheme? Could use structs throughout in MiniDst.h
    run_number = event1.run_number;
    event_number  = event1.event_number;
    event_number2 = event2.event_number;
    trigger   = event1.trigger;
    time_stamp = event1.time_stamp;
    rf_time1 = event1.rf_time1;
    rf_time2 = event2.rf_time1;  // This is a bit arbitrary. Possibly need something sensible here later.

    // Event1 gamma
    int e1_clus = event1.part.ecal_cluster[e1_gamma];
    ecal_cluster_energy.push_back(event1.ecal_cluster_energy[e1_clus]);
    ecal_cluster_time.push_back(event1.ecal_cluster_time[e1_clus]);
    ecal_cluster_x.push_back(event1.ecal_cluster_x[e1_clus]);
    ecal_cluster_y.push_back(event1.ecal_cluster_y[e1_clus]);
    ecal_cluster_z.push_back(event1.ecal_cluster_z[e1_clus]);
    ecal_cluster_seed_energy.push_back(event1.ecal_cluster_seed_energy[e1_clus]);

    vector<int> new_cluster_hits;
    for( auto hit: event1.ecal_cluster_hits[e1_clus]){
        ecal_hit_energy.push_back(event1.ecal_hit_energy[hit]);
        ecal_hit_time.push_back(event1.ecal_hit_time[hit]);
        ecal_hit_index_x.push_back(event1.ecal_hit_index_x[hit]);
        ecal_hit_index_y.push_back(event1.ecal_hit_index_y[hit]);
        new_cluster_hits.push_back(ecal_hit_energy.size()-1);
    }
    ecal_cluster_hits.push_back(new_cluster_hits);

    int new_cluster_id = ecal_cluster_energy.size()-1;
    part.Add(event1.part, e1_gamma,new_cluster_id,-1);

    // Event2 gamma
    int e2_clus = event2.part.ecal_cluster[e2_gamma];
    ecal_cluster_energy.push_back(event2.ecal_cluster_energy[e2_clus]);
    ecal_cluster_time.push_back(event2.ecal_cluster_time[e2_clus]);
    ecal_cluster_x.push_back(event2.ecal_cluster_x[e2_clus]);
    ecal_cluster_y.push_back(event2.ecal_cluster_y[e2_clus]);
    ecal_cluster_z.push_back(event2.ecal_cluster_z[e2_clus]);
    ecal_cluster_seed_energy.push_back(event2.ecal_cluster_seed_energy[e2_clus]);

    new_cluster_hits.clear();
    for( auto hit: event2.ecal_cluster_hits[e2_clus]){
        ecal_hit_energy.push_back(event2.ecal_hit_energy[hit]);
        ecal_hit_time.push_back(event2.ecal_hit_time[hit]);
        ecal_hit_index_x.push_back(event2.ecal_hit_index_x[hit]);
        ecal_hit_index_y.push_back(event2.ecal_hit_index_y[hit]);
        new_cluster_hits.push_back(ecal_hit_energy.size()-1);
    }
    ecal_cluster_hits.push_back(new_cluster_hits);

    new_cluster_id = ecal_cluster_energy.size()-1;
    part.Add(event2.part, e2_gamma,new_cluster_id,-2);

    return true;
}

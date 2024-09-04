//
// Created by Maurik Holtrop on 8/10/20.
//

#include "MollerMixer.h"

void MollerMixer::Start() {
   /// Setup the mixer.
   DefineBranchMap();

   auto br_names = GetBranchNames();

   // Set the branches not starting with ecal_ or part_ or track_ to false.
   std::for_each(br_names.begin(), br_names.end(), [this](string n){
      // If none of the strings are found, set the branch to false.
      if(n.compare(0,5,"part_")!=0 && n.compare(0,5,"ecal_")!=0 && n.compare(0,6,"track_")!=0 ) {
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
//   branch_map_try_emplace("two_electron_esum_original1", &two_electron_esum_original1);
//   branch_map_try_emplace("two_electron_esum_original2", &two_electron_esum_original2);
//   branch_map_try_emplace("two_electron_esum_mixed", &two_electron_esum_mixed);
//   branch_map_try_emplace("two_electron_mass_original1", &two_electron_mass_original1);
//   branch_map_try_emplace("two_electron_mass_original2", &two_electron_mass_original2);
   branch_map_try_emplace("two_electron_mass_mixed", &two_electron_mass_mixed);
   branch_map_try_emplace("n_evt1", &n_evt1);
   branch_map_try_emplace("n_evt2", &n_evt2);

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



long MollerMixer::Run(int max_event) {
   /// Run through max_event input events and mix them.

   long n_evt_out = 0;
   long n_evt_in1 = 0;
   long n_evt_in2 = 0;

   if(input_file_list.size()<2){
      input_file_list.push_back(input_file_list[0]); // Just duplicate same file.
   }

   for(int i_file=0; i_file < input_file_list.size()-1; ++i_file) {
      cout << "Opening file " << input_file_list[i_file] << " and " << input_file_list[i_file+1] << " for reading. \n";
      TFile infile1(input_file_list[i_file].c_str());
      tree1 = infile1.Get<TTree>("MiniDST");
      event1.SetBranchAddressesOnTree(tree1);  // Connect event1 to the tree1 branches.
      TFile infile2(input_file_list[i_file+1].c_str());
      tree2 = infile2.Get<TTree>("MiniDST");
      event2.SetBranchAddressesOnTree(tree2);  // Connect event2 to the tree2 branches.

      long max_entries1 = (long) tree1->GetEntries();
      long max_entries2 = (long) tree2->GetEntries();
      for (long i_evt1 = 0; i_evt1 < max_entries1; ++i_evt1) {
         tree1->GetEntry(i_evt1);
         n_evt_in1++;
         n_evt_in2 = 0;
         int number_of_evt1 = 0;  // Count how many times i_evt1 is re-used. Should be ~ mix_multiplyer.
         // Note: Cannot keep count in n_evt1 because it is set to zero in Clear().

         for (long i_evt2 = i_evt1+1; i_evt2 < max_entries2 && number_of_evt1 < mix_multiplyer; ++i_evt2) {
            tree2->GetEntry(i_evt2);

            n_evt_in2++;
            int e1_el_1 = 0;
            int e2_el_2 = 0;
            int number_of_evt2 = 0; // Count how many times i_evt2 is re-used.
            while (Find_Good_Moller_Pair(event1, event2, e1_el_1, e2_el_2)) {
               ///////////////////////////////////// Event1 - Event2 Moller Pair /////////////

               ROOT::Math::PxPyPzMVector p1(event1.part.px[e1_el_1], event1.part.py[e1_el_1], event1.part.pz[e1_el_1],
                                            0.000510998949996);
               ROOT::Math::PxPyPzMVector p2(event2.part.px[e2_el_2], event2.part.py[e2_el_2], event2.part.pz[e2_el_2],
                                            0.000510998949996);

               // Clear out the last output event.
               Clear();

               n_evt1 = number_of_evt1;
               n_evt2 = number_of_evt2;

               // Convenience: Compute the two photon invariant masses.
               auto psum = p1 + p2;
               two_electron_mass_mixed = psum.mass();

              Write_Mixed_Moller_Events(event1, e1_el_1, event2, e2_el_2);

               number_of_evt1++;
               number_of_evt2++;

               // Fill the output TTree
               md_output_tree->Fill();

               // Update output counter and inform user.
               n_evt_out++;

               if (n_evt_out % Counter_Freq == 0) {
                  printf("out: %'10ld  n1: %'8ld  n2: %'8ld  event1: %'8d.%02d  event2: %'8d.%02d\n",
                         n_evt_out, n_evt_in1, n_evt_in2, event1.event_number, e1_el_1, event2.event_number, e2_el_2);
               }

               // Enough already?
               if (max_event > 0 && n_evt_out >= max_event) return (n_evt_out);

               // Move on to the next electron(s).
               e2_el_2++;
               if(e2_el_2 >= event2.part.pdg.size()) {
                  e2_el_2 = 0;
                  e1_el_1++;
               }
            }
         }
      }
   }
return(0);
}

bool MollerMixer::Find_Good_Moller_Pair(MiniDst &event1, MiniDst &event2, int &found1, int &found2) {
   /// Starting from found1 and found2, search for a pair of good electrons^*.
   /// If found, return true and place the part indexes in found1 and found2.
   /// * A pair of good electrons here is simply PDG==11 and P_sum_min < Psum < P_sum_max, and each electron P< p_max.
   int i_electron1=found1;
   int i_electron2 = found2;
   while( i_electron1< event1.part.pdg.size() && i_electron2 < event2.part.pdg.size()) {
      if (event1.part.pdg[i_electron1] == 11 && event2.part.pdg[i_electron2] == 11) { // 2 electrons.
         //////////////////////////////// GOOD Moller PAIR CONDITIONS /////////////////////////////////
         ///
         ///  P_{1,2} < P_max
         ///  Psum_min < Psum < Psum_max;
         ///
         /////////////////////////////////////////////////////////////////////////////////////////////
         TVector3 p1(event1.part.px[i_electron1], event1.part.py[i_electron1], event1.part.pz[i_electron1]);
         TVector3 p2(event2.part.px[i_electron2], event2.part.py[i_electron2], event2.part.pz[i_electron2]);
         double psum = (p1 + p2).Mag();
         if (p1.Mag() < one_electron_p_max && p2.Mag() < one_electron_p_max &&
             two_electron_psum_min < psum && psum < two_electron_psum_max) {

            // Found a pair, so return.
            found1 = i_electron1;
            found2 = i_electron2;
            return (true);
         }
      }
      // Move on to the next electron(s).
      i_electron2++;
      if(i_electron2 >= event2.part.pdg.size()) {
         i_electron1++;
         i_electron2 = 0;
      }
   }

   return false;
}

bool MollerMixer::Write_Mixed_Moller_Events(MiniDst &event1, int e1_electron,
                                            MiniDst &event2, int e2_electron) {
/// Write out an event to the output TTree with two gamma's in it, one from event1 and the other from event2.
/// Also populate the run_number, event_numbers etc.

   run_number = event1.run_number;
   event_number  = event1.event_number;
   event_number2 = event2.event_number;
   trigger   = event1.trigger;
   time_stamp = event1.time_stamp;
   rf_time1 = event1.rf_time1;
   rf_time2 = event2.rf_time1;  // This is a bit arbitrary. Possibly need something sensible here later.

   // Copy the two particles into the output.
   Add_Particle(event1, e1_electron);
   Add_Particle(event2, e2_electron);

   return true;
}

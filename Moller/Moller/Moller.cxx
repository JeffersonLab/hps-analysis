//
// Created by Maurik Holtrop on 7/17/20.
//

#include "Moller.h"

Moller::Moller(std::string files) {
    if (!files.empty()) {
        ch.Add(files.c_str());
    } else {
        return;
    }
}



void Moller::Setup() {
    ///
    /// Setup a standard analysis.
    ///
    /// The analysis here is the same as the 2016 Bumphunt as documented in the analysis note.
    ///
    std::cout << "Setting up version " << __MOLLER_VERSION__ << "\n";
     // Step 1 - Select only events with 2 and only 2 electrons.
}

void Moller::Process() {
    ///
    /// Trigger Processing of all the filter steps and print a summary.
    ///
    std::cout << "Start processing.\n";
    TStopwatch timer;
    timer.Start();
    Print();
    timer.Stop();
    std::cout << "Process time was " << timer.RealTime() << " s  CPU: " << timer.CpuTime() << " s \n";
}

RNode Moller::Pass(RNode in, bool exact) {
   // Cut to select events that have two electrons. If exact = true (default) then cut on two and only
   // two electrons.
   // Parameters:  in     - Input RNode (RDateFrame)
   //              exact  - If true select exactly 2 electrons.

   return in.Filter("int n; for(int i=0;i<part_pdg.size();++i){if(part_pdg[i]==11) n++;} return n==2;");
}


RNode Moller::Cut_2_electrons(RNode in, bool exact) {
   // Cut to select events that have two electrons. If exact = true (default) then cut on two and only
   // two electrons.
   // Parameters:  in     - Input RNode (RDateFrame)
   //              exact  - If true select exactly 2 electrons.

   auto cut_two_electrons = [exact](RVec<int> part_pdg)->bool {
      int n=0;
      for(int i=0;i<part_pdg.size();++i){
         if(part_pdg[i]==11) n++;
      }
      return (exact && n==2) || (!exact && n>=2);
   };

   return in.Filter(cut_two_electrons, {"part_pdg"} );
}

RNode Moller::Select_El_Pairs(RNode in, double time_cut, std::string out_name){
   // Select all pairs of electrons from the event that have track_time within time_cut = 2ns.
   // Return the RDataFrame with an "electron_pairs" variable that are pair<int,int> type which are indexes to the part banks.

   auto select_el_pairs = [time_cut](RVec<int> part_pdg, RVec<int> part_track, RVec<double> track_time)-> std::vector< std::pair<int,int>> {
      std::vector< std::pair<int,int> > out;
      for(int i=0; i<part_pdg.size(); ++i)
         for(int j=i+1; j< part_pdg.size(); ++j){
            if(part_pdg[i] == 11 && part_pdg[j] == 11){
               int ii = part_track[i];
               int jj = part_track[j];
               if( abs(track_time[ii] - track_time[jj]) <  time_cut){
                  out.emplace_back( std::make_pair(i,j));
               }
            }
         }
      return out;
   };

   return in.Define(out_name, select_el_pairs, {"part_pdg", "part_track", "track_time"});
}

RNode Moller::Anti_Select_El_Pairs(RNode in, double time_cut, std::string out_name){
   // Select all pairs of electrons from the event that have track_time within time_cut = 2ns.
   // Return the RDataFrame with an "electron_pairs" variable that are pair<int,int> type which are indexes to the part banks.

   auto select_el_pairs = [time_cut](RVec<int> part_pdg, RVec<int> part_track, RVec<double> track_time)-> std::vector< std::pair<int,int>> {
      std::vector< std::pair<int,int> > out;
      for(int i=0; i<part_pdg.size(); ++i)
         for(int j=i+1; j< part_pdg.size(); ++j){
            if(part_pdg[i] == 11 && part_pdg[j] == 11){
               int ii = part_track[i];
               int jj = part_track[j];
               if( abs(track_time[ii] - track_time[jj]) > time_cut){
                  out.emplace_back( std::make_pair(i,j));
               }
            }
         }
      return out;
   };

   return in.Define(out_name, select_el_pairs, {"part_pdg", "part_track", "track_time"});
}



RNode Moller::Add_Momentum_Sum(RNode in, std::string in_name, std::string out_name){
   // Add the momentum sum (p1 + p2).Mag() for each pair in in_name as out_name

   auto add_psum = [](std::vector< std::pair<int, int> > pairs, RVec<double> part_px, RVec<double> part_py, RVec<double> part_pz ){
      std::vector<double> out;
      for(auto p: pairs){
         TVector3 p1(part_px[p.first], part_py[p.first], part_pz[p.first]);
         TVector3 p2(part_px[p.second], part_py[p.second], part_pz[p.second]);
         double psum = (p1+p2).Mag();
         out.push_back(psum);
      }
      return out;
   };

   return in.Define(out_name, add_psum, {in_name, "part_px", "part_py", "part_pz"});
}

RNode Moller::Refine_El_Pairs_1(RNode in, double pmin, double pmax, std::string in_name, std::string out_name){
   // Momentum sum cut for the electron pairs found in "in_name" collection of pairs: pmin < p1+p2 < pmax
   auto mom_sum_cut = [pmin, pmax](std::vector< std::pair<int, int> > pairs, RVec<double> part_px, RVec<double> part_py, RVec<double> part_pz ){
      std::vector< std::pair<int,int> > out;
      for(auto p: pairs){
         TVector3 p1(part_px[p.first], part_py[p.first], part_pz[p.first]);
         TVector3 p2(part_px[p.second], part_py[p.second], part_pz[p.second]);
         double psum = (p1+p2).Mag();
         if( pmin < psum && psum < pmax){
            out.push_back(p);
         }
      }
      return out;
   };

   return in.Define(out_name, mom_sum_cut, {in_name, "part_px", "part_py", "part_pz"});
};

RNode Moller::Refine_El_Pairs_2(RNode in, std::string in_name, std::string out_name){
   // Cut on the electron pairs, requiring one in the top and one in the bottom volume of HPS:  tan_lambda_1 * tan_lambra_2 < 0 => part_py_1 * part_py_2 < 0

   auto track_tan_lambda_cut = [](std::vector< std::pair<int, int> > pairs, RVec<double> part_py){
      std::vector< std::pair<int,int> > out;
      for(auto p: pairs){
         if( part_py[p.first]*part_py[p.second] < 0 ){
            out.push_back(p);
         }
      }
      return out;
   };

   return in.Define(out_name, track_tan_lambda_cut, {in_name, "part_py"});
};

RNode Moller::Refine_El_Pairs_3(RNode in, double min_theta, double max_theta,std::string in_name, std::string out_name){
   // Cut on the sum of the theta angle for the tracks between min_theta and max_theta, see Bradley's thesis.

   auto track_theta_cut = [min_theta, max_theta](std::vector< std::pair<int,int> > pairs, RVec<int> part_track, RVec<double> track_theta){
      std::vector< std::pair<int,int> > out;
      for(auto p: pairs){
         double angle_sum = track_theta[part_track[p.first]] + track_theta[part_track[p.second]];
         if( min_theta < angle_sum && angle_sum < max_theta){
            out.push_back(p);
         }
      }
      return out;
   };

   return in.Define(out_name, track_theta_cut, {in_name, "part_track", "track_theta"});
}


RNode Moller::Add_Moller_Inv_Mass(RNode in, std::string pair_name, std::string out_name){
   // Add the invariant mass for the two electrons in the event.
   // If there are more than two pairs, all combinatorics are used to compute the invariant mass.

   auto compute_invariant_mass = [pair_name](std::vector<std::pair<int,int> > &ele_pairs, RVec<double> &px, RVec<double> &py, RVec<double> &pz){
      std::vector<double> mass_out;
      for(const std::pair<int, int> &p: ele_pairs){
               TLorentzVector p1, p2;
               p1.SetXYZM(px[p.first],py[p.first],pz[p.first],0.0005109989499961642);
               p2.SetXYZM(px[p.second],py[p.second],pz[p.second],0.0005109989499961642);
               mass_out.push_back( (p1+p2).M());
      }
      return mass_out;
   };

   return in.Define(out_name,compute_invariant_mass,{pair_name, "part_px", "part_py", "part_pz"});
}




void Moller::Print(Option_t *opt) {

}
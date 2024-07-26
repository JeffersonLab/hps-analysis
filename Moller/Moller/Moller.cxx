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

RNode Moller::Add_Four_Vectors(RNode in, double y_rotation, std::string pair_name, std::string part_name,
                               std::string out_name_prefix) {
   // Add two four-vectors and related variables to the dataframe.

   auto add_4vector_first = [y_rotation](RVec<std::pair<int,int> > pairs, RVec<double> px, RVec<double> py, RVec<double> pz){
      RVec<TLorentzVector> out;
      for(int i=0; i<pairs.size(); ++i) {
         TLorentzVector v4;
         int item = pairs[i].first;
         v4.SetXYZM(px[item], py[item], pz[item], 0.000510998949996);
         v4.RotateY(y_rotation);
         out.emplace_back(v4);
      }
      return out;
   };

   auto add_4vector_second = [y_rotation](RVec<std::pair<int,int> > pairs, RVec<double> px, RVec<double> py, RVec<double> pz){
      RVec<TLorentzVector> out;
      for(int i=0; i<pairs.size(); ++i) {
         TLorentzVector v4;
         int item = pairs[i].second;
         v4.SetXYZM(px[item], py[item], pz[item], 0.000510998949996);
         v4.RotateY(y_rotation);
         out.emplace_back(v4);
      }
      return out;
   };


   auto add_minv = [](RVec<TLorentzVector> v1, RVec<TLorentzVector> v2){
      assert(v1.size() == v2.size());
      RVec<double> out;
      for(int i=0; i< v1.size(); ++i) {
         out.push_back((v1[i] + v2[i]).M());
      }
      return out;
   };

   auto add_theta = [](RVec<TLorentzVector> v1){
      RVec<double> out;
      for(int i=0; i< v1.size(); ++i) {
         out.push_back(v1[i].Theta());
      }
      return out;
   };

   auto add_energy = [](RVec<TLorentzVector> v1){
      RVec<double> out;
      for(int i=0; i< v1.size(); ++i) {
         out.push_back(v1[i].E());
      }
      return out;
   };

   auto add_beam_theta_x = [](RVec<TLorentzVector> v1, RVec<TLorentzVector> v2){
      assert(v1.size() == v2.size());
      RVec<double> out;
      for(int i=0; i< v1.size(); ++i) {
         out.push_back( atan2((v1[i]+v2[i]).X(), (v1[i]+v2[i]).Z()));
      }
      return out;
   };

   in = in.Define(out_name_prefix+"v1", add_4vector_first, {pair_name, part_name + "px", part_name + "py", part_name + "pz"});
   in = in.Define(out_name_prefix+"v2", add_4vector_second, {pair_name, part_name + "px", part_name + "py", part_name + "pz"});
   in = in.Define(out_name_prefix+"minv", add_minv, {out_name_prefix+"v1", out_name_prefix+"v2"});
   in = in.Define(out_name_prefix+"tht1", add_theta, {out_name_prefix+"v1"});
   in = in.Define(out_name_prefix+"tht2", add_theta, {out_name_prefix+"v2"});
   in = in.Define(out_name_prefix+"E1", add_energy, {out_name_prefix+"v1"});
   in = in.Define(out_name_prefix+"E2", add_energy, {out_name_prefix+"v2"});

   in = in.Define(out_name_prefix+"beam_theta_x", add_beam_theta_x, {out_name_prefix+"v1", out_name_prefix+"v2"});

   return in;
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

RNode Moller::Select_El_Pairs(RNode in, double time_cut_min, double time_cut_max, std::string out_name){
   // Select all pairs of electrons from the event that have track_time within time_cut_min (0ns) <= time_diff <= time_cut_max (2ns).
   // To anti-select electron pairs, set the time_cut_min > 2.ns and time_cut_max large.
   // Return the RDataFrame with an "electron_pairs" variable that are pair<int,int> type which are indexes to the part banks.

   auto select_el_pairs = [time_cut_min, time_cut_max](RVec<int> part_pdg, RVec<int> part_track, RVec<double> track_time)-> RVec< std::pair<int,int>> {
      RVec< std::pair<int,int> > out;
      for(int i=0; i<part_pdg.size(); ++i)
         for(int j=i+1; j< part_pdg.size(); ++j){
            if(part_pdg[i] == 11 && part_pdg[j] == 11){
               int ii = part_track[i];
               int jj = part_track[j];
               if( abs(track_time[ii] - track_time[jj]) >= time_cut_min && abs(track_time[ii] - track_time[jj]) <= time_cut_max){
                  out.emplace_back( std::make_pair(i,j));
               }
            }
         }
      return out;
   };

   return in.Define(out_name, select_el_pairs, {"part_pdg", "part_track", "track_time"});
}


RNode Moller::Select_El_Pairs_MC(RNode in, double time_cut_min, double time_cut_max, double z_cut, std::string out_name){
   // Select all pairs of electrons from the event that have track_time within time_cut = 2ns.
   // Return the RDataFrame with an "electron_pairs" variable that are pair<int,int> type which are indexes to the part banks.

   auto select_el_pairs_mc = [time_cut_min, time_cut_max, z_cut](RVec<int> pdg, RVec<double> time, RVec<double> z)-> RVec< std::pair<int,int>> {
      RVec< std::pair<int,int> > out;
      for(int i=0; i<pdg.size(); ++i)
         for(int j=i+1; j< pdg.size(); ++j){
            if(pdg[i] == 11 && pdg[j] == 11 && z[i] < z_cut && z[j] < z_cut){
               if( abs(time[i] - time[j]) >=  time_cut_min && abs(time[i] - time[j]) <=  time_cut_max){
                  out.emplace_back( std::make_pair(i,j));
               }
            }
         }
      return out;
   };

   return in.Define(out_name, select_el_pairs_mc, {"mc_part_pdg", "mc_part_time", "mc_part_z"});
}


RNode Moller::Add_Momentum_Sum(RNode in, std::string in_name, std::string out_name){
   // Add the momentum sum (p1 + p2).Mag() for each pair in in_name as out_name

   auto add_psum = [](RVec< std::pair<int, int> > pairs, RVec<double> part_px, RVec<double> part_py, RVec<double> part_pz ){
      RVec<double> out;
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
   // Momentum ( p_{1,2} > pmin) and momentum sum cut for the electron pairs found in "in_name" collection of pairs: pmin < p1+p2 < pmax
   auto mom_sum_cut = [pmin, pmax](RVec< std::pair<int, int> > pairs, RVec<double> part_px, RVec<double> part_py, RVec<double> part_pz ){
      RVec< std::pair<int,int> > out;
      for(auto p: pairs){
         TVector3 p1(part_px[p.first], part_py[p.first], part_pz[p.first]);
         TVector3 p2(part_px[p.second], part_py[p.second], part_pz[p.second]);
         double psum = (p1+p2).Mag();
         if( p1.Mag() < pmin && p2.Mag() < pmin && pmin < psum && psum < pmax){
            out.push_back(p);
         }
      }
      return out;
   };

   return in.Define(out_name, mom_sum_cut, {in_name, "part_px", "part_py", "part_pz"});
};

RNode Moller::Refine_El_Pairs_2(RNode in, std::string in_name, std::string out_name){
   // Cut on the electron pairs, requiring one in the top and one in the bottom volume of HPS:  tan_lambda_1 * tan_lambra_2 < 0 => part_py_1 * part_py_2 < 0

   auto track_tan_lambda_cut = [](RVec< std::pair<int, int> > pairs, RVec<double> part_py){
      RVec< std::pair<int,int> > out;
      for(auto p: pairs){
         if( part_py[p.first]*part_py[p.second] < 0 ){
            out.push_back(p);
         }
      }
      return out;
   };

   return in.Define(out_name, track_tan_lambda_cut, {in_name, "part_py"});
};

RNode Moller::Refine_El_Pairs_2(RNode in, double phi_cut, std::string in_name, std::string p4_prefix, std::string out_name){
   // Cut on the electron pairs requiring the tracks to be co-planar to better than phi_cut:  abs( t1.Phi() - t2.Phi() ) < phi_cut.

   auto track_tan_lambda_cut = [phi_cut](RVec<std::pair<int,int> > pairs, RVec<TLorentzVector> p4v1, RVec<TLorentzVector> p4v2){
      RVec< std::pair<int,int> > out;
      for(int i=0; i<pairs.size(); ++i){
         double phi_diff = p4v1[i].Phi()-p4v2[i].Phi(); phi_diff = (phi_diff>0)?phi_diff: phi_diff+TMath::Pi();
         if( abs(phi_diff - TMath::Pi()/2) < phi_cut ){
            out.push_back(pairs[i]);
         }
      }
      return out;
   };

   return in.Define(out_name, track_tan_lambda_cut, {in_name, p4_prefix+"p4v1", p4_prefix+"p4v2"});
};

RNode Moller::Refine_El_Pairs_3(RNode in, double theta_cut, std::string in_name, std::string p4_prefix, std::string out_name) {
   // Make a cut on the calculated beam direction, beam.Theta() = (p4v1+p4v2).Theta(). Note: this is different from theta_x.

   auto beam_theta_cut = [theta_cut](RVec<std::pair<int,int> > pairs,  RVec<TLorentzVector> p4v1, RVec<TLorentzVector> p4v2){
      RVec< std::pair<int,int> > out;
      for(int i=0; i<pairs.size(); ++i){
         if( (p4v1[i]+p4v2[i]).Theta() < theta_cut ){
            out.push_back(pairs[i]);
         }
      }
      return out;
   };

   return in.Define(out_name, beam_theta_cut, {in_name, p4_prefix+"p4v1", p4_prefix+"p4v2" });
}

RNode Moller::Refine_El_Pairs_X(RNode in, double min_theta, double max_theta, std::string angle1, std::string angle2, std::string in_name, std::string out_name){
   // Cut on the sum of the theta angle for the tracks between min_theta and max_theta, see Bradley's thesis.
   // This is probably a cut WE SHOULD NOT DO, since it will also select on invariant mass indirectly.

   auto track_theta_cut = [min_theta, max_theta](RVec< std::pair<int, int> > pairs, RVec<double> p4tht1, RVec<double> p4tht2){
      RVec< std::pair<int,int> > out;
      for(int i=0; i<pairs.size(); ++i){
         double angle_sum = p4tht1[i] + p4tht2[i];
         if( min_theta < angle_sum && angle_sum < max_theta){
            out.push_back(pairs[i]);
         }
      }
      return out;
   };

   return in.Define(out_name, track_theta_cut, {in_name, angle1, angle2});
}


RNode Moller::Add_Moller_Inv_Mass(RNode in, std::string pair_name, std::string out_name){
   // Add the invariant mass for the two electrons in the event.
   // If there are more than two pairs, all combinatorics are used to compute the invariant mass.

   auto compute_invariant_mass = [pair_name](RVec<std::pair<int,int> > &ele_pairs, RVec<double> &px, RVec<double> &py, RVec<double> &pz){
      RVec<double> mass_out;
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
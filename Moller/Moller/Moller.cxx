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
    /// Fill this in if you want code to process data.
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

   auto add_psum = [](RVec<TLorentzVector> v1, RVec<TLorentzVector> v2){
      assert(v1.size() == v2.size());
      RVec<double> out;
      for(int i=0; i< v1.size(); ++i) {
         out.push_back( (v1[i].Vect() + v2[i].Vect()).Mag());
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
   in = in.Define(out_name_prefix+"psum", add_psum, {out_name_prefix+"v1", out_name_prefix+"v2"});
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

   auto cut_two_electrons = [exact](const RVec<int> &part_pdg)->bool {
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

   auto select_el_pairs = [time_cut_min, time_cut_max]
         (const RVec<int> &part_pdg, const RVec<int> &part_track, const RVec<double> &track_time)-> RVec< std::pair<int,int>> {
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

   // Note that without cuts, the index is just a dummy!
   return in.Define(out_name, select_el_pairs, {"part_pdg", "part_track", "track_time"})
      .Define(out_name+"_idx","RVec<int> out;for(int i=0; i< "+out_name+".size(); ++i) out.push_back(i);return out;");
}

RNode Moller::Select_v0(RNode in, int type, double time_cut_min, double time_cut_max, std::string out_name){
   // Select vertexes from the v0 collection of type "type" that have a track time difference time_cut_min
   //  <= time_diff <= time_cut_max.
   // Return the RDataFrame with an "out_name" variable that are pair<int,int> of the particles in the v0_em_part and
   // v0_ep_part. This is indeed a duplication, and wasteful, but makes using the v0_ collection easier and more
   // compatible with the select_ele_pairs function.
   // The "out_name_idx" variable that are the indexes of the selected v0_xxxx variables that we want.

   auto add_electron_pairs = [](RVec<int> v0_em_part, RVec<int> v0_ep_part){
      RVec< std::pair<int,int> > out;
      for(int i=0; i<v0_em_part.size(); ++i){
         out.emplace_back( std::make_pair(v0_em_part[i], v0_ep_part[i]));
      }
      return out;
   };

   auto select_v0 = [type, time_cut_min, time_cut_max]
         (RVec<int> v0_type, const RVec<int>& v0_em_track, const RVec<int>& v0_ep_track,
               const RVec<double>& track_time)-> RVec<int> {
      RVec<int> out;
      for(int i=0; i<v0_type.size(); ++i)
         if(v0_type[i] == type){
            int ii = v0_em_track[i];
            int jj = v0_ep_track[i];
            if (abs(track_time[ii] - track_time[jj]) >= time_cut_min &&
            abs(track_time[ii] - track_time[jj]) <= time_cut_max) {
               out.push_back(i);
            }
         }
      return out;
   };

   return in.Define(out_name,add_electron_pairs,{"v0_em_part","v0_ep_part"})
   .Define(out_name+"_idx", select_v0, {"v0_type", "v0_em_track", "v0_ep_track", "track_time"});
}

RNode Moller::Select_El_Pairs_MC(RNode in, double time_cut_min, double time_cut_max, double z_cut, std::string out_name){
   // Select all pairs of electrons from the MCPartcile collection that are primary electrons and
   // have mc_part_time within time_cut_min <= time_diff <= time_cut_max.

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


//RNode Moller::Add_Momentum_Sum(RNode in, std::string in_name, std::string out_name){
//   // Add the momentum sum (p1 + p2).Mag() for each pair in in_name as out_name
//
//   auto add_psum = [](RVec< std::pair<int, int> > pairs, RVec<double> part_px, RVec<double> part_py, RVec<double> part_pz ){
//      RVec<double> out;
//      for(auto p: pairs){
//         TVector3 p1(part_px[p.first], part_py[p.first], part_pz[p.first]);
//         TVector3 p2(part_px[p.second], part_py[p.second], part_pz[p.second]);
//         double psum = (p1+p2).Mag();
//         out.push_back(psum);
//      }
//      return out;
//   };
//
//   return in.Define(out_name, add_psum, {in_name, "part_px", "part_py", "part_pz"});
//}

RNode Moller::Refine_El_Pairs_1(RNode in, double pcut, double pmin, double pmax,
                                std::string pairs_name, std::string in_name_idx, std::string out_name_idx){
   // Momentum ( p_{1,2} > pcut) and momentum sum cut (pmin < p1+p2 < pmax) for the electron pairs found in "in_name" index
   // Creates a new index "out_name" with the indexes pointing to the refined selection.
   auto mom_sum_cut = [pcut, pmin, pmax](RVec< std::pair<int, int> > pairs, RVec<int> in_idx, RVec<double> part_px, RVec<double> part_py, RVec<double> part_pz ){
      RVec<int> out;
      for(int idx: in_idx){
         int i_first = pairs[idx].first;
         int i_second = pairs[idx].second;
         TVector3 p1(part_px[i_first], part_py[i_first], part_pz[i_first]);
         TVector3 p2(part_px[i_second], part_py[i_second], part_pz[i_second]);
         double psum = (p1+p2).Mag();
         if( p1.Mag() < pcut && p2.Mag() < pcut && pmin < psum && psum < pmax){
            out.push_back(idx);
         }
      }
      return out;
   };

   return in.Define(out_name_idx, mom_sum_cut, {pairs_name, in_name_idx, "part_px", "part_py", "part_pz"});
};

RNode Moller::Refine_El_Pairs_2(RNode in, std::string pairs_name, std::string in_idx_name, std::string out_idx_name){
   // Cut on the electron pairs, requiring one in the top and one in the bottom volume of HPS:
   // tan_lambda_1 * tan_lambra_2 < 0 => part_py_1 * part_py_2 < 0

   auto track_tan_lambda_cut = [](RVec< std::pair<int, int> > pairs, RVec<int> in_idx, RVec<double> part_py){
      RVec<int> out;
      for(int idx: in_idx){
         if( part_py[pairs[idx].first]*part_py[pairs[idx].second] < 0 ){
            out.push_back(idx);
         }
      }
      return out;
   };

   return in.Define(out_idx_name, track_tan_lambda_cut, {pairs_name, in_idx_name, "part_py"});
};

RNode Moller::Refine_El_Pairs_2(RNode in, double phi_cut, std::string pairs_name, std::string in_idx_name, std::string p4_prefix, std::string out_idx_name){
   // Cut on the electron pairs requiring the tracks to be co-planar to better than phi_cut:  abs( t1.Phi() - t2.Phi() ) < phi_cut.
   // The p4_prefix is the TLorentzVectors collection for *all* the electron_pairs. The in_idx_name is the input indexes to the pairs.

   auto track_tan_lambda_cut = [phi_cut](RVec<std::pair<int,int> > pairs, RVec<int> in_idx, RVec<TLorentzVector> p4v1, RVec<TLorentzVector> p4v2){
      RVec<int> out;
      for(int idx: in_idx){
         double phi_diff = p4v1[idx].Phi()-p4v2[idx].Phi();
         phi_diff = phi_diff< -TMath::Pi()?phi_diff + 2*TMath::Pi(): phi_diff;
         if(abs(phi_diff) < phi_cut) out.push_back(idx);
      }
      return out;
   };

   return in.Define(out_idx_name, track_tan_lambda_cut, {pairs_name, in_idx_name, p4_prefix+"v1", p4_prefix+"v2"});
};

RNode Moller::Refine_El_Pairs_3(RNode in, double theta_cut, std::string in_idx_name, std::string p4_prefix, std::string out_idx_name) {
   // Make a cut on the calculated beam direction, beam.Theta() = (p4v1+p4v2).Theta(). Note: this is different from theta_x.

   auto beam_theta_cut = [theta_cut](RVec<int> in_idx, RVec<TLorentzVector> p4v1, RVec<TLorentzVector> p4v2){
      RVec<int> out;
      for(int idx: in_idx){
         if( (p4v1[idx].Vect()+p4v2[idx].Vect()).Theta() < theta_cut ){
            out.push_back(idx);
         }
      }
      return out;
   };

   return in.Define(out_idx_name, beam_theta_cut, {in_idx_name, p4_prefix+"v1", p4_prefix+"v2" });
}

RNode Moller::Refine_El_Pairs_X(RNode in, double min_theta, double max_theta, std::string angle1, std::string angle2,
                                std::string in_idx_name, std::string out_idx_name){
   // Cut on the sum of the theta angle for the tracks between min_theta and max_theta, see Bradley's thesis.
   // This is probably a cut WE SHOULD NOT DO, since it will also select on invariant mass indirectly.

   auto track_theta_cut = [min_theta, max_theta](RVec<int> in_idx, RVec<double> p4tht1, RVec<double> p4tht2){
      RVec<int> out;
      for(int ind: in_idx){
         double angle_sum = p4tht1[ind] + p4tht2[ind];
         if( min_theta < angle_sum && angle_sum < max_theta){
            out.push_back(ind);
         }
      }
      return out;
   };

   return in.Define(out_idx_name, track_theta_cut, {in_idx_name, angle1, angle2});
}

RNode Moller::Refine_El_Pairs_X2(RNode in, std::string pairs_name, std::string in_idx_name, std::string out_idx_name, bool isData){
   // This is the cut that uses the fiducial cut on the tracks as defined in Omar's Analysis Note, Appendix F.
   //
   auto TrackIsFiducial = [isData](double tr_x, double tr_y) -> bool {
      bool fiducial = false;
      if (isData) {
         if (tr_y > 0) {
            if (tr_y < 42. && tr_y < 13 - 0.26 * tr_x && tr_y > 18 - 0.08 * tr_x &&
                ((tr_x > -125. && tr_x < -95.) || (tr_x > -85 && tr_x < -55.))) {
               fiducial = true;
            }
         } else if (tr_y < 0) {
            if (tr_y < -23. && tr_y < -15. + 0.08 * tr_x && tr_y > -18. + 0.22 * tr_x &&
                ((tr_x > -75. && tr_x < -45) || (tr_x > -110. && tr_x < -95.))) {
               fiducial = true;
            }
         }
      } else {
         if (tr_y > 0) {
            if (tr_y > 23. && tr_y > 15 - 0.1 * tr_x && tr_y < 12 - 0.3 * tr_x && ((tr_x > -75. && tr_x
                                                                                                   < -50.) ||
                                                                                   (tr_x > -130. && tr_x < -95))) {
               fiducial = true;
            }
         } else if (tr_y < 0) {
            if (tr_y < -22 && tr_y < -15 + 0.1 * tr_x && tr_y > -15 + 0.25 * tr_x &&
                ((tr_x > -120. && tr_x < -94) || (tr_x > -75 && tr_x < -50.))) {
               fiducial = true;
            }
         }
      }
      return fiducial;
   };

   auto TracksFiducialTest = [TrackIsFiducial]( RVec<std::pair<int,int> > pairs, RVec<int> in_idx, RVec<double> track_x_at_ecal , RVec<double> track_y_at_ecal) {
      RVec<int> out;
      for(int ind: in_idx){
         int first = pairs[ind].first;
         int second = pairs[ind].second;
         if( TrackIsFiducial(track_x_at_ecal[first], track_y_at_ecal[first]) &&
             TrackIsFiducial(track_x_at_ecal[second], track_y_at_ecal[second]) ){
            out.push_back(ind);
         }
      }
      return out;
   };

   return in.Define(out_idx_name, TracksFiducialTest, {pairs_name, in_idx_name, "track_x_at_ecal", "track_y_at_ecal"});
}

RNode Moller::Refine_El_Pairs_chi2(RNode in, double chi2_cut, std::string pairs_name, std::string in_idx_name, std::string out_idx_name){
// Cut the electron pairs on the chi2 of each track.

      auto track_chi2_cut = [chi2_cut](RVec< std::pair<int, int> > pairs, RVec<int> in_indx, RVec<double> part_track_chi2){
         RVec<int> out;
         for(int ind: in_indx){
            if( part_track_chi2[pairs[ind].first] < chi2_cut && part_track_chi2[pairs[ind].second] < chi2_cut){
               out.push_back(ind);
            }
         }
         return out;
      };
      return in.Define(out_idx_name, track_chi2_cut, {pairs_name, in_idx_name, "part_track_chi2"});
}

RNode Moller::Refine_El_Pairs_nhit(RNode in, int nhit_cut, std::string in_idx_name, std::string out_idx_name){
// Cut the electron pairs on the number of hits in each track.
   auto track_nhit_cut = [nhit_cut](RVec< std::pair<int, int> > pairs, RVec<int> in_indx, RVec<int> track_nhits, RVec<int> part_track){
      RVec<int> out;
      for(int ind: in_indx){
         if( track_nhits[part_track[pairs[ind].first]] > nhit_cut && track_nhits[part_track[pairs[ind].second]] > nhit_cut){
            out.push_back(ind);
         }
      }
      return out;
   };
   return in.Define(out_idx_name, track_nhit_cut, {in_idx_name, "track_nhits", "part_track"});
}



//RNode Moller::Add_Moller_Inv_Mass(RNode in, std::string pair_name, std::string out_name){
//   // Add the invariant mass for the two electrons in the event.
//   // If there are more than two pairs, all combinatorics are used to compute the invariant mass.
//
//   auto compute_invariant_mass = [pair_name](RVec<std::pair<int,int> > &ele_pairs, RVec<double> &px, RVec<double> &py, RVec<double> &pz){
//      RVec<double> mass_out;
//      for(const std::pair<int, int> &p: ele_pairs){
//               TLorentzVector p1, p2;
//               p1.SetXYZM(px[p.first],py[p.first],pz[p.first],0.0005109989499961642);
//               p2.SetXYZM(px[p.second],py[p.second],pz[p.second],0.0005109989499961642);
//               mass_out.push_back( (p1+p2).M());
//      }
//      return mass_out;
//   };
//
//   return in.Define(out_name,compute_invariant_mass,{pair_name, "part_px", "part_py", "part_pz"});
//}

double Moller::mycrystalball(double x, double N, double alpha, double n, double sigma, double mean) {
   double z = (x - mean)/sigma;
   if (alpha < 0) z = -z;
   double abs_alpha = fabs(alpha);
   double C = n/fabs(alpha)*(1/(n-1)*exp(-0.5*abs_alpha*abs_alpha));
   double D = sqrt(M_PI/2)*(1 + erf(abs_alpha/sqrt(2)));
   double Nprime = N/(sigma*(C+D));
   if (z > -abs_alpha) {
      return Nprime*exp(-0.5*z*z);
   } else {
      double A = pow(n/abs_alpha,n) * exp(-0.5*abs_alpha*abs_alpha);
      double B = n/abs_alpha - abs_alpha;
      return Nprime*A/pow(B - z, n);
   }
}

void Moller::Print(Option_t *opt) {

}
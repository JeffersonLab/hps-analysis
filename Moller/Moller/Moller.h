//
// Created by Maurik Holtrop on 7/17/20.
//

// References:
// "Resonance search analysis of 2016 HPS spring run data." Bump hunt folks, https://confluence.slac.stanford.edu/download/attachments/146715820/HPS_2016_Bump_Hunt_Internal_Note.pdf
// "Search for a Heavy Photon in Electro-Produced e+e− Pairs with the HPS Experiment at JLab", Omar Moreno, Nathan Baltzell, Mathew Graham, John Jaros, https://confluence.slac.stanford.edu/download/attachments/146715820/engrun2015_resonance_search.pdf


#ifndef MOLLER_MOLLER_H
#define MOLLER_MOLLER_H

#define __MOLLER_VERSION__ "0.3.2"

#include "TObject.h"
#include "TChain.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDF/HistoModels.hxx"

using namespace ROOT;
using RNode = ROOT::RDF::RNode;

class Moller {

public:

   double e_time_cut_min{0.};
   double e_time_cut_max{2.};
   double E_beam{2.3};
   double e_maximum_momentum{0.75*E_beam};  //
   double ee_mom_sum_min{2.1};      // From 2016 note, table 8. Alternate is 0.8*E_beam
   double ee_mom_sum_max{2.45};       // ... Alternate is 1.18*E_beam
   double four_vector_rotation{-0.0302};


   Moller(std::string files = "");

   void Setup();
   void Process();
   void Print(Option_t *opt = "");

   RNode Select_El_Pairs(RNode in, double time_cut=4., double time_cut_max=20., std::string out_name="electron_pairs");
//   RNode Select_El_Pairs(RNode in, double time_cut=2., std::string out_name="electron_pairs");
   RNode Select_El_Pairs_MC(RNode in, double time_cut_min = 0., double time_cut_max=2., double z_cut = 0.001, std::string out_name="mc_electron_pairs");

   RNode Add_Four_Vectors(RNode in, double y_rotation=-0.0302, std::string pair_name="electron_pairs", std::string part_name="part_", std::string out_name_prefix="p4");
   RNode Add_Four_Vectors_V0(RNode in, double y_rotation=-0.0302, int v0_type=4, std::string out_name_prefix="v_uc");
   RNode Add_Moller_Inv_Mass(RNode in, std::string pair_name="electron_pairs", std::string out_name="moller_inv_mass");
   RNode Add_Momentum_Sum(RNode in, std::string pair_name="electron_pair", std::string out_name="el_pair_psum");
   RNode Refine_El_Pairs_1(RNode in, double pcut=2., double pmin=2.306*0.75, double pmax=2.306*1.15, std::string in_name="electron_pairs", std::string out_name="el_pairs_r1");
   RNode Refine_El_Pairs_2(RNode in, std::string in_name="el_pairs_r1", std::string out_name="el_pairs_r2");
   RNode Refine_El_Pairs_2(RNode in, double phi_cut, std::string in_name="el_pairs_r1", std::string p4_prefix="r1_", std::string out_name="el_pairs_r2");
   RNode Refine_El_Pairs_3(RNode in, double theta_cut, std::string in_name="el_pairs_r2", std::string p4_prefix="r2_", std::string out_name="el_pairs_r3");

   RNode Refine_El_Pairs_X(RNode in, double mint_heta=0.04, double max_theta=0.048, std::string angle1="r2_p4tht1", std::string angle2="r2_p4tht2", std::string in_name="el_pairs_r2", std::string out_name="el_pairs_r3");
   RNode Refine_El_Pairs_X2(RNode in,std::string in_name, std::string out_name, bool isData=true);
   RNode Cut_2_electrons(RNode in, bool exact = true);

   // Variadic function permits any number of arguments. Unfortunately, the RDataFrame.Define or .Filter do not work with this.
   template<typename ...Vecs>
   static ROOT::VecOps::RVec<bool>
   Combine_cuts_helper(ROOT::VecOps::RVec<bool> vec1, ROOT::VecOps::RVec<bool> vec2, Vecs ...vecs) {
      /// Make a logic and combination of the input vectors and return result.
      // The function makes use of C++17 variadic function templates and parameter packs, so that
      // any number of input variables can be given. It works, but bot with RDataFrame.Define() :-(

      RVec<bool> out;
      std::vector<RVec<bool>> inputs;
      (inputs.push_back(vecs), ...);  // C++17 parameter pack: push all argument vectors into a vector.

      for (int i = 0; i < inputs[0].size(); ++i) { // Combine
         // the inputs into the output.
         bool is_ok = vec1[i] && vec2[i] && (vecs[i] && ...); // parameter pack expansion combines all ...vecs
         out.push_back(is_ok);
      }
      return out;
   }

   static RVec<bool> Combine_cuts_helper2(RVec<bool> vec1, RVec<bool> vec2) { // Permit any number of arguments.
      RVec<bool> out;
      for (int i = 0; i < vec1.size(); ++i) { // Combe the inputs into the output.
         out.push_back(vec1[i] && vec2[i]);
      }
      return out;
   }

   static RVec<bool>
   Combine_cuts_helper3(RVec<bool> vec1, RVec<bool> vec2, RVec<bool> vec3) { // Permit any number of arguments.
      RVec<bool> out;
      for (int i = 0; i < vec1.size(); ++i) { // Combe the inputs into the output.
         out.push_back(vec1[i] && vec2[i] && vec3[i]);
      }
      return out;
   }

   static RVec<bool> Combine_cuts_helper4(RVec<bool> vec1, RVec<bool> vec2, RVec<bool> vec3,
                                          RVec<bool> vec4) { // Permit any number of arguments.
      RVec<bool> out;
      for (int i = 0; i < vec1.size(); ++i) { // Combe the inputs into the output.
         out.push_back(vec1[i] && vec2[i] && vec3[i] && vec4[i]);
      }
      return out;
   }

   static RVec<double> Vector_magnitude(RVec<double> x, RVec<double> y, RVec<double> z) {
      /// Return the magnitude of the x,y,z vector for each entry in x,y,z.
      RVec<double> out;
      for (int i = 0; i < x.size(); ++i) {
         out.push_back(sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]));
      }
      return out;
   }

   static RNode Combine_cuts(RNode in, std::string_view out_name, std::vector<std::string> cuts) {
      ///
      /// For the in RNode, Combine the named bool columns cuts into a single new named column 'out_name' using logical and (&&).
      /// Example:
      /// auto combo = Combine_cuts(dataframe, "v0_combo_cut", {"v0_clus_time_cut", "v0_d0_cut", "v0_track_clus_time_cut"} )
      /// this will combine "v0_clus_time_cut", "v0_d0_cut", "v0_track_clus_time_cut" item for item into "v0_combo_cut".
      ///
      // Note: this can also be implemented as a template<typename ...Vecs> Combine_cuts(RNode in, std::string_view out_name, ...vecs)
      // Unfortunately the RDataFrame.Define() seems to not like being passed a variadic function, so we need the switch either way.
      //
      switch (cuts.size()) {
         case 0:
         case 1:
            std::cout << "ERROR - You need at least two columns to combine.";
            break;
         case 2:
            return in.Define(out_name, Combine_cuts_helper2, cuts);
         case 3:
            return in.Define(out_name, Combine_cuts_helper3, cuts);
         case 4:
            return in.Define(out_name, Combine_cuts_helper4, cuts);
            break;
         default:
            std::cout << "I did not expect anyone to ever need to cut " << cuts.size() << " all at once!\n";
      }
      return in;
   }

   static std::string GetVersion(){
      return( std::string(__MOLLER_VERSION__));
   };

public:
   TChain ch{TChain("MiniDST")};
   RDataFrame dataframe{RDataFrame(ch)};

};

#endif //MOLLER_MOLLER_H

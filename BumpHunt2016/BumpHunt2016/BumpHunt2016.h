//
// Created by Maurik Holtrop on 7/17/20.
//

#ifndef BUMPHUNT2016_BUMPHUNT2016_H
#define BUMPHUNT2016_BUMPHUNT2016_H

#include "TObject.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDF/HistoModels.hxx"

#include "FilterStore.h"

using namespace ROOT;
using RNode = ROOT::RDF::RNode;


class BumpHunt2016 {
public:
    enum ParticleType {  // Copy of Vertex and Particle types from hps-dst.
        FINAL_STATE_PARTICLE = 0,
        UC_V0_CANDIDATE      = 1,
        BSC_V0_CANDIDATE     = 2,
        TC_V0_CANDIDATE      = 3,
        UC_MOLLER_CANDIDATE  = 4,
        BSC_MOLLER_CANDIDATE = 5,
        TC_MOLLER_CANDIDATE  = 6,
        OTHER_ELECTRONS      = 7,
        UC_VC_CANDIDATE      = 8
    };

public:
    BumpHunt2016(std::string files="");
    void Setup();
    void Process();
    void Print(Option_t *opt="");
    RNode GetLastFilterNode(){ return filters.back().node;};
    FilterStore &GetLastFilter(){ return filters.back();};
    FilterStore &GetFilter(unsigned int i){ return filters[i];};
    FilterStore &GetFilter(std::string_view name){ return filters[filtermap[name]];};
    RNode GetFilterNode(unsigned int i){ return filters[i].node;};
    RNode GetFilterNode(std::string_view name){ return filters[filtermap[name]].node;};

    // Helper Function
    RNode FilterStep(RNode df, const std::string_view name="",const std::string_view v0_status_name=""){
        std::string step_name, status_name;
        if (name.empty()) {
            step_name = df.GetFilterNames().back();
        }else{
            step_name = name;
        }
        filters.emplace_back(df, step_name, v0_status_name);
        filtermap[name]=filters.size()-1;
        return df;
    }
    // Define specific filter steps.
    RNode Select_Pair1(RNode in);
 //   RNode Correct_all_ecal_times(RNode in);
    RNode Correct_v0_clus_times(RNode in);
    RNode Cut_on_reduced_chi1(RNode in, double em_chi2_max=12., double ep_chi2_max = 12.);
    RNode Cut_on_good_pid(RNode in, double good_pid_max=10.);
    RNode Simple_track_cluster_time_cut(RNode in, double em_dt_max= 6., double ep_dt_max= 6.);
    RNode Momentum_dependent_track_cluster_time_cut(RNode in);
    RNode Cut_on_min_max_vertex_momentum(RNode in, double psum_min=1.9, double psum_max=2.4);
    RNode Cut_on_cluster_in_time(RNode in);
    RNode Cut_on_cluster_time_diff(RNode in, double dt_max=1.43);
    RNode Combine_valid_v0_and_reduce(RNode in, std::vector<std::string> args);

    // Variadic function permits any number of arguments. Unfortunately, the RDataFrame.Define or .Filter do not work with this.
    template<typename ...Vecs>
    static ROOT::VecOps::RVec<bool> Combine_cuts_helper(ROOT::VecOps::RVec<bool> vec1, ROOT::VecOps::RVec<bool> vec2, Vecs ...vecs){
        /// Make a logic and combination of the input vectors and return result.
        // The function makes use of C++17 variadic function templates and parameter packs, so that
        // any number of input variables can be given. It works, but bot with RDataFrame.Define() :-(

        RVec<bool> out;
        std::vector<RVec<bool>> inputs;
        (inputs.push_back(vecs), ...);  // C++17 parameter pack: push all argument vectors into a vector.

        for( int i=0; i< inputs[0].size(); ++i){ // Combe the inputs into the output.
            bool is_ok = vec1[i] && vec2[i] && (vecs[i] && ...); // parameter pack expansion combines all ...vecs
            out.push_back(is_ok);
        }
        return out;
    }

    static RVec<bool> Combine_cuts_helper2(RVec<bool> vec1, RVec<bool> vec2){ // Permit any number of arguments.
        RVec<bool> out;
        for( int i=0; i< vec1.size(); ++i){ // Combe the inputs into the output.
            out.push_back(vec1[i] && vec2[i]);
        }
        return out;
    }

    static RVec<bool> Combine_cuts_helper3(RVec<bool> vec1, RVec<bool> vec2, RVec<bool> vec3){ // Permit any number of arguments.
        RVec<bool> out;
        for( int i=0; i< vec1.size(); ++i){ // Combe the inputs into the output.
            out.push_back(vec1[i] && vec2[i] && vec3[i]);
        }
        return out;
    }

    static RVec<bool> Combine_cuts_helper4(RVec<bool> vec1, RVec<bool> vec2, RVec<bool> vec3, RVec<bool> vec4){ // Permit any number of arguments.
        RVec<bool> out;
        for( int i=0; i< vec1.size(); ++i){ // Combe the inputs into the output.
            out.push_back(vec1[i] && vec2[i] && vec3[i] && vec4[i]);
        }
        return out;
    }

    static RVec<double> Vector_magnitude(RVec<double> x, RVec<double> y, RVec<double> z){
        /// Return the magnitude of the x,y,z vector for each entry in x,y,z.
        RVec<double> out;
        for(int i=0; i< x.size(); ++i){
            out.push_back( sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]));
        }
        return out;
    }

    static RNode Combine_cuts(RNode in, std::string_view out_name, std::vector<std::string> cuts){
        ///
        /// For the in RNode, Combine the named bool columns cuts into a single new named column 'out_name' using logical and (&&).
        /// Example:
        /// auto combo = Combine_cuts(dataframe, "v0_combo_cut", {"v0_clus_time_cut", "v0_d0_cut", "v0_track_clus_time_cut"} )
        /// this will combine "v0_clus_time_cut", "v0_d0_cut", "v0_track_clus_time_cut" item for item into "v0_combo_cut".
        ///
        // Note: this can also be implemented as a template<typename ...Vecs> Combine_cuts(RNode in, std::string_view out_name, ...vecs)
        // Unfortunately the RDataFrame.Define() seems to not like being passed a variadic function, so we need the switch either way.
        //
        switch(cuts.size()){
            case 0:
            case 1:
                std::cout << "ERROR - You need at least two columns to combine.";
                break;
            case 2:
                return in.Define(out_name,Combine_cuts_helper2, cuts);
            case 3:
                return in.Define(out_name,Combine_cuts_helper3, cuts);
            case 4:
                return in.Define(out_name,Combine_cuts_helper4, cuts);
                break;
            default:
                std::cout << "I did not expect anyone to ever need to cut " << cuts.size() << " all at once!\n";
        }
        return in;
    }

public:
    TChain ch{TChain("MiniDST")};
    RDataFrame dataframe{RDataFrame(ch)};
    std::map<std::string_view,unsigned int> filtermap;
    std::vector<FilterStore> filters;

    double clust_track_time_offset{55.};   // Value used by Rafo in setting_2016_pass1.h:21 CL_trk_time_Offset_Data = 55.;
    std::vector<double> Tr_Cl_Time_Cut_Top_Upper_Lim{-6.56042, 27.9583, -25.605, 6.9747};     // Coefficients for 3 order poly for time cuts.
    std::vector<double> Tr_Cl_Time_Cut_Top_Lower_Lim{-1.55686, -3.51391, 4.19438, -1.08659};
    std::vector<double> Tr_Cl_Time_Cut_Bot_Upper_Lim{-6.15725, 25.4919, -22.3805, 5.73573};
    std::vector<double> Tr_Cl_Time_Cut_Bot_Lower_Lim{-3.33916, -0.537733, 2.06764, -0.629115};
    std::vector<double> Cluster_Time_Cut_Bot_Upper_Lim{58.5, 3.40282, -1.00306};
    std::vector<double> Cluster_Time_Cut_Bot_Lower_Lim{45.51, 7.55268, -1.89745};

     inline const RVec<double> xxx_ecal_time_correction(RVec<double> ecal_time, RVec<int> ecal_ix, RVec<int> ecal_iy){
        RVec<double> out;
        for(int i = 0; i< ecal_time.size(); ++i){
            out.push_back( ecal_time[i] - Ecal_Time_Correction_values[Ecal_Time_Corr_bin(ecal_ix[i],ecal_iy[i])]);
        }
        return out;
    };

    inline static unsigned int Ecal_Time_Corr_bin(int ix, int iy) { return (ix+23)*11+(iy+5); };
    static std::vector<double> Ecal_Time_Correction_values;
};

#endif //BUMPHUNT2016_BUMPHUNT2016_H

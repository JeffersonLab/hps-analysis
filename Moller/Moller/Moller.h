//
// Created by Maurik Holtrop on 7/17/20.
//

#ifndef MOLLER_MOLLER_H
#define MOLLER_MOLLER_H

#define __MOLLER_VERSION__ "0.1.1"

#include "TObject.h"
#include "TChain.h"
#include "TROOT.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDF/HistoModels.hxx"

#include "FilterStore.h"

using namespace ROOT;
using RNode = ROOT::RDF::RNode;


class Moller {

public:
    Moller(std::string files = "");

    void Setup();
    void Process();
    void Print(Option_t *opt = "");

   RNode Pass(RNode in, bool exact = true);
   RNode Cut_2_electrons(RNode in, bool exact = true);



    FilterStore &GetLastFilter() { return filters.back(); };
    FilterStore &GetFilter(unsigned int i) { return filters[i]; };
    FilterStore &GetFilter(std::string_view name) { return filters[filtermap[name]]; };
    RNode GetLastFilterNode() { return filters.back().node; };
    RNode GetFilterNode(unsigned int i) { return filters[i].node; };
    RNode GetFilterNode(std::string_view name) { return filters[filtermap[name]].node; };

    // Helper Function
    RNode FilterStep(RNode df, const std::string_view name = "", const std::string_view v0_status_name = "") {
        std::string step_name, status_name;
        if (name.empty()) {
            step_name = df.GetFilterNames().back();
        } else {
            step_name = name;
        }
        filters.emplace_back(df, step_name, v0_status_name);
        filtermap[name] = filters.size() - 1;
        return df;
    }
    // Define specific filter steps.

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
    std::map<std::string_view, unsigned int> filtermap;
    std::vector<FilterStore> filters;

};

#endif //MOLLER_MOLLER_H

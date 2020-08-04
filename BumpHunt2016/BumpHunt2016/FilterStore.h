//
// Created by Maurik Holtrop on 7/22/20.
//

#ifndef BUMPHUNT2016_FILTERSTORE_H
#define BUMPHUNT2016_FILTERSTORE_H

#include "TObject.h"
#include "TChain.h"
#include "TROOT.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDF/HistoModels.hxx"

using namespace ROOT;
using RNode = ROOT::RDF::RNode;


class FilterStore{

public:
    std::string name;
    std::string v0_status_column; // This is a hack because we cannot, yet, overwrite a column. So cuts on the v0 array
    // need a new colunm for each cut. Not a problem, but we need to keep track.
    ROOT::RDF::RNode node;

    // For convenience, we count the number of events and the number of v0 that survived after this step
    ROOT::RDF::RResultPtr<ULong64_t> Count;
    ROOT::RDF::RResultPtr<ULong64_t> V0_Count;

    // Some convenient histograms that we want to fill for each of the steps.
    const ROOT::RDF::TH1DModel M_inv_def{"Minv","Invariant Mass;M_{inv} [GeV];counts",6000,0.,0.3};
    const ROOT::RDF::TH1DModel P_v0_def{"Pv0","Momentum Sum;P [GeV];counts",600,0.7,2.7672};
    const ROOT::RDF::TH1DModel N_v0_def{"Nv0","Number of V0 Vertexes;N;counts",30,-0.5,29.5};
    ROOT::RDF::RResultPtr<TH1D> N_v0;
    ROOT::RDF::RResultPtr<TH1D> M_inv;
    ROOT::RDF::RResultPtr<TH1D> P_v0;

public:
    FilterStore(): node(RDataFrame(0)){};
    FilterStore(ROOT::RDF::RNode df, const std::string_view filter_name="", const std::string_view v0_status_name=""):
            v0_status_column(v0_status_name), name(filter_name), node(df) {
        if (name.empty()) {
            name = node.GetFilterNames().back();
            std::cout << "No name specified. Setting to: " << name << std::endl;
        }
        Fill();
    }

    static ULong64_t count_true(RVec<bool> v){
        /// Count the number of times v has a true entry, and return result.
        /// Use: auto h = df.Define("ntrue",count_true,"v0_valid_test").Histo1D("ntrue") // Or Sum("ntrue")
        return (ULong64_t)std::count(v.begin(),v.end(),true);
    };

    static ULong64_t count_false(RVec<bool> v){
        /// Count the number of times v has a false entry, and return result. See also count_true.
        return (ULong64_t)std::count(v.begin(),v.end(),true);
    };

    static RVec<double> copy_if_true(RVec<bool> v,RVec<double> m){
        /// Copy vector m to the output if the corresponding entry in v is true.
        /// Example:
        /// auto tmp = df.Define("minv",copy_if_true,{"v0_status_good","v0_mass"})
        /// auto h_minv = tmp.Histo1D("minv")
        /// auto h_nv0 = tmp.Define("nv0","minv.size()").Histo1D("nv0")
        /// auto h_nv0_again = df.Define("nv0",count_true,"v0_status_good").Histo1D("nv0") // same result as h_nv0
        RVec<double> out;
        for(int i=0;i<v.size();++i)
            if( v[i] ) out.push_back(m[i]);
        return(out);
    };

    void Fill(){
        /// Fill the histograms and counts.
        if( v0_status_column.empty() ){
            Count = node.Count();
            auto tmp_node = node.Define("nv0","(ULong64_t)v0_mass.size()");
            N_v0  = tmp_node.Histo1D(N_v0_def,"nv0");
            V0_Count = tmp_node.Sum<ULong64_t>("nv0");
            M_inv = node.Histo1D(M_inv_def,"v0_mass");
            P_v0 = node.Histo1D(P_v0_def,"v0_p");
        }else {
            auto tmp_node = node.Define("nv0",count_true,{v0_status_column});
            Count = tmp_node.Filter("nv0>0").Count();
            V0_Count = tmp_node.Sum<ULong64_t>("nv0");
            N_v0 = tmp_node.Histo1D(N_v0_def,"nv0");
            M_inv = node.Define("minv",copy_if_true,{v0_status_column,"v0_mass"}).Histo1D(M_inv_def,"minv");
            P_v0 = node.Define("vp",copy_if_true,{v0_status_column,"v0_p"}).Histo1D(P_v0_def,"vp");
        }

    };
};

#endif //BUMPHUNT2016_FILTERSTORE_H

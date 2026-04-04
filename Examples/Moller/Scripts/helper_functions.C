//
// Created by Maurik Holtrop on 7/28/20.
//
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"

using namespace ROOT;

void helper_functions(){
    cout << "Loaded helper functions.\n";
}

RVec<int> get_ix( RVec<int> v0_p, RVec<int> part_ecal, RVec<int> ecal_seed_idx, RVec<int> ecal_hit_index_x){
    RVec<int> out;
    for(unsigned long i=0; i< v0_p.size(); ++i){
        int i_part = v0_p[i];
        if(i_part < (int)part_ecal.size() ) {
            int i_clus = part_ecal[i_part];
            if(i_clus < (int)ecal_seed_idx.size() ) {
                int i_hit = ecal_seed_idx[i_clus];
                if(i_hit < (int)ecal_hit_index_x.size() ){
                    int ix = ecal_hit_index_x[i_hit];
                    out.push_back( ix);
                }else{
                    cout << "Mis-index for ecal_hit_index_x: i_hit = " << i_hit << endl;
                    out.push_back(-25);
                }
            }else{
                cout << "Mis-index for ecal_seed_idx : i_clus = " << i_clus << endl;
                out.push_back(-26);
            }
        }else{
            cout << "Mis-index for part_ecal, i_part = " << i_part << endl;
            out.push_back(-27);
        }
    }
    return out;
}

template<typename T>
RVec<T> diff( RVec<T> in1, RVec<T> in2){
    RVec<T> out;
    for(unsigned long i=0; i<in1.size(); ++i){
        out.push_back(in1[i]-in2[i]);
    }
    return out;
}
//
// Created by Maurik Holtrop on 8/10/20.
//

#ifndef HPS_ANALYSIS_GAMMAMIXER_H
#define HPS_ANALYSIS_GAMMAMIXER_H

#include "MiniDst.h"

class GammaMixer : public MiniDst {
public:
    explicit GammaMixer(const string &input_file=""){
        if(!input_file.empty()) input_file_list.push_back(input_file);
    };
    explicit GammaMixer(const vector<string> &infile_list): input_file_list(infile_list){};
    ~GammaMixer() override = default;

    void Start() override;
    long Run(int max_event=0) override;
    bool Find_Good_Photon_Pair(MiniDst &event, int &found1, int &found2);
    bool Write_Mixed_Photon_Events(MiniDst &event1, int e1_gamma,
                                   MiniDst &event2, int e2_gamma);
    void SetDebugLevel(const int level) override { md_Debug = level; event1.md_Debug=level; event2.md_Debug=level;};
public:

    vector<string> input_file_list{};
    int  mix_multiplyer{100};
    TTree  *tree1;
    MiniDst event1;
    TTree  *tree2;
    MiniDst event2;
    MiniDst event_out;
    int event_number2;

    double ecal_cluster_sum_max{2.5};
};


#endif //HPS_ANALYSIS_GAMMAMIXER_H

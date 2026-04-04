//
// Created by Maurik Holtrop on 8/10/20.
//

#ifndef HPS_ANALYSIS_MOLLERMIXER_H
#define HPS_ANALYSIS_MOLLERMIXER_H

#include "MiniDst.h"
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"

class MollerMixer : public MiniDst {
public:
    explicit MollerMixer(const string &input_file=""): MiniDst("mixed_moller.root"){
        if(!input_file.empty()) input_file_list.push_back(input_file);
    };
    explicit MollerMixer(const vector<string> &infile_list): input_file_list(infile_list), MiniDst("mixed_moller.root"){};
    ~MollerMixer() override = default;

    void Start() override;
    long Run(int max_event=0) override;
    bool Find_Good_Moller_Pair(MiniDst &event1, MiniDst &event2, int &found1, int &found2);
    bool Write_Mixed_Moller_Events(MiniDst &event1, int e1_e,
                                   MiniDst &event2, int e2_e);
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

    double one_electron_p_max{2.};      /// Maximum momentum for good electron. Cut out FEE.
    double two_electron_psum_min{1.9};  /// Minimum Psum for good pair.
    double two_electron_psum_max{2.6};  /// Maximum Psum for good pair.

    double two_electron_psum_mixed;      /// Energy sum of mixed photon pair.
    double two_electron_mass_mixed;    /// Invariant mass of mixed photon pair.

    int n_evt1;     /// Number of times evt1 is used in output events.
    int n_evt2;     /// Number of times evt2 is used in output events.

};


#endif //HPS_ANALYSIS_MOLLERMIXER_H

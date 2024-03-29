//
// Created by Maurik Holtrop on 8/10/20.
//

#ifndef HPS_ANALYSIS_GAMMAMIXER_H
#define HPS_ANALYSIS_GAMMAMIXER_H

#include "MiniDst.h"
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"

class GammaMixer : public MiniDst {
public:
    explicit GammaMixer(const string &input_file=""): MiniDst("mixed_photons.root"){
        if(!input_file.empty()) input_file_list.push_back(input_file);
    };
    explicit GammaMixer(const vector<string> &infile_list): input_file_list(infile_list), MiniDst("mixed_photons.root"){};
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

    double two_photon_esum_original1;  /// Energy sum of original photon pair1.
    double two_photon_esum_original2;  /// Energy sum of original photon pair2.
    double two_photon_esum_mixed;      /// Energy sum of mixed photon pair.
    double two_photon_mass_original1;  /// Invariant mass of original photon pair1.
    double two_photon_mass_original2;  /// Invariant mass of original photon pair2.
    double two_photon_mass_mixed;    /// Invariant mass of mixed photon pair.
    int n_evt1;     /// Number of times evt1 is used in output events.
    int n_evt2;     /// Number of times evt2 is used in output events.

    double ecal_cluster_sum_max{2.5};
    double ecal_cluster_delta_t_max{2.};
    double delta_esum_tolerance{0.1};
};


#endif //HPS_ANALYSIS_GAMMAMIXER_H

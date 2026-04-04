//
// Created by Maurik Holtrop on 7/17/20.
//

#include "BumpHunt2016.h"

//
// The values below come from the histogram h_cl_dtMeans1, from the file ECalTimeCorrections.root
// i.e.
// curl https://raw.githubusercontent.com/rafopar/BumpHunt_2016/master/pass4/codes/ECalTimeCorrections.root --output EcalTimeCorrections.root
// Then process the file with:
//
//    TFile f("EcalTimeCorrections.root");
//    auto hh = f.Get<TH2D>("h_cl_dtMeans1");
//    vector<double> Ecal_Time_Correction_values;
//    for(int i= -23; i<=23; ++i){
//        for(int j= -5; j<=5 ; ++j) {
//            int bin = hh->FindBin(i, j);
//            double corr = hh->GetBinContent(bin);
//            Ecal_Time_Correction_values.push_back(corr);
//        }
//    };
//    Ecal_Time_Correction_values
//
std::vector<double> BumpHunt2016::Ecal_Time_Correction_values{
-0.75026002, -0.68307181, -0.59611512, -0.64478648, -0.56925647, 0.0000000, -0.52642221, -0.53895966,
-0.68219937, -0.62116834, -0.70899636, -0.66042258, -0.52627154, -0.50428485, -0.43001173, -0.37586268,
0.0000000, -0.46692702, -0.49700639, -0.61713363, -0.75655868, -0.71113362, -0.58005049, -0.46624460,
-0.47811817, -0.38888681, -0.45690412, 0.0000000, -0.49818310, -0.51551135, -0.51935191, -0.54776497,
-0.52833296, -0.55364188, -0.44020204, -0.37748451, -0.42402051, -0.48374532, 0.0000000, -0.50064093,
-0.41440679, -0.49889730, -0.55616775, -0.55188745, -0.52250515, -0.47751180, -0.35538342, -0.34443177,
-0.28772372, 0.0000000, -0.33731854, -0.36951066, -0.24644737, -0.49418063, -0.46424499, -0.52822127,
-0.43006187, -0.38724883, -0.31373552, -2.3720402, 0.0000000, -0.40992466, -0.37391507, -0.36787960,
-0.52615783, -0.42051336, -0.38255943, -0.41243702, -0.28052009, -0.24989208, -0.28082649, 0.0000000,
-0.43903420, -0.28419604, -0.33742880, -0.33175124, -0.40711425, 0.37842711, -0.34413345, -0.20736298,
-0.20045778, -0.15312953, 0.0000000, -0.24630717, -0.17338650, -0.29094550, -0.33588279, -0.43919996,
-0.32151398, -0.13251877, -0.18678496, -0.16914229, -0.19443224, 0.0000000, -0.17779541, -0.22003293,
-0.20575450, -0.25986622, -0.33046959, -0.25271797, -0.15543545, -0.14521106, -0.086179860, -0.087964705,
0.0000000, -0.13418359, -0.12478361, -0.15542092, -0.17324335, -0.31610667, -0.22460311, -0.13417091,
-0.041118922, 0.050629063, -0.037567849, 0.0000000, -0.10422521, -0.017451536, -0.13891715, -0.17785429,
-0.22952653, -0.11467069, 0.034367627, 0.025715598, 0.0047735034, 0.0020989856, 0.0000000, -0.018882476,
0.052054538, -0.027546697, -0.076399580, -0.069685099, -0.17449679, 0.081496940, -0.11959442, 0.20397207,
0.10549869, 0.0000000, 0.062282264, 0.083412750, 0.066474813, -0.17639195, -0.14378892, -0.011365322,
0.10320838, 0.17158459, 0.15117064, 0.0000000, 0.0000000, 0.0000000, 0.068950094, 0.072955139,
-0.0077184647, -0.066690920, 0.079038826, 0.24400419, 0.25637221, 0.21495079, 0.0000000, 0.0000000,
0.0000000, 0.12113915, 0.19530752, 0.13414216, -0.095208788, 0.0081697113, 0.16052895, 0.40679956,
0.37526563, 0.0000000, 0.0000000, 0.0000000, 0.18161875, 0.16750495, 0.21482197, 0.0033187753,
0.044005563, 0.19116058, 0.24408135, 0.24724697, 0.0000000, 0.0000000, 0.0000000, 0.23623452,
0.20320798, 0.23529795, 0.070705029, -0.026638946, 0.15878922, -1.8813170, 0.33735747, 0.0000000,
0.0000000, 0.0000000, 0.24239273, 0.17266395, 0.048263587, 0.059984458, -0.15585617, 0.084470771,
0.10631026, 0.27691884, 0.0000000, 0.0000000, 0.0000000, 0.20695881, 0.21904097, 0.066242115,
0.040668700, -0.17022700, -0.0023344546, 0.070298693, 0.25747002, 0.0000000, 0.0000000, 0.0000000,
0.19578895, -0.067534133, 0.018941531, -0.16726274, -0.24794994, -0.18128575, -0.020565481, 0.12229376,
0.0000000, 0.0000000, 0.0000000, 0.17811920, -0.12236778, -0.21328774, -0.21863784, -0.17919689,
-0.22814638, 0.035966078, 0.12182934, 0.0000000, 0.0000000, 0.0000000, 0.025374759, -0.073682494,
-0.24766737, -0.36753158, -0.19467036, -0.11713873, 0.056021890, 0.17381649, -0.033793395, 0.0000000,
0.029397404, 0.046598848, -0.020761068, -0.14671359, -0.24689577, 0.0000000, 0.0000000, 0.0000000,
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
-0.23109031, -0.17656154, -0.097064838, 0.057874013, 0.074590542, 0.0000000, 0.086474275, 0.077808233,
0.0063331243, -0.13402992, -0.27318083, -0.31746407, -0.10547297, -0.0035540375, 0.17858878, 0.15726234,
0.0000000, 0.11720833, 0.065230428, 0.12783665, -0.15444367, -0.26830053, -0.27201492, -0.22048321,
-0.076012360, 0.059289663, 0.13631512, 0.0000000, 0.10446771, 0.026587143, -0.067188883, -0.14862435,
-0.21732551, -0.28516616, -0.13621088, 0.044318176, 0.087440834, 0.16096904, 0.0000000, 0.097479749,
0.060152523, -0.060386873, -0.19403209, -0.28959813, -0.30110422, -0.17647879, -0.041553816, 0.074485467,
0.078234246, 0.0000000, 0.10085945, 0.026360611, -0.032304605, -0.22026764, -0.30366253, -0.23042522,
-0.19244052, -0.022979824, 0.011346341, 0.017019872, 0.0000000, -0.14494251, 0.035052863, -0.079114704,
-0.16077045, -0.26824599, -0.36470480, -0.16517379, -0.10017535, 0.012827673, 0.052277784, 0.0000000,
0.035425544, 0.0095804258, -0.11456773, -0.17034741, -0.25600002, -0.26883126, -0.12311157, -0.070118603,
0.0047733377, -0.011481393, 0.0000000, 0.035902616, -0.0047106854, -0.14008511, -0.25384297, -0.24125654,
-0.21226010, -0.23259435, -0.20524395, -0.022366865, 0.049753875, 0.0000000, 0.024173285, -0.038587459,
-0.098852941, -0.19078475, -0.33332038, -0.28802034, -0.19208771, -0.11704047, 0.0081070349, 0.062783717,
0.0000000, 0.43136617, -0.42787377, -0.11561982, -0.17796796, -0.30722452, -0.42524328, -0.22287108,
-0.16715438, -0.050331699, -0.089436565, 0.0000000, -0.022914526, -0.055531631, -0.14672225, -0.18012780,
-0.31073504, -0.29674560, -0.19767986, -0.10317627, -0.047677983, -0.028249768, 0.0000000, -0.088072622,
-0.080666068, -0.15225426, -0.22288246, -0.23568127, -0.32840664, -0.16793883, -0.25937249, -0.089759596,
-0.080627606, 0.0000000, -0.090053177, -0.13021164, -0.12667625, -0.30886168, -0.34996117, -0.30713922,
-0.29936515, -0.19115804, -0.11180507, -0.092170679, 0.0000000, -0.19566378, -0.14460990, -0.25148830,
-0.34321460, -0.32441979, -0.36396225, -0.25401105, -0.30506204, -0.33224990, -0.11715130, 0.0000000,
-0.26343853, -0.020136777, -0.24119883, -0.30098379, -0.28676175, -0.21753017, -0.26334399, -0.26068267,
-0.11097404, -0.28708684, 0.0000000, -0.15238664, -0.25213766, -0.17854512, -0.18860816, -0.25672936,
-0.36047962, -0.37820041, -0.18927729, -0.35255372, -0.19213386, 0.0000000, -0.16307235, -0.21663737,
-0.25116925, -0.26940815, -0.34731667, -0.49903613, -0.28458297, -0.38627261, -0.47116022, -0.28459192,
0.0000000, -0.35342644, -0.31046296, -0.38556618, -0.37519874, -0.40939309, -0.58824693, -0.56770292,
-0.29905121, -0.49216849, -0.13123238, 0.0000000, -0.36832417, -0.34828206, -0.30448403, -0.51716441,
-0.44250883, -0.69462212, -0.46341327, -0.22133034, -0.34785473, -0.28406102, 0.0000000, -0.23767368,
-0.54057508, -0.14932041, -0.19918182, -0.32803825, -0.62760300, -0.61223266, -0.54376651, -0.58133528,
-0.31691209, 0.0000000, -0.18486808, -0.096299102, -0.52916715, -0.37889366, -0.51342784, 0.0000000,
-1.0841206, -0.52346845, -0.28761965, -0.34782407, 0.0000000, -0.39253779, -0.44088768, -0.15958038,
-0.82053283, -0.54027944, 0.0000000, -0.71669285, -1.0184029, -0.46685955, -0.0028577400, 0.0000000,
0.13244554, -0.082701878, -0.58102489, -0.35462361, 0.0000000 };


BumpHunt2016::BumpHunt2016(std::string files) {
    if (!files.empty()) {
        ch.Add(files.c_str());
    } else {
        return;
    }
}

RNode BumpHunt2016::Select_Pair1(RNode in) {
    ///
    /// First filter step that prepares the most basic parts of the data.
    /// Selects the pair1 trigger events.
    ///   - Pair1 is bit9 in the 2019 data, and the MiniDST maker translates the 2016 bits to the 2019 definition.
    /// Computes the total momentum for the v0 particles.
    ///

    // You can also use: .Filter("(trigger & 0x200)","pair1")
    return FilterStep(in.Define("v0_p",Vector_magnitude,{"v0_px","v0_py","v0_pz"})
    .Filter([](unsigned int trig){return (trig & (1<<9));}, {"trigger"}, "pair1"));
};

RNode BumpHunt2016::Correct_v0_clus_times(RNode in){
    ///
    /// Computes the corrected times for the ECal, according to Rafo's corrections.
    /// The corrections are stored in a 2-D histogram, by crystal index.
    ///
    /// This step only defines a new variable, without any filtering, so there is no need for a FilterStep.
    ///

    // std::vector<double>  &tmp_cp = BumpHunt2016_Ecal_Time_Correction_values;
    auto ecal_time_correction = [this](RVec<double> ecal_time, RVec<int> ecal_ix, RVec<int> ecal_iy){
        RVec<double> out;
        for(int i = 0; i< ecal_time.size(); ++i){
            int bin = this->Ecal_Time_Corr_bin(ecal_ix[i],ecal_iy[i]);
            double correction = this->Ecal_Time_Correction_values[bin];
            double newtime  = ecal_time[i] - correction;
            out.push_back(newtime);
        }
        return out;
    };

    return FilterStep(in.Define("v0_em_clus_time_corr",ecal_time_correction,{"v0_em_clus_time", "v0_em_clus_ix", "v0_em_clus_iy"}).
    Define("v0_ep_clus_time_corr",ecal_time_correction,{"v0_ep_clus_time", "v0_ep_clus_ix", "v0_ep_clus_iy"}),
    "clus_time_corr");
}
//
//RNode BumpHunt2016::Correct_all_ecal_times(RNode in){
//    ///
//    /// Computes the corrected times for the ECal, according to Rafo's corrections.
//    /// The corrections are stored in a 2-D histogram, by crystal index.
//    ///
//    /// This step only defines a new variable, without any filtering, so there is no need for a FilterStep.
//    ///
//    auto ecal_time_correction = [this](RVec<double> ecal_time, RVec<int> ecal_ix, RVec<int> ecal_iy){
//        RVec<double> out;
//        for(int i = 0; i< ecal_time.size(); ++i){
//            out.push_back( ecal_time[i] - this->Ecal_Time_Correction_values[this->Ecal_Time_Corr_bin(ecal_ix[i],ecal_iy[i])]);
//        }
//        return out;
//    };
//    return in.Define("ecal_cluster_time_corr",ecal_time_correction,{"ecal_cluster_time", "ecal_cluster_seed_ix", "ecal_cluster_seed_iy"});
//}
//

RNode BumpHunt2016::Cut_on_reduced_chi1(RNode in, double em_chi2_max, double ep_chi2_max){
    ///
    /// Filter on the vertex type.
    ///
    /// First define lambda function to calculate the rchi2 and then one to cut on Chi2 for electron and positron.
    /// Note: This needs to be the reduced Chi2, which is different for 5 hit and 6 hit tracks. This requires the
    /// v0_em_n_hits and v0_ep_n_hits, which come from the associated tracks. NDF = 2*n_hits - 5;
    ///
    auto add_rchi = [](RVec<double> v0_ex_chi, RVec<int> v0_ex_track, RVec<int> track_n_hits)->RVec<double>{
        RVec<double> out;
        for(int i=0;i< v0_ex_track.size(); ++i){
            int n_hits = track_n_hits[v0_ex_track[i]];
            int NDF = 2*n_hits - 5;
            double rchi = v0_ex_chi[i]/NDF;
            out.push_back(rchi);
        }
        return out;
    };

    auto do_chi2_cut = [em_chi2_max, ep_chi2_max](RVec<double> v0_em_rchi2, RVec<double> v0_ep_rchi2){
        RVec<bool> out;
        for(int i=0;i<v0_em_rchi2.size();++i){
            out.push_back((v0_em_rchi2[i] < em_chi2_max ) && (v0_ep_rchi2[i] < ep_chi2_max) );
        }
        return out;
    };

    // Define steps can be removed if items are added to the MiniDST writer.
    return FilterStep(
            in.Define("v0_em_rchi2",add_rchi,{"v0_em_chi2", "v0_em_track","track_n_hits"}).
            Define("v0_ep_rchi2",add_rchi,{"v0_ep_chi2", "v0_ep_track","track_n_hits"}).
            Define("v0_valid_chi2", do_chi2_cut,{"v0_em_rchi2", "v0_ep_rchi2"}),
            "chi2_cut","v0_valid_chi2");
}

RNode BumpHunt2016::Cut_on_good_pid(RNode in, double good_pid_max){
    ///
    /// Cut on the "goodness of pid" which relates to how closely the track points to the cluster.
    ///
    auto do_good_pid_cut = [](RVec<double> v0_em_good_pid, RVec<double> v0_ep_good_pid){
        RVec<bool> out;
        for(int i=0; i<v0_em_good_pid.size(); ++i){
            out.push_back( (v0_em_good_pid[i]< 10) && v0_ep_good_pid[i] < 10);
        }
        return out;
    };
    return FilterStep(in.Define("v0_valid_good_pid", do_good_pid_cut,
                                                {"v0_em_good_pid", "v0_ep_good_pid"}),
                                   "good_pid","v0_valid_good_pid");

}

RNode BumpHunt2016::Simple_track_cluster_time_cut(RNode in, double em_dt_max, double ep_dt_max) {
    ///
    /// Do a simple time cut: | cluster_time - offset - track_time | < dt_max
    ///
    double offset = clust_track_time_offset; // So we can capture value.
    auto do_time_cut = [offset, em_dt_max, ep_dt_max](RVec<double> v0_em_track_time, RVec<double> v0_em_clus_time,
                              RVec<double> v0_ep_track_time, RVec<double> v0_ep_clus_time,
                              RVec<bool> v0_valid_good_pid) {
        RVec<bool> out;
        for (int i = 0; i < v0_em_track_time.size(); ++i) {
            out.push_back((abs(v0_em_clus_time[i] - offset - v0_em_track_time[i]) < em_dt_max) &&
                          (abs(v0_ep_clus_time[i] - offset - v0_ep_track_time[i]) < ep_dt_max ) &&
                          v0_valid_good_pid[i]);
        }
        return out;
    };
    return FilterStep( in.Define("v0_valid_time", do_time_cut,
            {"v0_em_track_time", "v0_em_clus_time", "v0_ep_track_time", "v0_ep_clus_time"}),
                     "valid_time", "v0_valid_time");
}

RNode BumpHunt2016::Momentum_dependent_track_cluster_time_cut(RNode in){
///
/// Cut on the track cluster time difference with limits that depend on the particle momentum.
/// This is Rafo's cut, as described in the 2016 Bump Hunt Analysis Note.
///
/// Constants for the polynomials come from the top level class, so "this" is captured.
///
    auto do_time_cut = [this](RVec<double> track_time, RVec<double> clus_time, RVec<double> pos_ecal_y, RVec<double> mom){
        RVec<bool> out;
        std::vector<double> c_low;
        std::vector<double> c_hi;
        for(int i=0; i<track_time.size(); ++i){  // For each V0
            double dt = abs(clus_time[i] - this->clust_track_time_offset - track_time[i]);
            double p = mom[i];
            if( pos_ecal_y[i]>0){
                c_low = this->Tr_Cl_Time_Cut_Top_Lower_Lim;
                c_hi  = this->Tr_Cl_Time_Cut_Top_Upper_Lim;
            }else{
                c_low = this->Tr_Cl_Time_Cut_Bot_Lower_Lim;
                c_hi  = this->Tr_Cl_Time_Cut_Bot_Upper_Lim;
            }
            double low_lim = (c_low[0] + p*(c_low[1] + p*(c_low[2] + p*c_low[3])));
            double hi_lim = (c_hi[0] + p*(c_hi[1] + p*(c_hi[2] + p*c_hi[3])));
            out.push_back( low_lim < dt < hi_lim );
        }
        return out;
    };
    auto combine_cuts = [](RVec<bool> cut1, RVec<bool> cut2, RVec<bool> passthrough){
        RVec<bool> out;
        for(int i=0; i<cut1.size(); ++i){
            out.push_back( cut1[i] && cut2[i] && passthrough[i]);
        }
        return(out);
    };

    // The time cut is applied in 3 steps: the em time, the ep time, and the combination.
    return FilterStep( in.Define("v0_em_valid_time",do_time_cut,{"v0_em_track_time","v0_em_clus_time", "v0_em_pos_ecal_y", "v0_em_p"}).
    Define("v0_ep_valid_time",do_time_cut,{"v0_ep_track_time","v0_ep_clus_time", "v0_ep_pos_ecal_y", "v0_ep_p"}).
    Define("v0_valid_time",combine_cuts,{"v0_em_valid_time","v0_ep_valid_time","v0_valid_good_pid"}),
    "valid_time","v0_valid_time" );
}

RNode  BumpHunt2016::Cut_on_min_max_vertex_momentum(RNode in, double psum_min, double psum_max){
    ///
    /// Cut on the minimum and maximum momentum for the vertexed particles.
    /// This is the equivalent of cutting on the magnitude of the 3-vector sum of the electron and positron momentum.
    ///
    auto cut_v0_mom = [psum_min, psum_max](RVec<double> px, RVec<double> py, RVec<double> pz)->RVec<bool>{
        RVec<bool> out;
        for(int i=0; i< px.size(); ++i){
            double p = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            out.push_back( psum_min < p && p < psum_max);
        }
        return out;
    };

    return FilterStep( in.Define("v0_valid_pv",cut_v0_mom,{"v0_px","v0_py","v0_pz"}),
            "valid_pv","v0_valid_pv");
}

RNode BumpHunt2016::Cut_on_cluster_in_time(RNode in){
    ///
    /// Cut on the time of the Ecal cluster to be in a specific window.
    /// For top the window is just 30 - 70, for bottom it is an energy dependent cut, documented in the analysis note.
    ///
    auto cut_on_cluster_time = [this](RVec<double> cluster_time, RVec<double> cluster_energy, RVec<double> cluster_pos_y)->RVec<bool>{
        RVec<bool> out;
        double low=0;
        double hi =0;
        for(int i=0;i< cluster_time.size(); ++i){
            if(cluster_pos_y[i] > 0){
                low = 30.;
                hi  = 70.;
            }else{
                double e = cluster_energy[i];
                low = this->Cluster_Time_Cut_Bot_Lower_Lim[0] +
                        e*(this->Cluster_Time_Cut_Bot_Lower_Lim[1] +
                        e*(this->Cluster_Time_Cut_Bot_Lower_Lim[2]));
                hi = this->Cluster_Time_Cut_Bot_Upper_Lim[0] +
                      e*(this->Cluster_Time_Cut_Bot_Upper_Lim[1] +
                         e*(this->Cluster_Time_Cut_Bot_Upper_Lim[2]));
            }
            out.push_back( low < cluster_time[i] && cluster_time[i] < hi);
        }
        return out;
    };

    return FilterStep(in.Define("v0_em_valid_clus_time", cut_on_cluster_time,
            {"v0_em_clus_time","v0_em_clus_energy","v0_em_clus_pos_y"})
            .Define("v0_ep_valid_clus_time", cut_on_cluster_time,
                    {"v0_ep_clus_time","v0_ep_clus_energy","v0_ep_clus_pos_y"})
                    .Define("v0_valid_clus_time",Combine_cuts_helper2,
                            {"v0_em_valid_clus_time","v0_ep_valid_clus_time"}),
                    "valid_clus_time","v0_valid_clus_time");
};

RNode BumpHunt2016::Cut_on_cluster_time_diff(RNode in, double dt_max) {
    ///
    /// Cut on the time difference between the two clusters.
    ///
    auto cut_dt_clus = [dt_max](RVec<double> v0_em_clus_time, RVec<double> v0_ep_clus_time)->RVec<bool>{
        RVec<bool> out;
        for(int i=0; i< v0_em_clus_time.size(); ++i){
            out.push_back( abs(v0_em_clus_time[i] - v0_ep_clus_time[i]) < dt_max );
        }
        return out;
    };

    return FilterStep(in.Define("v0_valid_clus_dt",cut_dt_clus,{"v0_em_clus_time","v0_ep_clus_time"}),
            "valid_clus_dt","v0_valid_clus_dt");
}

RNode BumpHunt2016::Combine_valid_v0_and_reduce(RNode in, std::vector<std::string> args){
    ///
    /// Take all the cuts on v0 and combine the logic into a "v0_valid", then remove events where
    /// all the v0_valid are false.
    ///
    auto fin_cut = [](RVec<bool> v){return std::any_of(v.begin(),v.end(),[](bool x){ return x;}); };
    return FilterStep(Combine_cuts(in,"v0_valid",args).Filter(fin_cut,{"v0_valid"},"final_selection"),
               "combine_cuts","v0_valid");
}

void BumpHunt2016::Setup() {
    ///
    /// Setup a standard analysis.
    ///
    /// The analysis here is the same as the 2016 Bumphunt as documented in the analysis note.
    ///
    std::cout << "Setting up version 2.1 \n";
     // Step 1 - Select only events with the pair1 trigger == true.
    auto pair1 = Select_Pair1(dataframe);
    auto ecal_time_corr = Correct_v0_clus_times(pair1);
    auto valid_clus_times = Cut_on_cluster_in_time(ecal_time_corr);
    auto valid_pv = Cut_on_min_max_vertex_momentum(valid_clus_times, 1.9, 2.4);
    auto valid_clus_time_diff = Cut_on_cluster_time_diff(valid_pv, 1.43);
    auto combo = Combine_valid_v0_and_reduce(valid_clus_time_diff,{"v0_valid_clus_time","v0_valid_pv","v0_valid_clus_dt"});
    auto final = FilterStep(combo.Define("n_v0",FilterStore::count_true,{"v0_valid"})
            .Filter("n_v0==1"),"final","v0_valid");
}

void BumpHunt2016::Process() {
    ///
    /// Trigger Processing of all the filter steps and print a summary.
    ///
    if( filters.size() == 0){
        std::cout << "Please run Setup() first.\n";
        return;
    }
    std::cout << "Start processing.\n";
    TStopwatch timer;
    timer.Start();
    Print();
    timer.Stop();
    std::cout << "Process time was " << timer.RealTime() << " s  CPU: " << timer.CpuTime() << " s \n";
}

void BumpHunt2016::Print(Option_t *opt) {
    unsigned long long top_count = *(filters[0].Count);
    unsigned long long top_v0_count = *(filters[0].V0_Count);
    unsigned long long bot_count = *(filters.back().Count);
    unsigned long long bot_v0_count = *(filters.back().V0_Count);

    std::cout << "Filter summary:  step (  v0 valid column )       events           v0 vertexes \n";
    for (auto f : filters) {
        printf(" %20s (%20s) = %9llu (%8.4f %%)  %9llu  (%8.4f %%)\n", f.name.c_str(), f.v0_status_column.c_str(),
               (*f.Count), 100. * double(*f.Count) / top_count, (*f.V0_Count),
               100. * double(*f.V0_Count) / top_v0_count);
    }
    std::cout << "Final counts: " << bot_count << " events, and " << bot_v0_count << " vertexes.\n";
}
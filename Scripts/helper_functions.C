//
// Created by Maurik Holtrop on 7/28/20.
//
#include "TStyle.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include <locale.h>

using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;

void SetStyle(void);

void helper_functions(){
    cout << "Loaded helper functions.\n";
    setlocale(LC_NUMERIC, "");
    SetStyle();
    ROOT::EnableImplicitMT();
}

RVec<double> compute_svt_raw_hit_delta_t(RVec<double> svt_hit_time,
                                         RVec<double> svt_raw_hit_t0,
                                         RVec<RVec<int>> svt_hit_raw_index){
    /// Computes the difference in time between the raw hit and the svt (composite) hit.
    RVec<double> out;
    for(int i_hit = 0; i_hit < (int)svt_hit_raw_index.size(); i_hit++){
        for(int i_raw=0; i_raw < (int)svt_hit_raw_index[i_hit].size(); i_raw++){
            // printf("%3d:%1d  Tsvt= %6.3f  Traw = %6.3f DeltaT = %6.3f \n",
            // i_hit,i_raw,md->svt_hit_time[i_hit],
            // md->svt_raw_hit_t0[md->svt_hit_raw_index[i_hit][i_raw]],
            // md->svt_raw_hit_t0[md->svt_hit_raw_index[i_hit][i_raw]]-md->svt_hit_time[i_hit]);
            double dtime = svt_raw_hit_t0[svt_hit_raw_index[i_hit][i_raw]] - svt_hit_time[i_hit];
            out.push_back(dtime);
        }
    }
    return(out);
}

RVec<double> compute_svt_hit_time_ext(RVec<double> svt_hit_time,
                                      RVec<RVec<int>> svt_hit_raw_index){
    /// Extends the svt_hit_time array to be the same size as the svt_raw_hit_delta array.
    /// This permits a 2d histogram of dt vs svt_time.
    RVec<double> out;
    for(int i_hit = 0; i_hit < (int)svt_hit_raw_index.size(); i_hit++){
        for(int i_raw=0; i_raw < (int)svt_hit_raw_index[i_hit].size(); i_raw++){
            out.push_back(svt_hit_time[i_hit]);
        }
    }
    return(out);
}

struct Compute_Track_SVT_dt {

    int type_to_use{1};
    RVec<double> operator()(RVec<double> track_time,
                            RVec<int> track_type,
                            RVec<RVec<int>> track_svt_hits,
                            RVec<double> svt_hit_time) {
        /// Compute an array of time differences between the track time and the times of the
        /// individual SVT hits that make up the track. Use 'type_to_use' to select track type.
        RVec<double> out;
        for (int i = 0; i < (int) track_time.size(); ++i) {
            if (type_to_use >= 0 && track_type[i] == type_to_use) {
                for (int j = 0; j < (int) track_svt_hits[i].size(); ++j) {
                    int hit_index = track_svt_hits[i][j];
                    double dtime = svt_hit_time[hit_index] - track_time[i];
                    out.push_back(dtime);
                }
            }
        }
        return (out);
    }
};

template <typename T> RVec<T> diff_vec( RVec<T> in1, RVec<T> in2){
    RVec<T> out;
    for(unsigned long i = 0; i < in1.size(); ++i){
        out.push_back(in1[i] - in2[i]);
    }
    return out;
}

int    select_int_value{0};
template <typename T>
RVec<T> select_on_int( RVec<int> sel, RVec<T> var ){
    RVec<T> out;
    for(unsigned long i = 0; i < sel.size(); ++i) {
        if (sel[i] == select_int_value) {
            out.push_back(var[i]);
        }
    }
    return out;
}

double select_less_than_value{0.1};
RVec<double> select_less_than(RVec<double> sel, RVec<double> var){
    RVec<double> out;
    for(unsigned long i = 0; i < sel.size(); ++i) {
        if (sel[i] < select_less_than_value) {
            out.push_back(var[i]);
        }
    }
    return out;
}

RVec<double> add_part_chi2(RVec<double> track_chi, RVec<int> part_track){
    RVec<double> out;
    for(unsigned long i = 0; i < part_track.size(); ++i){
        if(part_track[i]>= 0) {
            out.push_back(track_chi[part_track[i]]);
        }else{
            out.push_back(-3.);
        }
    }
    return out;
}

RVec<double> calc_minv_1(RVec<int> selector, RVec<double> e, RVec<double> px, RVec<double> py, RVec<double> pz){
    /// The selector must be true for those particles that are to be considered for calculating the invariant mass.
    /// The selector vector must be the same size as e (and thus px, py, pz).
    /// The output is only filled if at least two such particles exist.
    ///
    /// The selector can be created using the RVec operations, i.e.
    ///       df.Define("selector","(part_pdg == 22)")
    /// will select photons.
    /// Additional conditions can be combined, or the selector can be computed in another function.
    ///
    RVec<double> out;
    for(unsigned long i=0; i< selector.size(); ++i) {
        if (selector[i] == 0) continue;
        Math::PxPyPzEVector p1(px[i], py[i], pz[i], e[i]);
        out.push_back(p1.M());
    }
    return out;
}

double calc_minv_use_particle_mass = 0;
RVec<double> calc_minv_2(RVec<int> selector, RVec<double> px, RVec<double> py, RVec<double> pz){
    /// The selector must be true for those particles that are to be considered for calculating the invariant mass.
    /// The selector vector must be the same size as e (and thus px, py, pz).
    /// The output is only filled if at least two such particles exist.
    ///
    /// The selector can be created using the RVec operations, i.e.
    ///       df.Define("selector","(part_pdg == 22)")
    /// will select photons.
    /// Additional conditions can be combined, or the selector can be computed in another function.
    ///
    RVec<double> out;
    if( selector.size()<2) return out;
    for(unsigned long i=0; i< selector.size(); ++i){
        if( selector[i] == 0) continue;
        for(unsigned long j=i+1; j< selector.size(); ++j){
            if( selector[j] == 0) continue;
            Math::PxPyPzMVector p1(px[i], py[i], pz[i], calc_minv_use_particle_mass);
            Math::PxPyPzMVector p2(px[j], py[j], pz[j], calc_minv_use_particle_mass);
            auto pf = p1 + p2;
            out.push_back( pf.mass());
        }
    }
    return out;
}

double calc_minv_cluster_time_cut = 2.;
double calc_minv_cluster_energy_min_cut = 0.3;
double calc_minv_cluster_energy_max_cut = 1.7;
double calc_minv_cluster_energy_sum_cut = 2.2;
RVec<double> calc_minv_2gt(RVec<int> pdg, RVec<int> ecal_cluster, RVec<double> cluster_time,
                           RVec<double> cluster_energy, RVec<double> px, RVec<double> py, RVec<double> pz){
     /// same as calc_minv_2, but adds a cut on the ecal cluster time and selects for photons without a selector.
     /// These additional cuts could possibly be done in the selector, but it seems trickier if there are more than
     /// two clusters in the event, with different time combinations.
    RVec<double> out;
    if( pdg.size()<2) return out;
    for(unsigned long i=0; i< pdg.size(); ++i){
        if( pdg[i] != 22) continue;
        for(unsigned long j=i+1; j< pdg.size(); ++j){
            if( pdg[j] != 22) continue;
            // Find the ecal clusters.
            int clus_id_i = ecal_cluster[i];
            int clus_id_j = ecal_cluster[j];
            if( clus_id_i >= (int)cluster_time.size() || clus_id_i < 0) continue; // No good cluster found.
            if( clus_id_j >= (int)cluster_time.size() || clus_id_j < 0) continue; // No good cluster found.
            double clus_time_i = cluster_time[clus_id_i];
            double clus_time_j = cluster_time[clus_id_j];
            if( abs(clus_time_i - clus_time_j) > calc_minv_cluster_time_cut ) continue;  // Out of time.
            double clus_energy_i = cluster_energy[clus_id_i];
            double clus_energy_j = cluster_energy[clus_id_j];
            if(clus_energy_i < calc_minv_cluster_energy_min_cut || clus_energy_i > calc_minv_cluster_energy_max_cut ) continue;
            if(clus_energy_j < calc_minv_cluster_energy_min_cut || clus_energy_j > calc_minv_cluster_energy_max_cut ) continue;
            if( clus_energy_i + clus_energy_j > calc_minv_cluster_energy_sum_cut) continue; // Too much energy.
            Math::PxPyPzMVector p1(px[i], py[i], pz[i], calc_minv_use_particle_mass);
            Math::PxPyPzMVector p2(px[j], py[j], pz[j], calc_minv_use_particle_mass);
            auto pf = p1 + p2;
            out.push_back( pf.mass());
        }
    }
    return out;
}



RVec<double> calc_minv_2e(RVec<int> selector, RVec<double> e, RVec<double> px, RVec<double> py, RVec<double> pz){
    /// The selector must be true for those particles that are to be considered for calculating the invariant mass.
    /// The selector vector must be the same size as e (and thus px, py, pz).
    /// The output is only filled if at least two such particles exist.
    ///
    /// The selector can be created using the RVec operations, i.e.
    ///       df.Define("selector","(part_pdg == 22)")
    /// will select photons.
    /// Additional conditions can be combined, or the selector can be computed in another function.
    ///
    RVec<double> out;
    if( selector.size()<2) return out;
    for(unsigned long i=0; i< selector.size(); ++i){
        if( selector[i] == 0) continue;
        for(unsigned long j=i+1; j< selector.size(); ++j){
            if( selector[j] == 0) continue;
            Math::PxPyPzEVector p1(px[i], py[i], pz[i], e[i]);
            Math::PxPyPzEVector p2(px[j], py[j], pz[j], e[j]);
            auto pf = p1 + p2;
            out.push_back(pf.mass());
//            double minv2 = (e[i]+e[j])*(e[i]+e[j]) - (px[i]+px[j])*(px[i]+px[j]) - (py[i]+py[j])*(py[i]+py[j]) - (pz[i]+pz[j])*(pz[i]+pz[j]);
//            if(minv2>0) out.push_back(  sqrt( minv2));
//            else        out.push_back( -sqrt(-minv2));
        }
    }
    return out;
}

RVec<double> Ecal_hit_tdif(double rf_time, RVec<double> ecal_hit_time){
    RVec<double> out;
    for(unsigned long i=0; i<ecal_hit_time.size(); ++i){
        out.push_back(ecal_hit_time[i]-rf_time);
    }
    return out;
}

RVec<double> Ecal_hit_tdif2(RVec<double> ecal_hit_time){
    RVec<double> out;
    for(unsigned long i=1; i<ecal_hit_time.size(); ++i){
        out.push_back(ecal_hit_time[i] - ecal_hit_time[0]);
    }
    return out;
}

void SetStyle(void){

    auto hpsStyle= new TStyle("HPS","HPS style");
    // use plain black on white colors
    int icol=0;
    hpsStyle->SetFrameBorderMode(icol);
    hpsStyle->SetCanvasBorderMode(icol);
    hpsStyle->SetPadBorderMode(icol);
    hpsStyle->SetPadColor(icol);
    hpsStyle->SetCanvasColor(icol);
    hpsStyle->SetStatColor(icol);
    // hpsStyle->SetFillColor(icol);

    // set the paper & margin sizes
    hpsStyle->SetPaperSize(20,26);
    hpsStyle->SetPadTopMargin(0.05);
    hpsStyle->SetPadRightMargin(0.05);
    hpsStyle->SetPadBottomMargin(0.18);
    hpsStyle->SetPadLeftMargin(0.14);

    // use large fonts
    // int font=72;
    int font=42;
    double tsize=0.08;
    double tzsize = 0.045;
    hpsStyle->SetTextFont(font);

    hpsStyle->SetTextSize(tsize);
    hpsStyle->SetLabelFont(font,"x");
    hpsStyle->SetTitleFont(font,"x");
    hpsStyle->SetLabelFont(font,"y");
    hpsStyle->SetTitleFont(font,"y");
    hpsStyle->SetLabelFont(font,"z");
    hpsStyle->SetTitleFont(font,"z");

    hpsStyle->SetLabelSize(tsize,"x");
    hpsStyle->SetTitleSize(tsize,"x");
    hpsStyle->SetLabelSize(tsize,"y");
    hpsStyle->SetTitleSize(tsize,"y");
    hpsStyle->SetLabelSize(tzsize,"z");
    hpsStyle->SetTitleSize(tzsize,"z");

    hpsStyle->SetTitleOffset(0.7,"y");
    hpsStyle->SetTitleOffset(1.15,"x");

    // use bold lines and markers
    // hpsStyle->SetMarkerStyle(20);
    hpsStyle->SetMarkerSize(1.0);
    hpsStyle->SetHistLineWidth(3);
    hpsStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

    // get rid of X error bars and y error bar caps
    // hpsStyle->SetErrorX(0.001)

    // do not display any of the standard histogram decorations
    hpsStyle->SetOptTitle(0);
    // hpsStyle->SetOptStat(1111);
    hpsStyle->SetOptStat(0);
    // hpsStyle->SetOptFit(1111);
    hpsStyle->SetOptFit(0);

    // put tick marks on top and RHS of plots
    hpsStyle->SetPadTickX(1);
    hpsStyle->SetPadTickY(1);

    gROOT->SetStyle("Plain");
    // gStyle->SetPadTickX(1);
    // gStyle->SetPadTickY(1);
    gROOT->SetStyle("HPS");
    gROOT->ForceStyle();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // overwrite hps styles
    hpsStyle->SetPadLeftMargin(0.14);
    hpsStyle->SetPadRightMargin(0.06);
    hpsStyle->SetPadBottomMargin(0.11);
    hpsStyle->SetPadTopMargin(0.05);
    hpsStyle->SetFrameFillColor(0);

    int NCont = 255;
    std::vector<double> stops{0.00, 0.34, 0.61, 0.84, 1.00};
    std::vector<double> red{0.00, 0.00, 0.87, 1.00, 0.51};
    std::vector<double> green{0.00, 0.81, 1.00, 0.20, 0.00};
    std::vector<double> blue{0.51, 1.00, 0.12, 0.00, 0.00};
    TColor::CreateGradientColorTable(stops.size(), stops.data(), red.data(), green.data(), blue.data(), NCont);
    gStyle->SetNumberContours(NCont);
}

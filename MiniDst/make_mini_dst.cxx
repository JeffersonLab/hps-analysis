//
// Created by Maurik Holtrop on 7/11/20.
//

#include "cxxopts.hpp"
#include "TChain.h"
#include "Dst2016.h"
#include "LcioReader.h"
#include <locale.h>

using namespace std;

int main(int argc, char **argv){

    setlocale(LC_NUMERIC, "");

    cxxopts::Options options(argv[0], " - Write a ROOT MiniDst for HPS data.\n");
    options
            .positional_help(" infile1 infile2 ...")
            .show_positional_help();

    bool apple = false;

    options
            .add_options()
                    ("d,debug", "Increase debug level")
                    ("q,quiet", "Run quiet.")
                    ("a,all", "Store all known values in file. Equivalent to -c -e -s",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("c,ecal_clusters", "Store Ecal Clusters",
                            cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("e,ecal_hits", "Store Ecal Hits",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("s,svt_hits", "Store SVT Hits",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("g,gbl_tracks_only", "Only store GBL tracks",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("T,no_tracks", "DO NOT Store tracks",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("M,no_mc_particles", "DO NOT Store MCParticles, even if they are in the input",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("m,mc_particles", "Store MCParticles, even if they are not in the input (so all zeros).",
                     cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("o,output", "Output file", cxxopts::value<std::string>()->default_value("mini_dst.root"))
                    ("i,inputfiles",
                     "List of input files which will be concatenated into a single output mini dst file. The -i "
                     "is optional. ",
                     cxxopts::value<std::vector<std::string>>())
                    ("n,num_evt","Only process num_evt event", cxxopts::value<long>()->default_value("0"))
                    ("h,help", "Print help")
                    ;


    options.parse_positional({"inputfiles"});
    auto args = options.parse(argc, argv);
    if (args.count("help") || (args.count("inputfiles") == 0))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    auto& infiles = args["inputfiles"].as<std::vector<std::string>>();
    auto& outfile = args["output"].as<std::string>();
    auto& store_all = args["all"].as<bool>();
    auto& store_ecal_clusters = args["ecal_clusters"].as<bool>();
    auto& store_ecal_hits = args["ecal_hits"].as<bool>();
    auto& store_svt_hits = args["svt_hits"].as<bool>();
    auto& no_tracks = args["no_tracks"].as<bool>();
    auto& gbl_tracks_only = args["gbl_tracks_only"].as<bool>();
    auto& no_mc_particles = args["no_mc_particles"].as<bool>();
    auto& store_mc_particles = args["no_mc_particles"].as<bool>();
    auto& num_evt = args["num_evt"].as<long>();
    int debug = 0;
    if(args.count("quiet") == 0 ){
        debug = args.count("debug") + 1;
    }else {
        debug = 0;
    }
    if(debug) {
        cout << "Debug was set to: " << debug << endl;
        cout << "Input: ";
        for (auto &v : infiles) {
            cout << " " << v;
        }
        cout << endl;
        cout << "Output: " << outfile << endl;
        if(store_ecal_clusters) cout << "Storing the ECAL clusters.\n";
        if(store_ecal_hits) cout << "Storing the ECAL hits.\n";
    }

    MiniDst *dst{nullptr};
    bool is_dst_type = false;
    if( infiles.size()>0 && infiles[0].find(".root") != string::npos ) {
        // The first file in the list has .root extension.
        is_dst_type = true;
        auto chain = new TChain("HPS_Event");
        for (auto &v : infiles) {
            chain->Add(v.c_str());
        }
        auto dst2016 = new Dst2016(chain);
        dst = static_cast<MiniDst*>(dst2016);
        if(dst2016->event->getNumberOfMCParticles() > 0 && !dst->write_mc_particles && !no_mc_particles){
            cout << "Warning: Input is MC Data, but write_mc_particles is not set. Turning on write_mc_particles!\n";
            dst->write_mc_particles = true;
        }
    }else if( infiles.size()>0 && infiles[0].find(".slcio") != string::npos ) {
        auto dstlcio = new LcioReader(infiles);
        dst = static_cast<MiniDst*>(dstlcio);
    }
    int debug_code = 0;
    if( debug <= 0){
        debug_code = 0;
    }else if (debug == 1){
        debug_code = 0x0A;
    }else{
        debug_code = 0x0A + ( (debug-1) << 4);
        printf("Debug code = 0x%02X \n",debug_code);
    }
    if(debug>0) cout << "Debug code = " << debug_code << endl;
    dst->SetDebugLevel(debug_code);
    dst->write_mc_particles = store_mc_particles;
    dst->write_ecal_cluster = store_ecal_clusters || store_all;
    dst->write_ecal_hits = store_ecal_hits || store_all;
    dst->write_svt_hits = store_svt_hits || store_all;
    dst->write_tracks = !no_tracks;
    dst->write_only_gbl_tracks = gbl_tracks_only;
    dst->SetOutputFileName(outfile);
#ifdef DEBUG
    cout << "Extra debug code compiled.\n";
    dst->Counter_Freq =10000;
#endif
    dst->Start();
    dst->Run(num_evt);
    dst->End();
}
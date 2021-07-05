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

    // This is a nicer way to do options in C++. See cxxopts.hpp file.
    string help_string = "Write a ROOT MiniDst for HPS data.\nVersion: 1.0.0, using MiniDst.h version "+MiniDst::_version_()+
                         "\nCompiled with "+__VERSION__+"\n";
    cxxopts::Options options(argv[0], help_string);
    options
            .positional_help(" infile1 infile2 ...")
            .show_positional_help();

    options.add_options()
            ("d,debug", "Increase debug level")
            ("q,quiet", "Run quiet.")
            ("a,all", "Store all known values in file, except the raw stuff. Equivalent to -c -e -s",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("c,ecal_clusters", "Store Ecal Clusters",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("e,ecal_hits", "Store Ecal Hits",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("s,svt_hits", "Store SVT 3D and/or strip Hits",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("r,raw_svt_hits", "Store SVT RAW Hits",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("g,gbl_tracks", "Store GBL tracks",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("G,gbl_kinks", "Store GBL track kink data",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("k,kf_tracks", "Store KF tracks",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("t,matched_tracks", "Store Matched tracks",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("T,all_tracks", "Store all tracks (set gbl, kf and matched track to true).",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("no_kf_particles","Do NOT store the KF particles or vertexes.",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("no_gbl_particles","Do NOT store the GBL particles or vertexes.",
             cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
            ("M,no_mc_particles", "DO NOT Store MCParticles, even if they are in the input",
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
    try{
        auto args = options.parse(argc, argv);
        if (args.count("help") || (args.count("inputfiles") == 0))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        auto& infiles = args["inputfiles"].as<std::vector<std::string>>();
        auto& outfile = args["output"].as<std::string>();
        auto& store_all = args["all"].as<bool>();
        auto& store_svt_hits = args["svt_hits"].as<bool>();
        auto& all_tracks = args["all_tracks"].as<bool>();
        auto& kf_tracks = args["kf_tracks"].as<bool>();
        auto& matched_tracks = args["matched_tracks"].as<bool>();
        auto& no_mc_particles = args["no_mc_particles"].as<bool>();
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
            if(dst2016->event->getNumberOfMCParticles() > 0 && !dst->use_mc_particles && !no_mc_particles){
                cout << "Warning: Input is MC Data, but write_mc_particles is not set. Turning on write_mc_particles!\n";
                dst->use_mc_particles = true;
            }
        }else if( infiles.size()>0 && infiles[0].find(".slcio") != string::npos ) {
            auto dstlcio = new LcioReader(infiles);
            dst = static_cast<MiniDst*>(dstlcio);
            if(!no_mc_particles) dst->use_mc_particles = true;  // will be set to false if no MCParticle collection in LcioReader.
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
        dst->use_ecal_cluster = store_all || args["ecal_clusters"].as<bool>();
        dst->use_ecal_hits = store_all || args["ecal_hits"].as<bool>();
        dst->use_svt_hits = store_svt_hits || store_all;
        dst->use_svt_raw_hits = args["raw_svt_hits"].as<bool>();
        dst->use_kf_tracks = kf_tracks || all_tracks || store_all;
        dst->use_gbl_tracks = all_tracks || store_all || args["gbl_tracks"].as<bool>();
        dst->use_gbl_kink_data = store_all || args["gbl_kinks"].as<bool>();
        dst->use_matched_tracks = matched_tracks || all_tracks || store_all;
        dst->use_kf_particles = !args["no_kf_particles"].as<bool>();
        dst->use_gbl_particles = !args["no_gbl_particles"].as<bool>();
        dst->SetOutputFileName(outfile);

        if(args["no_kf_particles"].as<bool>()){  // Erase and and all KF particle types in the output list.
            vector<int> copy_single(dst->particle_types_single); // make a copy
            dst->particle_types_single.clear();
            for(int p: copy_single){
                if(p >= dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_single.push_back(p);
            }
            vector<int> copy_double(dst->particle_types_double); // make a copy
            dst->particle_types_double.clear();
            for(int p: copy_double){
                if(p >= dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_single.push_back(p);
            }
        }

        if(args["no_gbl_particles"].as<bool>()){  // Erase and and all GBL particle types in the output list.
            vector<int> copy_single(dst->particle_types_single); // make a copy
            dst->particle_types_single.clear();
            for(int p: copy_single){
                if(p < dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_single.push_back(p);
            }
            vector<int> copy_double(dst->particle_types_double); // make a copy
            dst->particle_types_double.clear();
            for(int p: copy_double){
                if(p < dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_single.push_back(p);
            }
        }

#ifdef DEBUG
        cout << "Extra debug code compiled.\n";
        dst->Counter_Freq =10000;
#endif
        dst->Start();
        dst->Run(num_evt);
        dst->End();
    }catch(const cxxopts::OptionException& e){
        std::cout << "Error: " << e.what() << std::endl;
        std::cout << options.help() << std::endl;
        exit(0);
    }
}
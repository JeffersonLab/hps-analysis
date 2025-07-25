//
// Created by Maurik Holtrop on 7/11/20.
//

#include "cxxopts.hpp"
#include "TChain.h"
#include "MiniDst.h"
#include "Dst2016.h"

#include "LcioReader.h"
#include <locale.h>

using namespace std;

int main(int argc, char **argv){

   setlocale(LC_NUMERIC, "");

   // This is a nicer way to do options in C++. See cxxopts.hpp file.
   string help_string = string("Write a ROOT MiniDst for HPS data.\n") +
                        "Version: 1.1, using MiniDst.h version " + MiniDst::_version_() +
                        ", LcioReader version " + LcioReader::_version_() +
                        "\nCompiled with "+__VERSION__+"\n";
   cxxopts::Options options(argv[0], help_string);
   options
         .positional_help(" infile1 infile2 ...")
         .show_positional_help();

   options.add_options()
         ("d,debug", "Increase debug level")
         ("q,quiet", "Run quiet.")
         ("a,all", "Store all known values in file, except the raw stuff or uncorrected. Equivalent to -c -e -s -h",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("c,all_clusters", "Store Ecal and Hodo Clusters",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("E,ecal_clusters", "Store Ecal Clusters",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("U,ecal_clusters_uncor", "Store UNCORRECTED Ecal Clusters",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("e,ecal_hits", "Store Ecal Hits",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("H,hodo_clusters", "Store Hodo Clusters",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("h,hodo_hits", "Store Hodo Hits",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("s,svt_hits", "Store SVT 3D and/or strip Hits",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("r,raw_hits", "Store RAW Hits for SVT, ECal and Hodo",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("g,gbl_tracks", "Store GBL tracks",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("G,gbl_kinks", "Store GBL track kink data",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("k,kf_tracks", "Store KF tracks",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("t,matched_tracks", "Store Matched tracks",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("T,all_tracks", "Store all tracks (set gbl, kf and matched track to true)",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("X,extra_track","Store extra information from tracks, i.e. track states",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("no_kf_particles","Do NOT store the KF particles or vertexes.",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("gbl_particles","Also store the GBL particles or vertexes if available.",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("kf_has_postfix","The KF collections have the _KF postfix. (default is no postfix)",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("gbl_has_no_postfix","The GBL collections do not have the _GBL postfix. (default is _GBL)",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("m,magfield","Magnetic field strength to use in Tesla to compute px,py,pz",
               cxxopts::value<double>()->default_value("0."))
         ("M,no_mc_particles", "DO NOT Store MCParticles, even if they are in the input",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("x,use_mc_scoring", "Store extra MC output from scoring planes.",
          cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
         ("o,output", "Output file", cxxopts::value<std::string>()->default_value("mini_dst.root"))
         ("i,inputfiles",
          "List of input files which will be concatenated into a single output mini dst file. The -i "
          "is optional. ",
          cxxopts::value<std::vector<std::string>>())
         ("n,num_evt","Only process num_evt event", cxxopts::value<long>()->default_value("0"))
         ("help", "Print help")
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
      auto& all_tracks = args["all_tracks"].as<bool>();
      auto& kf_tracks = args["kf_tracks"].as<bool>();
      auto& matched_tracks = args["matched_tracks"].as<bool>();
      auto& num_evt = args["num_evt"].as<long>();
      int debug = 0;
      if(args.count("quiet") == 0 ){
         debug = args.count("debug") + 1;
      }else {
         debug = 0;
      }
      if(debug > 1) {
         cout << "Debug was set to: " << debug << endl;
      }

      if(debug) {
         cout << "Input: ";
         for (auto &v : infiles) {
            cout << " " << v;
         }
         cout << endl;
         cout << "Output: " << outfile << endl;
      }

      int debug_code = 0;
      if( debug <= 0) debug_code = 0;
      if (debug == 1) debug_code = MiniDst::kDebug_Error + MiniDst::kDebug_Warning;
      if (debug >= 2) debug_code += MiniDst::kDebug_Error + MiniDst::kDebug_Warning + MiniDst::kDebug_Info;
      if (debug >= 3) debug_code += MiniDst::kDebug_L1;
      if (debug >= 4) debug_code += MiniDst::kDebug_L2;
      if (debug >= 5) debug_code += MiniDst::kDebug_L3;
      if (debug >= 6) debug_code += MiniDst::kDebug_L4;

      MiniDst *dst{nullptr};
      bool is_dst_type = false;
      if( !infiles.empty() && infiles[0].find(".root") != string::npos ) {
         // The first file in the list has .root extension.
         is_dst_type = true;

         auto chain = new TChain("HPS_Event");
         for (auto &v : infiles) {
            chain->Add(v.c_str());
         }
         auto dst2016 = new Dst2016(chain);
         dst = static_cast<MiniDst *>(dst2016);
         if (dst2016->event->getNumberOfMCParticles() > 0 && !dst->use_mc_particles &&
             !args["no_mc_particles"].as<bool>()) {
            cout << "Warning: Input is MC Data, but write_mc_particles is not set. Turning on write_mc_particles!\n";
            dst->use_mc_particles = true;
         }
      }else if( infiles.size()>0 && infiles[0].find(".slcio") != string::npos ) {
         auto dstlcio = new LcioReader(infiles, debug_code);
         if(args["magfield"].as<double>() > 0.) dstlcio->magnetic_field = args["magfield"].as<double>();
         if(args["kf_has_postfix"].as<bool>()) dstlcio->kf_has_no_postscript = false;
         if(args["gbl_has_no_postfix"].as<bool>()) dstlcio->gbl_has_no_postscript = true;
         dst = static_cast<MiniDst*>(dstlcio);

//         // Slightly "expensive", but it is really nice to know ahead of time if we need MCParticle in the DST.
//         dstlcio->lcio_reader->open(infiles[0]);
//         auto lcio_event = dstlcio->lcio_reader->readNextEvent();
//         const vector<string> *col_names = lcio_event->getCollectionNames();
//
//         if( std::find(col_names->begin(), col_names->end(), "MCParticle") != col_names->end()
//             and !no_mc_particles) dst->use_mc_particles = true;
//         dstlcio->lcio_reader->close();
      }else{
         cout << "We need either an SLCIO file, or a DST2016 root file for input. \n Abort. \n";
         exit(1);
      }

      dst->SetDebugLevel(debug_code);
      dst->use_ecal_cluster = store_all || args["ecal_clusters"].as<bool>() || args["all_clusters"].as<bool>();
      dst->use_ecal_cluster_uncor = args["ecal_clusters_uncor"].as<bool>();
      dst->use_hodo_clusters = store_all || args["hodo_clusters"].as<bool>() || args["all_clusters"].as<bool>();
      dst->use_ecal_hits = store_all || args["ecal_hits"].as<bool>();
      dst->use_hodo_hits = store_all || args["hodo_hits"].as<bool>();
      dst->use_svt_hits =  store_all || args["svt_hits"].as<bool>();
      dst->use_svt_raw_hits = args["raw_hits"].as<bool>();
      dst->use_ecal_raw_hits = args["raw_hits"].as<bool>();
      dst->use_hodo_raw_hits = args["raw_hits"].as<bool>();
      dst->use_kf_tracks = kf_tracks || all_tracks || store_all;
      dst->use_gbl_tracks = all_tracks || args["gbl_tracks"].as<bool>();
      dst->use_gbl_kink_data = args["gbl_kinks"].as<bool>();
      dst->use_matched_tracks = matched_tracks || all_tracks;
      dst->use_extra_tracks = args["extra_track"].as<bool>();
      dst->use_kf_particles = !args["no_kf_particles"].as<bool>();
      dst->use_gbl_particles = args["gbl_particles"].as<bool>();
      dst->use_mc_particles = !args["no_mc_particles"].as<bool>();
      dst->use_mc_scoring = args["use_mc_scoring"].as<bool>();
      dst->SetOutputFileName(outfile);

      if(args["no_kf_particles"].as<bool>()){  // Erase all KF particle types in the output list.
         vector<int> copy_single(dst->particle_types_single); // make a copy
         dst->particle_types_single.clear();
         for(int p: copy_single){
            if(p >= dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_single.push_back(p);
         }
         vector<int> copy_double(dst->particle_types_double); // make a copy
         dst->particle_types_double.clear();
         for(int p: copy_double){
            if(p >= dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_double.push_back(p);
         }
      }

      if(!args["gbl_particles"].as<bool>()){  // Erase all GBL particle types in the output list.
         vector<int> copy_single(dst->particle_types_single); // make a copy
         dst->particle_types_single.clear();
         for(int p: copy_single){
            if(p < dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_single.push_back(p);
         }
         vector<int> copy_double(dst->particle_types_double); // make a copy
         dst->particle_types_double.clear();
         for(int p: copy_double){
            if(p < dst->FINAL_STATE_PARTICLE_GBL) dst->particle_types_double.push_back(p);
         }
      }

#ifdef DEBUG
      cout << "Extra debug code compiled.\n";
      dst->Counter_Freq =100;
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
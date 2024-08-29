//
//

#include "cxxopts.hpp"
#include <locale.h>
#include "GammaMixer.h"
#include "MollerMixer.h"

using namespace std;

int main(int argc, char **argv) {
    setlocale(LC_NUMERIC, "");

    cxxopts::Options options(argv[0], " - Write a ROOT MiniDst for HPS data.\n");
    options
            .positional_help(" infile1 infile2 ...")
            .show_positional_help();

    bool apple = false;

    options
            .add_options()
                    ("d,debug", "Increase debug level",cxxopts::value<int>()->default_value("1"))
                    ("C,counter", "Counter frequency",cxxopts::value<int>()->default_value("1000"))
                    ("q,quiet", "Run quiet.")
                    ("o,output", "Output file", cxxopts::value<std::string>()->default_value("mini_dst.root"))
                    ("i,inputfiles",
                     "List of input files which will be concatenated into a single output mini dst file. The -i "
                     "is optional. ",
                     cxxopts::value<std::vector<std::string>>())
                    ("n,num_evt","Only process num_evt event to output.", cxxopts::value<long>()->default_value("0"))
                    ("x,mix_mult","Mixing multiplier number", cxxopts::value<int>()->default_value("10"))
                    ("p,p_min","Minimum momentum for electron that is not FEE", cxxopts::value<double>()->default_value("2.0"))
                    ("m,esum_min","Energy difference tolerance for esum (1.9)", cxxopts::value<double>()->default_value("1.9"))
                    ("M,esum_max","Energy difference tolerance for esum (2.5)", cxxopts::value<double>()->default_value("2.5"))
                    ("g,gammas","Mix the photons of the events",cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    ("e,moller","Mix the electrons of the events",cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
                    // ("ep,trident","Mix the electron positron of the events",cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
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
    auto& num_evt = args["num_evt"].as<long>();
    auto& mix_mult= args["mix_mult"].as<int>();
    int debug = 0;
    if(args.count("quiet") == 0 ){
        debug = args["debug"].as<int>();
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

   if( !args["gammas"].as<bool>() && !args["moller"].as<bool>() ){
      cout << "Choose either gammas -g or mollers -ee to mix!\n";
   }

    if(args["gammas"].as<bool>()) {
       auto gamma_mixer = GammaMixer(infiles);
       gamma_mixer.SetDebugLevel(debug);
       gamma_mixer.SetOutputFileName(outfile);
       gamma_mixer.mix_multiplyer = mix_mult;
       gamma_mixer.Counter_Freq = args["counter"].as<int>();
       gamma_mixer.delta_esum_tolerance = 0.1;
       // args["esum_tolerance"].as<double>();
       gamma_mixer.Start();
       gamma_mixer.Run(num_evt);
       gamma_mixer.End();
    }
    if(args["moller"].as<bool>() ){
       auto moller_mixer = MollerMixer(infiles);
       moller_mixer.SetDebugLevel(debug);
       moller_mixer.SetOutputFileName(outfile);
       moller_mixer.mix_multiplyer = mix_mult;
       moller_mixer.Counter_Freq = args["counter"].as<int>();
       moller_mixer.one_electron_p_max = args["p_min"].as<double>();
       moller_mixer.two_electron_psum_min = args["esum_min"].as<double>();
       moller_mixer.two_electron_psum_max = args["esum_max"].as<double>();
       moller_mixer.Start();
       moller_mixer.Run(num_evt);
       moller_mixer.End();
    }
}
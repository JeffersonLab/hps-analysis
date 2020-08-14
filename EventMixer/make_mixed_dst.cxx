//
//

#include "cxxopts.hpp"
#include <locale.h>
#include "GammaMixer.h"

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
                    ("q,quiet", "Run quiet.")
                    ("o,output", "Output file", cxxopts::value<std::string>()->default_value("mini_dst.root"))
                    ("i,inputfiles",
                     "List of input files which will be concatenated into a single output mini dst file. The -i "
                     "is optional. ",
                     cxxopts::value<std::vector<std::string>>())
                    ("n,num_evt","Only process num_evt event to output.", cxxopts::value<long>()->default_value("0"))
                    ("x,mix_mult","Mixing multiplier number", cxxopts::value<int>()->default_value("10"))
                    ("t,esum_tolerance","Energy difference tolerance for esum (0.1)", cxxopts::value<double>()->default_value("0.1"))
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

    auto gamma_mixer = GammaMixer(infiles);
    gamma_mixer.SetDebugLevel(debug);
    gamma_mixer.SetOutputFileName(outfile);
    gamma_mixer.mix_multiplyer = mix_mult;
    gamma_mixer.Counter_Freq = 10000;
    gamma_mixer.delta_esum_tolerance=args["esum_tolerance"].as<double>();
    gamma_mixer.Start();
    gamma_mixer.Run(num_evt);
    gamma_mixer.End();
}
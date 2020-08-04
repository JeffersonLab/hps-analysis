//
// Created by Maurik Holtrop on 7/17/20.
//
#include "cxxopts.hpp"
#include "BumpHunt2016.h"

int main(int argc, char **argv) {

    cxxopts::Options options(argv[0], " - Write a ROOT MiniDst for HPS data.\n");
    options
            .positional_help(" infile1 infile2 ...")
            .show_positional_help();

    bool apple = false;

    options
            .add_options()
                    ("d,debug", "Increase debug level")
                    ("q,quiet", "Run quiet.")
                    ("m,multithread", "Turn on implicit multithreading.")
                    ("o,output", "Output file", cxxopts::value<std::string>()->default_value("mini_dst.root"))
                    ("i,inputfiles",
                     "List of input files which will be concatenated into a single output mini dst file. The -i "
                     "is optional. ",
                     cxxopts::value<std::vector<std::string>>())
                    ("n,num_evt", "Only process num_evt event", cxxopts::value<long>()->default_value("0"))
                    ("h,help", "Print help");


    options.parse_positional({"inputfiles"});
    auto args = options.parse(argc, argv);
    if (args.count("help") || (args.count("inputfiles") == 0)) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    auto &infiles = args["inputfiles"].as<std::vector<std::string>>();
    auto &outfile = args["output"].as<std::string>();
    int debug = args.count("debug");
    int quiet = args.count("quiet");
    if (debug) {
        std::cout << "Debug was set to: " << debug << std::endl;
        std::cout << "Input: ";
        for (auto &v : infiles) {
            std::cout << " " << v;
        }
        std::cout << std::endl;
        std::cout << "Output: " << outfile << std::endl;
    }
    if(args.count("multithread")){
        std::cout << "Turn on Multi-Threading \n";
        ROOT::IsImplicitMTEnabled();
    }
    auto bh = new BumpHunt2016(infiles[0]);
    bh->Setup();
    bh->Process();
}
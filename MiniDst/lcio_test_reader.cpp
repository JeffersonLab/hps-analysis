//
// Created by Maurik Holtrop on 8/6/20.
//

#include <vector>
#include <iostream>

#include "cxxopts.hpp"
#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>

using namespace std;

int main(int argc, char **argv) {

    cxxopts::Options options(argv[0], " - Write a ROOT MiniDst for HPS data.\n");
    options
            .positional_help(" infile1 infile2 ...")
            .show_positional_help();

    bool apple = false;

    options
            .add_options()
                    ("d,debug", "Increase debug level")
                    ("i,inputfiles",
                     "List of input files which will be concatenated into a single output mini dst file. The -i "
                     "is optional. ",
                     cxxopts::value<std::vector<std::string>>());

    options.parse_positional({"inputfiles"});
    auto args = options.parse(argc, argv);

    vector<string> infiles;
    if(args["inputfiles"].count()) {
        infiles = args["inputfiles"].as<std::vector<std::string>>();
    }else{
        infiles = {"/data/HPS/data/physrun2019/slcio/hps_010735_00000.slcio"};
    }
    int debug = args.count("debug");

    IO::LCReader* lcio_reader{IOIMPL::LCFactory::getInstance()->createLCReader()};
    EVENT::LCEvent* lcio_event{nullptr};
    int run_number;
    int event_number;

    for(auto file: infiles) {

        cout << "-------------------------------------------------------------------------------\n";
        cout << "LCIO::Open - " << file << "\n\n";
        lcio_reader->open(file);

        bool first_event_read{false};

        while ((lcio_event = lcio_reader->readNextEvent())) {
            run_number = lcio_event->getRunNumber();
            event_number = lcio_event->getEventNumber();

            if(!first_event_read) {
                const vector<string> *col_names = lcio_event->getCollectionNames();
                cout << "Run Number: " << run_number << " first event: " << event_number << "\n\n";
                for( auto col: (*col_names) ){
                    cout << " " << col;
                }
                cout << "\n";
                first_event_read = true;
            }


        }
        cout << "\nRun " << run_number << " last event: " << event_number << "\n\n";
        lcio_reader->close();
    }
}
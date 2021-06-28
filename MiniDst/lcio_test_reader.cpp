//
// Created by Maurik Holtrop on 8/6/20.
//

#include <vector>
#include <iostream>

#include "cxxopts.hpp"
#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCGenericObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/ParticleID.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerRawData.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <UTIL/BitField64.h>
#include <locale.h>

using namespace std;

void print_collection_id(EVENT::LCEvent* lcio_event, string name){
    // Print collection ids for tracker like collection name.
    EVENT::LCCollection* collection = lcio_event->getCollection(name.c_str());
    unsigned long long raw_value{0};
    for(int i=0; i < collection->getNumberOfElements(); ++i){
        auto track_hit = dynamic_cast<IMPL::TrackerHitImpl*>(collection->getElementAt(i));
        unsigned long long value = ((unsigned long long)(track_hit->getCellID0()) & 0xffffffff) |
                                   ( (unsigned long long)(track_hit->getCellID1()) << 32);

        EVENT::LCObjectVec raw_hits = track_hit->getRawHits();
        if(!raw_hits.empty()) {
            auto lcio_raw_hit = dynamic_cast<EVENT::TrackerRawData *>(raw_hits.at(0));
            raw_value = ((unsigned long long) (lcio_raw_hit->getCellID0()) & 0xffffffff) |
                                       ((unsigned long long) (lcio_raw_hit->getCellID1()) << 32);
        }else{
            raw_value = 0;
        }

        printf("ID value from %35s: %'20llu =? %'20llu\n", name.c_str(), value, raw_value);
    }
}

int main(int argc, char **argv) {

    setlocale(LC_NUMERIC, "");
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

        int count{0};
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

            if( ++ count < 10) {
                print_collection_id(lcio_event, "StripClusterer_SiTrackerHitStrip1D");
                print_collection_id(lcio_event, "HelicalTrackHits");
                print_collection_id(lcio_event, "RotatedHelicalTrackHits");
            }
        }
        cout << "\nRun " << run_number << " last event: " << event_number << "\n\n";
        lcio_reader->close();
    }
}
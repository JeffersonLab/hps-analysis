//
// Created by Maurik Holtrop on 7/24/20.
//

#ifndef MINIDST_LCIOREADER_H
#define MINIDST_LCIOREADER_H

#include <string>
#include <vector>
#include <iostream>

#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCGenericObject.h>
#include <EVENT/LCCollection.h>

#include "MiniDst.h"

using namespace std;

class LcioReader : public MiniDst {

public:
    LcioReader(const string input_file=""){
        if(!input_file.empty()) input_file_list.push_back(input_file);
    };
    LcioReader(const vector<string> infile_list): input_file_list(infile_list){};
    ~LcioReader(){};


    virtual void Start() override;
    virtual long Run(int nevt=0) override;
    virtual void End() override;

public:
    IO::LCReader* lcio_reader{IOIMPL::LCFactory::getInstance()->createLCReader()};
    EVENT::LCEvent* lcio_event{nullptr};
    vector<string> input_file_list{};

    unsigned long evt_count{0}; // Event sequence number.
    unsigned long Counter_Freq{10}; // How often to print a status line.
    const vector<string> *col_names; // Store the collection names from the first event.
    bool data_type_is_known{false};  // The LCIO data is different between 2015/2016 and 2019. This is true when that is known.
    bool is_2016_data{false};  // True for 2015 and 2016 data: i.e. there is Trigger info in the TriggerBank
    bool is_2019_data{false};  // True for 2019 data: i.e. there is a TSBank and a VTPBank.
    bool is_MC_data{false};    // True is there is an MCParticles bank.

};


#endif //MINIDST_LCIOREADER_H

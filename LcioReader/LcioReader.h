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

#include "MiniDst.h"

using namespace std;

class LcioReader : public MiniDst {

public:
    LcioReader(const string_view input_file=""){
        if(!input_file.empty()) input_file_list.push_back(input_file);
    };
    LcioReader(const vector<string_view> infile_list): input_file_list(infile_list){};
    ~LcioReader(){};


    virtual void Start() override;
    virtual long Run(int nevt=0) override;
    virtual void End() override;

public:
    IO::LCReader* lcio_reader{IOIMPL::LCFactory::getInstance()->createLCReader()};
    EVENT::LCEvent* lcio_event{nullptr};
    vector<string_view> input_file_list{};

};


#endif //MINIDST_LCIOREADER_H

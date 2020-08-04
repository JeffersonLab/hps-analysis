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

#include "MiniDst.h"

using namespace std;

class LcioReader : public MiniDst {

    const vector<string> Type_to_Collection{
        "FinalStateParticles",  // FINAL_STATE_PARTICLE: 0
        "UnconstrainedV0Candidates", // UC_V0_CANDIDATE: 1,
        "BeamspotConstrainedV0Candidates", // BSC_V0_CANDIDATE: 2,
        "TargetConstrainedV0Candidates", // TC_V0_CANDIDATE: 3,
        "UnconstrainedMollerCandidates", // UC_MOLLER_CANDIDATE: 4,
        "BeamspotConstrainedMollerCandidates", // BSC_MOLLER_CANDIDATE: 5,
        "TargetConstrainedMollerCandidates", // TC_MOLLER_CANDIDATE: 6,
        "OtherElectrons", // OTHER_ELECTRONS: 7,
        "UnconstrainedVcCandidates", // UC_VC_CANDIDATE: 8
    };
    const vector<string> Type_to_VertexCollection{
            "",  // FINAL_STATE_PARTICLE: 0
            "UnconstrainedV0Vertexes", // UC_V0_CANDIDATE: 1,
            "BeamspotConstrainedV0Vertices", // BSC_V0_CANDIDATE: 2,
            "TargetConstrainedV0Vertices", // TC_V0_CANDIDATE: 3,
            "UnconstrainedMollerVertices", // UC_MOLLER_CANDIDATE: 4,
            "BeamspotConstrainedMollerVertices", // BSC_MOLLER_CANDIDATE: 5,
            "TargetConstrainedMollerVertices", // TC_MOLLER_CANDIDATE: 6,
            "", // OTHER_ELECTRONS: 7,
            "UnconstrainedVcVertices", // UC_VC_CANDIDATE: 8
    };

public:
    LcioReader(const string input_file=""){
        if(!input_file.empty()) input_file_list.push_back(input_file);
    };
    LcioReader(const vector<string> infile_list): input_file_list(infile_list){};
    ~LcioReader(){};


    virtual void Start() override;
    virtual long Run(int nevt=0) override;
    virtual void End() override;

    virtual void Fill_Part_From_LCIO(Basic_Particle_t *bp,EVENT::ReconstructedParticle *lcio_part);
    virtual void Fill_SubPart_From_LCIO(Sub_Particle_t *sub,EVENT::ReconstructedParticle *daughter,
                                        EVENT::LCEvent *lcio_event);

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

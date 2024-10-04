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
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <UTIL/BitField64.h>

#include "MiniDst.h"

#define LCIOReader__Version__ "1.2.1"
using namespace std;

// Collection in May 2021 versions of the 2019 data with both KF and GBL tracking.
//
//  BeamspotConstrainedV0Candidates BeamspotConstrainedV0Candidates_KF BeamspotConstrainedV0Vertices
//  BeamspotConstrainedV0Vertices_KF EcalCalHits EcalClusters EcalClustersCorr EcalReadoutHits EcalUncalHits
//  FADCGenericHits FinalStateParticles FinalStateParticles_KF GBLKinkData GBLKinkDataRelations GBLTracks
//  HelicalTrackHitRelations HelicalTrackHits HodoCalHits HodoGenericClusters HodoReadoutHits KFGBLStripClusterData
//  KFGBLStripClusterDataRelations KFTrackData KFTrackDataRelations KalmanFullTracks MatchedToGBLTrackRelations
//  MatchedTracks OtherElectrons OtherElectrons_KF RFHits RotatedHelicalTrackHitRelations RotatedHelicalTrackHits
//  SVTFittedRawTrackerHits SVTRawTrackerHits SVTShapeFitParameters StripClusterer_SiTrackerHitStrip1D TSBank
//  TargetConstrainedV0Candidates TargetConstrainedV0Candidates_KF TargetConstrainedV0Vertices
//  TargetConstrainedV0Vertices_KF TrackData TrackDataRelations TriggerBank UnconstrainedV0Candidates
//  UnconstrainedV0Candidates_KF UnconstrainedV0Vertices UnconstrainedV0Vertices_KF UnconstrainedVcCandidates
//  UnconstrainedVcCandidates_KF UnconstrainedVcVertices UnconstrainedVcVertices_KF VTPBank


class LcioReader : public MiniDst {
//
// Type_to_collection and Type_to_VertexCollection are used to map the type of collection to the collection name.
// Unfortunately, in the *data* these collection types are controlled in the lcsim steering file, and are thus
// not consistent between various data sets.
// In older data sets the _GBL marker will not be present. In newer data sets the _KF marker will be omitted.
// Which one is which can be switched by the kf_has_not_postscript flag.
//
    vector<string> Type_to_Collection{
        "FinalStateParticles_KF",  // FINAL_STATE_PARTICLE: 0
        "UnconstrainedV0Candidates_KF", // UC_V0_CANDIDATE: 1,
        "BeamspotConstrainedV0Candidates_KF", // BSC_V0_CANDIDATE: 2,
        "TargetConstrainedV0Candidates_KF", // TC_V0_CANDIDATE: 3,
        "UnconstrainedMollerCandidates_KF", // UC_MOLLER_CANDIDATE: 4,
        "BeamspotConstrainedMollerCandidates_KF", // BSC_MOLLER_CANDIDATE: 5,
        "TargetConstrainedMollerCandidates_KF", // TC_MOLLER_CANDIDATE: 6,
        "OtherElectrons_KF", // OTHER_ELECTRONS: 7,
        "UnconstrainedVcCandidates_KF", // UC_VC_VERTICES_GBL: 8
        "FinalStateParticles_GBL",  // FINAL_STATE_PARTICLE_GBL: 0+9
        "UnconstrainedV0Candidates_GBL", // UC_V0_VERTICES_GBL: 1+9
        "BeamspotConstrainedV0Candidates_GBL", // BSC_V0_VERTICES_GBL: 2+9
        "TargetConstrainedV0Candidates_GBL", // TC_V0_VERTICES_GBL: 3+9,
        "UnconstrainedMollerCandidates_GBL", // UC_MOLLER_VERTICES_GBL: 4+9,
        "BeamspotConstrainedMollerCandidates_GBL", // BSC_MOLLER_VERTICES_GBL: 5+9,
        "TargetConstrainedMollerCandidates_GBL", // TC_MOLLER_VERTICES_GBL: 6+9,
        "OtherElectrons_GBL", // OTHER_ELECTRONS_GBL: 7+9,
        "UnconstrainedVcCandidates_GBL", // UC_VC_VERTICES_GBL: 8+9
    };
    vector<string> Type_to_VertexCollection{
            "",  // FINAL_STATE_PARTICLE: 0
            "UnconstrainedV0Vertices_KF", // UC_V0_CANDIDATE: 1,
            "BeamspotConstrainedV0Vertices_KF", // BSC_V0_CANDIDATE: 2,
            "TargetConstrainedV0Vertices_KF", // TC_V0_CANDIDATE: 3,
            "UnconstrainedMollerVertices_KF", // UC_MOLLER_CANDIDATE: 4,
            "BeamspotConstrainedMollerVertices_KF", // BSC_MOLLER_CANDIDATE: 5,
            "TargetConstrainedMollerVertices_KF", // TC_MOLLER_CANDIDATE: 6,
            "", // OTHER_ELECTRONS: 7,
            "UnconstrainedVcVertices_KF", // UC_VC_CANDIDATE: 8
            "",  // FINAL_STATE_PARTICLE: 0
            "UnconstrainedV0Vertices_GBL", // UC_V0_VERTICES_GBL: 1,
            "BeamspotConstrainedV0Vertices_GBL", // BSC_V0_VERTICES_GBL: 2,
            "TargetConstrainedV0Vertices_GBL", // TC_V0_VERTICES_GBL: 3,
            "UnconstrainedMollerVertices_GBL", // UC_MOLLER_VERTICES_GBL: 4,
            "BeamspotConstrainedMollerVertices_GBL", // BSC_MOLLER_VERTICES_GBL: 5,
            "TargetConstrainedMollerVertices_GBL", // TC_MOLLER_VERTICES_GBL: 6,
            "", // OTHER_ELECTRONS_GBL: 7,
            "UnconstrainedVcVertices_GBL", // UC_VC_VERTICES_GBL: 8
    };

public:
    explicit LcioReader(const string &input_file="", int debug_level=0x07);
    explicit LcioReader(const vector<string> &infile_list, int debug_level=0x07);
    ~LcioReader() override = default;

    static string _version_(){return(LCIOReader__Version__);};
    void Clear() override;
    void Start() override;
    void SetupLcioDataType();
    bool Process(Long64_t entry) override;
    long Run(int max_event) override;
    void End() override;

    void Fill_Basic_Particle_From_LCIO(Basic_Particle_t *bp, EVENT::ReconstructedParticle *lcio_part,
                                               bool fill_momentum=true);
    void Fill_Single_Particle_From_LCIO(Single_Particle_t *bp,EVENT::ReconstructedParticle *lcio_part, int type);
    void Fill_Vertex_From_LCIO(Vertex_Particle_t *bp,EVENT::Vertex *lcio_vert, int type);
    void Fill_SubPart_From_LCIO(Sub_Particle_t *sub,EVENT::ReconstructedParticle *daughter, int type);

    bool has_collection(const char *name) const{
        return(std::find(col_names->begin(), col_names->end(), name) != col_names->end());}

public:
    IO::LCReader* lcio_reader{IOIMPL::LCFactory::getInstance()->createLCReader()};
    EVENT::LCEvent* lcio_event{nullptr};
    vector<string> input_files{};

    unsigned long evt_count{0}; // Event sequence number.
    const vector<string> *col_names{nullptr}; // Store the LCIO collection names from the first event.
    vector<string> svt_hit_collections;
    bool data_type_is_known{false};  // The LCIO data is different between 2015/2016 and 2019. This is true when that is known.
    bool is_2016_data{false};  // True for 2015 and 2016 data: i.e. there is Trigger info in the TriggerBank
    bool is_2019_data{false};  // True for 2019 data: i.e. there is a TSBank and a VTPBank.
    bool is_MC_data{false};    // True is there is an MCParticles bank.
    bool has_rf_hits{true};
    bool kf_has_not_postscript{true}; // Set to true if the Vertex etc collections have no postscript of _KF.
    bool gbl_has_no_postscript{false}; // GBL has not postscript of _GBL.

    double magnetic_field{0.};

    /// Maps to help navigate the event.
    /// Set them up here so that they are available in different sections.
    /// Note: Each map must be cleared at the start of the event to avoid potential cross event contamination.
    /// Clearing is done in Clear().
    // map<IMPL::CalorimeterHitImpl*, int> hodo_hit_to_index_map; // Map to link hodoscope hit to index.
    map<IMPL::CalorimeterHitImpl*, int> ecal_hit_to_index_map; // Map to link calorimeter hit to index.
    map<int, int> ecal_id0_to_hit_index;                      // Map the hit cellid0 to the index
    map<IMPL::ClusterImpl*, int> ecal_cluster_to_index_map; // Map to link calorimeter hit to index.
    map<IMPL::ClusterImpl*, int> ecal_cluster_uncor_to_index_map; // Map to link uncor calorimeter hit to index.
    map<IMPL::TrackerHitImpl*, int> svt_hit_to_index_map; // Map to link svt hit to index.
    map<EVENT::TrackerRawData *, pair<int,int>> svt_raw_hit_to_index_map; // Map of raw hit to indexes of stored raw hit and fits.
    map<EVENT::Track*, int> kf_track_to_index_map;           // Map to link GBL only track to index.
    map<EVENT::Track*, int> gbl_track_to_index_map;           // Map to link GBL only track to index.
    map<EVENT::Track*, int> matched_track_to_index_map;       // Map to link Matched only track to index.
    map<EVENT::Track*, int> any_track_to_index_map;           // Map to link any track to index.
    map<EVENT::ReconstructedParticle*, int> any_particle_to_index_map; // Map to link any particle to the particle index.

private:
   /// LCIO Decoders. We need these instantiated only once.
   ///
   UTIL::BitField64 hodo_hit_field_decoder{"system:6,barrel:3,layer:4,ix:4,iy:-3,hole:-3"};
   UTIL::BitField64 ecal_hit_field_decoder{"system:6,layer:2,ix:-8,iy:-6"};
   UTIL::BitField64 raw_svt_hit_decoder{"system:6,barrel:3,layer:4,module:12,sensor:1,side:32:-2,strip:12"};

};


#endif //MINIDST_LCIOREADER_H

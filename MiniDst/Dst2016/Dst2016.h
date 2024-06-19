///
/// An interface class to read the 2016 DST format from hps-dst and convert it to the MiniDST format.
///

#ifndef HPS_ANALYSIS_DST2016_H
#define HPS_ANALYSIS_DST2016_H
#include <cstdio>
#include <iostream>
#include <variant>
#include "TObject.h"
#include "TTree.h"
#include "TVector3.h"

#include "BaseAna.h"
#include "MiniDst.h"

///
/// We derive from BaseAna, with implements a ROOT TSelector. For this code, that is a bit of
/// overkill, but it was the quickest way for me to implement a translator from the old 2016 DST to
/// the MiniDST format. The whole BaseAna and associated classes are far more complicated than needed,
/// but that is what happens when a format gets abandoned.
///
class Dst2016 :  public MiniDst, public BaseAna {

public:
    explicit Dst2016(TTree *tree=nullptr, string out_file_name="minidst2016.root");
    ~Dst2016(){};

    void Clear() override;
    virtual void Start() override;
    virtual long Run(int nevt) override;
    virtual void End() override;
    virtual void SlaveBegin(TTree *tree= nullptr) override;
    virtual Bool_t  Process(Long64_t entry) override;
    virtual void    Terminate() override;
    virtual void SetDebugLevel(const int level) override {
        md_Debug = level;
        fDebug=level;};

   vector<int> particle_types_single_dst2016{HpsParticle::FINAL_STATE_PARTICLE, HpsParticle::OTHER_ELECTRONS};

   // Yes vertex.
   vector<int> particle_types_double_dst2016{HpsParticle::UC_V0_CANDIDATE, HpsParticle::BSC_V0_CANDIDATE,
                                             HpsParticle::TC_V0_CANDIDATE, HpsParticle::UC_VC_CANDIDATE,
                                             HpsParticle::UC_MOLLER_CANDIDATE, HpsParticle::BSC_MOLLER_CANDIDATE,
                                             HpsParticle::TC_MOLLER_CANDIDATE};

   const map<int, int> particle_type_translate_dst2016_to_minidst{
         {HpsParticle::FINAL_STATE_PARTICLE, MiniDst::FINAL_STATE_PARTICLE_GBL}, // Note, no KF in dst2016
         {HpsParticle::UC_V0_CANDIDATE, MiniDst::UC_V0_VERTICES_GBL},
         {HpsParticle::BSC_V0_CANDIDATE, MiniDst::BSC_V0_VERTICES_GBL},
         {HpsParticle::TC_V0_CANDIDATE, MiniDst::TC_V0_VERTICES_GBL},
         {HpsParticle::UC_MOLLER_CANDIDATE, MiniDst::UC_MOLLER_VERTICES_GBL},
         {HpsParticle::BSC_MOLLER_CANDIDATE, MiniDst::BSC_MOLLER_VERTICES_GBL},
         {HpsParticle::TC_MOLLER_CANDIDATE, MiniDst::TC_MOLLER_VERTICES_GBL},
         {HpsParticle::OTHER_ELECTRONS, MiniDst::OTHER_ELECTRONS_GBL},
         {HpsParticle::UC_VC_CANDIDATE,MiniDst::UC_VC_VERTICES_GBL}
   };


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
ClassDef(Dst2016, 1);
#pragma clang diagnostic pop
};


#endif //HPS_ANALYSIS_DST2016_H

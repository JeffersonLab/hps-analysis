///
/// An interface class to read the 2016 DST format from hps-dst and convert it to the MiniDST format.
///

#ifndef HPS_ANALYSIS_DST2016_H
#define HPS_ANALYSIS_DST2016_H
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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
ClassDef(Dst2016, 1);
#pragma clang diagnostic pop
};


#endif //HPS_ANALYSIS_DST2016_H

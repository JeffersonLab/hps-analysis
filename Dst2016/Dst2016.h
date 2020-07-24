//
// Created by Maurik Holtrop on 7/6/20.
//

#ifndef HPS_ANALYSIS_DST2016_H
#define HPS_ANALYSIS_DST2016_H
#include <iostream>
#include <variant>
#include "TObject.h"
#include "TTree.h"
#include "TVector3.h"

#include "BaseAna.h"
#include "MiniDst.h"

class Dst2016 :  public MiniDst, public BaseAna {

public:
    explicit Dst2016(TTree *tree=nullptr, string out_file_name="minidst2016.root");
    ~Dst2016(){};

    virtual void Start() override;
    virtual long Run(int nevt) override;
    virtual void End() override;
    virtual void SlaveBegin(TTree *tree= nullptr) override;
    virtual Bool_t  Process(Long64_t entry) override;
    virtual void    Terminate() override;

    void            SetOutputFileName(const string& outfile){md_output_file_name=outfile;};

    void clear();  // Clear all the vectors.
    //std::map<std::string, vector<double>* > &Get_brmap(){return branch_map;};
    // #define FULL_CLEAR(xxx)  { for( auto p: xxx){delete p;}; xxx.clear(); };

    template<typename T> inline void FULL_CLEAR(T& xxx){
        for(auto p: xxx){ delete p; };
        xxx.clear();
    }

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
ClassDef(Dst2016, 1);
#pragma clang diagnostic pop
};


#endif //HPS_ANALYSIS_DST2016_H

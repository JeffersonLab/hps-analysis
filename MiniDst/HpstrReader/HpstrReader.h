//
// Created by Maurik Holtrop on 8/30/24.
//

#ifndef HPSTRREADER_H
#define HPSTRREADER_H

#include "MiniDst.h"

#include "Event.h"
#include "HodoHit.h"
#include "HodoCluster.h"
#include "CalCluster.h"
#include "CalHit.h"


#define HpstrReader__Version__ "0.0.1"

using namespace std;

class HpstrReader: public MiniDst{

public:
   explicit HpstrReader(const string &input_file="", int debug_level=0x07);
   explicit HpstrReader(const vector<string> &infile_list, int debug_level=0x07);
   explicit HpstrReader(TChain *chain, int debug_level=0x07);
   ~HpstrReader() override = default;

   static string _version_(){return(HpstrReader__Version__);};

   void Clear() override;
   void Start() override;
   void Init();
   bool Process(Long64_t entry) override;

public:
   vector<string> input_files{};

   TChain *fChain;   //!pointer to the analyzed TTree or TChain

   EventHeader *event_header{nullptr};
   TBranch *b_event_header{nullptr};
   HodoHit *hodo_hit{nullptr};
   TBranch *b_hodo_hit{nullptr};
   HodoCluster *hodo_cluster{nullptr};
   TBranch *b_hodo_cluster{nullptr};
};


#endif //HPSTRREADER_H

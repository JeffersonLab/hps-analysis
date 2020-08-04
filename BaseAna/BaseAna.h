/*
 *  BaseAna.h
 */
/*
 * This class helps with the reading of HPS DSTs by providing a convenient way to access the
 * data from either a TTree, a TChain or when using PROOF (Parallel ROOT).
 *
 * For basic usage, consult the GitHub pages.
 *
 */
// TODO:
//       * Create copy constuctor?
//

#ifndef BaseAna_h
#define BaseAna_h

#include <iostream>

using namespace std;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#include <TROOT.h>
#include <TProofServ.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include "HpsEvent.h"
#pragma clang diagnostic pop

//
// I *think* we cannot derive from HpsEvent, and still make all the streaming work, since
// there would be an issue setting the branch address, given that we we want to add some
// contant in the BaseAna class.
// Having HpsEvent as a member, and creating pass-through calls, is safer.
//
class BaseAna : public TSelector {
public:
  TTree  *fChain{nullptr};         //! Chain holding events.
  Int_t   fChain_type{0};    //!  Controls the behavior of GetEntry: 1=TTree, 2= TChain, 3=Proof
  Int_t   fCurrent{0};       //! Current Tree number is a TChain
  Int_t   fDebug{0};
  int     fCounter_Freq{1<<30};
  bool    is_process{false};

  HpsEvent *event{nullptr};   // Stores the HpsEvent
  TBranch  *b_event{nullptr}; // Branch to HpsEvent

  long    current_event{0}; // Current event in file
  long    evt_count{0};     // Current sequence number for event.

  string  output_file_name;
  TFile   *output_file{nullptr};

  double b_field{0.};
  
  enum Debug_codes {
    kDebug_Quiet   = 0x00,
    kDebug_Error   = 0x01,
    kDebug_Warning = 0x02,
    kDebug_Info    = 0x04,
    kDebug_L1      = 0x10,  // Debug levels
    kDebug_L2      = 0x20,
    kDebug_L3      = 0x30,
    kDebug_L4      = 0x40
  };
  
  enum Chain_types {
    kIs_TTree  = 1,
    kIs_TChain = 2,
    kIs_TChain_Process = 3,
    kIs_TProof = 4
  };
  
  static const int ecal_nx{23};
  static const int ecal_ny{5};
  
public:
  BaseAna(TTree *tree=NULL,string out_file_name="Baseana.root");
  virtual ~BaseAna();
  virtual Int_t Version(void) const { return 1; }

  virtual void    Begin(TTree *tree=nullptr);
  virtual void    SlaveBegin(TTree *tree=nullptr /* tree */);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0);
  virtual Int_t   Next(void);
  virtual void    SetOption(const char *option) { fOption = option; };
  virtual void    SetObject(TObject *obj) { fObject = obj; };
  virtual void    SetInputList(TList *input) { fInput = input; };
  virtual TList  *GetOutputList() const{ return fOutput; };
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  virtual void    Start() { Begin(); SlaveBegin();};
  virtual void    End() { SlaveTerminate(); Terminate();};
// Debug Manipulators:
  virtual void    SetDebugLevel(const int level){ fDebug = level;};
  virtual void    SetDebugLevelBit(const int level){ fDebug = fDebug | level;};
  virtual void    ClearDebugLevelBit(const int level){ fDebug = fDebug & ~level;};
  virtual void    SetDebugLevelVerbose(void){ fDebug = kDebug_Error +kDebug_Warning +kDebug_Info;};
  virtual int     GetDebugLevel(void){ return fDebug; };

  
// Useful extra methods not needed for PROOF
  
  virtual long             Run(int nevent=0);
  virtual void            PrintParticle(const HpsParticle *part,int i=0) const;
  virtual void            Print(Option_t *opt="") const;
//
// HpsEvent pass through getters, with some extras.
//
  virtual HpsParticle    *GetFSParticle(int n);
  virtual double         GetAbsMomentum(HpsParticle *part);
  virtual double         GetFSAbsMomentum(int n);
  
  virtual EcalCluster*   GetEcalCluster(int n) const {return event->getEcalCluster(n);};
  virtual EcalHit*       GetEcalHit(int n) const {return event->getEcalHit(n);};
  virtual MCParticle*    GetMCParticle(int n) const {return event->getMCParticle(n);};
  virtual int            GetEventNumber() const { return event->getEventNumber(); };
  virtual long           GetEventTime() const { return event->getEventTime(); };
  virtual int            GetNumberOfEcalHits() const { return event->getNumberOfEcalHits(); };
  virtual int            GetNumberOfEcalClusters() const {return event->getNumberOfEcalClusters();};
  virtual int            GetNumberOfMCParticles() const { return event->getNumberOfMCParticles(); };
  virtual int            GetNumberOfParticles(HpsParticle::ParticleType type) const {return event->getNumberOfParticles(type);};
  virtual int            GetNumberOfSvtHits() const { return event->getNumberOfSvtHits(); };
  virtual int            GetNumberOfTracks()  const { return event->getNumberOfTracks();  };
  virtual int            GetNumberOfGblTracks() const { return event->getNumberOfGblTracks(); };
  virtual HpsParticle *  GetParticle(HpsParticle::ParticleType type, int n) const {return event->getParticle(type, n);};
  virtual double         GetRfTime(const int channel) const { return event->getRfTime(channel); };
  virtual int            GetRunNumber() const  { return event->getRunNumber(); };
  virtual SvtHit*        GetSvtHit(int n) const {return event->getSvtHit(n);};
  virtual SvtTrack*      GetTrack(int n)  const {return event->getTrack(n);};
  virtual GblTrack*      GetGblTrack(int n) const {return event->getGblTrack(n);};
  virtual bool           IsPair0Trigger() const { return event->isPair0Trigger(); };
  virtual bool           IsPair1Trigger() const { return event->isPair1Trigger(); };
  virtual bool           IsPulserTrigger() const { return event->isPulserTrigger(); };
  virtual bool           IsSingle0Trigger() const { return event->isSingle0Trigger(); };
  virtual bool           IsSingle1Trigger() const { return event->isSingle1Trigger(); };
  virtual bool           IsSvtBiasOn() const { return event->isSvtBiasOn(); };
  virtual bool           IsSvtClosed() const { return event->isSvtClosed(); };
  virtual bool           IsSvtLatencyGood() const { return event->isSvtLatencyGood(); };
  virtual bool           HasSvtBurstModeNoise() const { return event->hasSvtBurstModeNoise(); };
  virtual bool           HasSvtEventHeaderErrors() const { return event->hasSvtEventHeaderErrors(); };
  
  virtual HpsEvent *     GetEvent() const {return event;};
  
  void            SetOutputFileName(const string& outfile){output_file_name=outfile;};
  string          GetOutputFileName(void){return(output_file_name);};
  void            WriteList(TList *ll);

  void            SetBField(double field){b_field=field;};
  double          GetBField(void){return(b_field);};
  
  void            DrawEcal(int n=2);
  TH2F           *EcalHitMap(void);
  TH2F           *ClusterMap(int n_cl=0);
  void            FancyPlot(TH2F *histo,int opt);
  inline int      TrXhit(const int x) const{ return (x>=0?x:x+1); };
  
  ClassDef(BaseAna,1);
};

#endif

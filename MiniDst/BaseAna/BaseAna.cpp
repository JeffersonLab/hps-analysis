/*!
 * \class  BaseAna
 * \brief  The base analysis class for hps-dst based analysis projects.
 *
 * The following methods are defined in this file:
 *  -  Begin():          called every time a loop on the tree starts,
 *                       a convenient place to create your histograms.
 *  -  SlaveBegin():     called after Begin(), when on PROOF called only on the
 *                       slave servers.
 *  -  Process():        called for each event, in this function you decide what
 *                       to read and fill your histograms.
 *  -  SlaveTerminate(): called at the end of the loop on the tree, when on PROOF
 *                       called only on the slave servers.
 *  -  Terminate():      called at the end of the loop on the tree,
 *                       a convenient place to draw/fit your histograms.
 *
 * ##Use:
 *
 *  The intent is to be able to use this class with PROOF, however, for practicality, it can also be used directly with a TChain (or TTree),
 *  or given a file and a subsequent Run().
 *
 *  Use with "Run()":
 *  This example assumes a class derived from BaseAna called MyAna, which has a Process() that actually does something.
 *
 *  auto ch = new TChain("HPS_Event");
 *  ch->Add("*.root");
 *  auto ma = new MyAna(ch);
 *  ch->Start()
 *  ch->Run()
 *  ch->End()
 *
 *  Use with TChain and TProof:
 *  auto ch = new TChain("HPS_Event");
 *  ch->Add("*.root");
 *  TProof::Open("workers=8"); // Skip this and the following for single threaded.
 *  gProof->Exec("R__LOAD_LIBRARY(/data/HPS/lib/libMyAna.dylib)")
 *  ch->SetProof();            // Last line to Skip
 *  auto ma = new MyAna();
 *  ch->Process(ma)
*/

#include <iostream>
#include "BaseAna.h"
#include "HpsParticle.h"

ClassImp(BaseAna);

BaseAna::BaseAna(TTree *tree, string out_file_name) : fChain(NULL), fDebug(kDebug_Info + kDebug_Warning + kDebug_Error),
                                                      is_process(false), event(NULL), b_event(NULL),
                                                      output_file_name(out_file_name), b_field(0) {
    // Constructor. If tree is passed, initialize with that tree.
    if (tree) {
        Init(tree);
    } else {
        if (fDebug & kDebug_Warning)
            cout << "BaseAnna::BaseAna() Constructor called without a tree. Not initialized.\n";
    }

//  if (gROOT->IsProofServ()) {
//    if (gProofServ->IsMaster()) {
//      if(fDebug & kDebug_Info) printf("Macro running on master server\n");
//      // single remote init
//    } else {
//      if(fDebug & kDebug_Info)printf("Macro running on %d of %d\n", gProofServ->GetGroupId(),
//             gProofServ->GetGroupSize());
//      // parallel remote init
//    }
//  }else{
//    if(fDebug & kDebug_Info) cout << "Running outside of PROOF\n";
//  }

}

BaseAna::~BaseAna(){
    // Clean up. Warning!!! Careful here with the PROOF system. Deleting histograms causes errors.
}

void BaseAna::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
 
  if(fDebug & kDebug_L1) cout << "      BaseAna::Init(): \n";
 
  if(!event){
    event = new HpsEvent();
  }
  current_event=0;
  
  if(tree){
    fChain = tree;
    //
    // Determine the type of Chain, and if it is within Proof, since
    // unfortunately, the behavior of GetEntry() must depend on the type of chain.
    //
    if(false ){ // gROOT->IsProofServ()
      fChain_type = kIs_TProof;
      if(fDebug & kDebug_L1)   cout << "            -- PROOF type chain. \n";
    }else{
      if(fChain->InheritsFrom(TChain::Class())){
        if(is_process){
          fChain_type=kIs_TChain_Process;
          if(fDebug & kDebug_L1) cout << "            -- CHAIN::Process type chain. \n";
        }else{
          fChain_type=kIs_TChain;
          Notify();
          if(fDebug & kDebug_L1) cout << "            -- CHAIN type chain. \n";
          fChain->GetEntry(0);
        }
      }else{
        fChain_type=kIs_TTree;
        Notify();
        if(fDebug & kDebug_L1) cout << "            -- TREE type chain. \n";
      }
      
    }
   
  }else{
    if(fDebug & kDebug_Error) cout << "**** ERROR **** - init without a chain.\n";
  }
}

Bool_t BaseAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  if(fDebug & kDebug_L1) cout << "BaseAna::Notify():\n";
  if( !fChain ) cout << "Notify does not have a chain!!!!!\n";
  else{
    b_event = fChain->GetBranch("Event");
    if(fDebug & kDebug_L2) cout << "Notify(): got branch? " << (b_event!=NULL? "Yes" : "NO" ) << endl;
    
    if(b_event != NULL){
      if(!event){
        if(fDebug & kDebug_Error) cout << "ERROR ***** Notify has no event object.\n";
        return kFALSE;
      }
      b_event->SetAddress(&event);  // This is needed when running through a chain, not really for Proof but does no harm.
      if(fDebug & kDebug_L2) cout << "Notify(): Address set. \n";
      
    }else{
      if(fDebug & kDebug_Error) cout << "**** ERROR **** File in TChain does not contain Event branch: "<< fChain->GetTreeNumber() <<"\n";
      return kFALSE;
    }
  }
  return kTRUE;
}

void BaseAna::Begin(TTree *tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  if(fDebug & kDebug_L1) cout << "BaseAna::Begin(): \n";
  
  if(tree){
    if(!fChain){
      if(fDebug & kDebug_L1) cout << "      -- Begin called with tree and not already initialized, calling Init(tree) with is_process. \n";
      is_process=true;
      Init(tree);  // With Proof, we get instantiated with a tree, with a TChain, we may not.
    }
  }
}

void BaseAna::SlaveBegin(TTree */*tree */)
{
  // The SlaveBegin() function is called after the Begin() function when processing a chain.
  // When running with PROOF SlaveBegin() is called on each slave server, but Begin is NOT called on each slave,
  // so SlaveBegin() needs to setup the essentials, i.e. Histograms etc.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  if(fDebug & kDebug_L1)cout << "BaseAna::SlaveBegin(): \n";

  if(!event){
    if(fDebug & kDebug_L2)cout << "BaseAna::SlaveBegin(): create a new event object.\n";
    event = new HpsEvent();
  }
}

Int_t BaseAna::GetEntry(Long64_t entry, Int_t getall){
  // Retrieve the entry from the chain.
  // This routine simply gets the entry, it does not process the entry for you. It will check to see if the entry is in the current Tree
  // and if not link to the correct tree and call Notify()
  current_event=entry;
  if(fChain){
    if(fChain_type == kIs_TProof || fChain_type == kIs_TChain_Process ){ // TProof or TChain::Process()
      return (int)(fChain->GetTree() ? fChain->GetTree()->GetEntry(entry, getall) : 0);
    }else if(fChain_type == kIs_TChain){              // Standard Chain, needs to get correct tree first.
       Int_t treeNumInChain = fChain->GetTreeNumber();
      Long64_t centry = fChain->LoadTree(entry);
      if(centry < 0){
        if(fDebug & kDebug_Warning) cout << "Centry was negative. \n";
        return( -1 );
      }
      Int_t currentTreeNumInChain = fChain->GetTreeNumber();
      if (treeNumInChain != currentTreeNumInChain) {
        if(fDebug & kDebug_L1) cout << "Old Tree: "<< treeNumInChain << "  New Tree from the chain: " << currentTreeNumInChain << endl;
        Notify();
      }
      return (int)(fChain->GetTree() ? fChain->GetTree()->GetEntry(centry, getall) : 0);
      
    }else{ // Normal TTree, just get the entry.
      return (int)fChain->GetEntry(entry,getall);
    }
  }
  
  return 0;
}

Bool_t BaseAna::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  if(fDebug & kDebug_L2) cout << "Process: " << entry << " Tree: " << fChain->GetTreeNumber() << endl;
  int stat =GetEntry(entry);
  if(  stat <= 0 ){
    if(fDebug & kDebug_Error){
      cout << "GetEntry("<< entry << ") returned with status "<<  stat << endl;
      printf("i: %9ld  event: %9d  clust: %2d  track: %2d \n", evt_count, GetEventNumber(), GetNumberOfEcalClusters(), GetNumberOfTracks());
    }
    Abort("Bad event");
    return false;
  }
  if( (evt_count++ % fCounter_Freq ) == 0 && fDebug>0 ) {
    printf("i: %9ld  event: %9d  clust: %2d  track: %2d  tree: %2d\n", evt_count, GetEventNumber(), GetNumberOfEcalClusters(), GetNumberOfTracks(),fChain->GetTreeNumber());
  }
  return(kTRUE);
}


Int_t BaseAna::Next(void){
  // Get the next event and process it.
  int stat=Process( current_event); // Note that current_event is updated in GetEntry(), so ++current_event is NOT correct!
  current_event++;
  return(stat);
}

void BaseAna::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  if(fDebug & kDebug_L1) cout << "BaseAna::SlaveTerminate() \n";
}

//
// Helper method to allow itterative writing of nested lists to a file.
// It is expected that the file has already been opened.
//
void BaseAna::WriteList(TList *ll){
  TIter nn(ll);
  while(TObject *o=nn()){
    if(o->IsFolder()){
      if(fDebug & kDebug_L1) cout << "Create folder: "<< o->GetName() << " \n";
      output_file->mkdir(o->GetName());
      output_file->cd(o->GetName());
      if( o->IsA() == TList::Class() ){
        WriteList( (TList *)o);
      }else if(o->IsA() == TDirectory::Class()){
        TDirectory *dir=(TDirectory *)o;
        TList *ld = dir->GetList();
        WriteList(ld);
      }else{
        cout << "ERROR -- I don't know this type of folder object :" << o->ClassName() << endl;
      }
      output_file->cd("..");
    }else{
      if(fDebug & kDebug_L1) cout << "Write Object: "<< o->GetName() << " \n";
      o->Write();
    }
  }
}

void BaseAna::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  if(fDebug & kDebug_L1) cout << "BaseAna::Terminate() \n";
  TList *list=GetOutputList();
  if(fDebug & kDebug_L1) cout << "The list has "<< list->GetEntries() << " entries \n";
  
  // Write the contends to a file/
  if(fDebug & kDebug_L1) cout << "Writing output file: " << output_file_name << endl;
  if( !output_file) output_file = new TFile(output_file_name.data(),"RECREATE");

  WriteList(list);
  
  output_file->Close();
  
}

long BaseAna::Run(int nevent){
  // Run through nevents, if nevent <=0 then run until the end of the chain.
  long max_event = fChain->GetEntries();
  if( nevent>0 && nevent+current_event < max_event ) max_event = nevent;
  
  int stat=1;
  while(current_event+1 < max_event && stat>0){
    stat=Next();
  }
  return(current_event);
}

void BaseAna::PrintParticle(const HpsParticle *part,int i) const{
  if(part){
    // Identify which cluster is associated with this particle.
    int cluster_num=-2;
    TRefArray *cl_ref = part->getClusters();
    if(cl_ref->GetEntries()>0){
      EcalCluster *p_clus = static_cast<EcalCluster *>(cl_ref->At(0));
      for(int nc=0;nc<GetNumberOfEcalClusters();++nc){
        if( GetEcalCluster(nc) == p_clus){
          cluster_num = nc;
        }
      }
    }
    // Identify which track is associated with this particle.
    int track_num=-1;
    TRefArray *tr_ref = part->getTracks();
    if(tr_ref->GetEntries()>0){
      SvtTrack *p_track = static_cast<SvtTrack *>(tr_ref->At(0));
      if( part->getType() <= 32){
        for(int nt=0;nt<GetNumberOfTracks();++nt){
          if( GetTrack(nt) == p_track){
            track_num = nt;
            break;
          }
        }
      }else{
        for(int ng=0;ng< GetNumberOfGblTracks();++ng){
          if(GetGblTrack(ng) == p_track){
            track_num = 100+ng;
            break;
          }
        }
      }
    }
    vector<double> mom=part->getMomentum();
    printf("%2d: PID: %4d (%8.1f) Charge: %2d Type: %2d -- E= %6.4f P=(%5.2f,%5.2f,%5.2f) Cluster: %d  Track:%3d\n",
           i,part->getPDG(),part->getGoodnessOfPID(),part->getCharge(),part->getType(),
           part->getEnergy(),mom[0],mom[1],mom[2],cluster_num+1,track_num);
    std::vector<float> covmat  = part->getCovariantMatrix();
    printf("    CovMat: ");
    for(int i=0;i<10;++i) printf("%5.2f ",covmat[i]);
    printf("\n");
  }
}

void BaseAna::Print(Option_t *options) const{
  
  string opts(options);
  char evt_time[71];
  time_t evt_time_n;
  long evt_time_frac;
  if(event == NULL){
    cout <<  "Event is NULL \n";
    return;
  }
  
  long evt_time_raw = GetEventTime();
  if( evt_time_raw > 1000000000 ){
    evt_time_n = evt_time_raw/1000000000;
    evt_time_frac =evt_time_raw%1000000000;
  }else{
    evt_time_n = 0;
    evt_time_raw = GetEventTime();
    evt_time_frac=0;
  }
  bool do_all = opts.find("all")!=std::string::npos;

  struct tm *ti = localtime(&evt_time_n);
  strftime(evt_time,49,"%c",ti);
  
  // Header is always printed
  printf("Run: %6d  Event: %7d  time: %ld = %s %ld evt_in_chain: %ld\n",GetRunNumber(),GetEventNumber(),evt_time_n,evt_time,evt_time_frac,current_event);
  
  if(do_all || (opts.find("trig")!=std::string::npos)){
    printf("Trigger");
    if( IsPulserTrigger() )printf(" pulser ");
    if( IsSingle0Trigger())printf(" single0 ");
    if( IsSingle1Trigger())printf(" single1 ");
    if( IsPair0Trigger() )printf(" pair0 ");
    if( IsPair1Trigger() )printf(" pair1 ");
    printf("\n");
  
    printf("SVT status: ");
    if(IsSvtBiasOn() && IsSvtClosed() && IsSvtLatencyGood()) printf(" OK \n");
    else printf("NOT OK: Bias %1d Closed: %1d Latency: %1d \n",IsSvtBiasOn(),IsSvtClosed(),IsSvtLatencyGood());
    printf("N_Ecal:%3d  N_SVT:%3d N_Track:%3d \n",GetNumberOfEcalClusters(),GetNumberOfSvtHits(),GetNumberOfTracks());
    printf("N_Part:%2d  UC_V0:%2d UC_VC:%2d BC_V0:%2d  TC_V0:%2d  OE:%2d\n",
           GetNumberOfParticles(HpsParticle::FINAL_STATE_PARTICLE),
           GetNumberOfParticles(HpsParticle::UC_V0_CANDIDATE),
           GetNumberOfParticles(HpsParticle::UC_VC_CANDIDATE),
           GetNumberOfParticles(HpsParticle::BSC_V0_CANDIDATE),
           GetNumberOfParticles(HpsParticle::TC_V0_CANDIDATE),
           GetNumberOfParticles(HpsParticle::OTHER_ELECTRONS));
  }
  
  if(do_all || (opts.find("vert")!=std::string::npos)){
    printf("Constrained Particles:");
    for(auto particle_type: {HpsParticle::UC_V0_CANDIDATE,HpsParticle::UC_VC_CANDIDATE,HpsParticle::BSC_V0_CANDIDATE,HpsParticle::TC_V0_CANDIDATE}){
      if(GetNumberOfParticles(particle_type)){
        if(particle_type == HpsParticle::UC_V0_CANDIDATE)       cout << "Unconstrained Vertex: ";
        else if(particle_type == HpsParticle::UC_VC_CANDIDATE)  cout << "Unconst same side Vx: ";
        else if(particle_type == HpsParticle::BSC_V0_CANDIDATE) cout << "Beamspot Constrained: ";
        else if(particle_type == HpsParticle::TC_V0_CANDIDATE)  cout << "Target Constrained:   ";
      }
      
      for(int i=0; i<GetNumberOfParticles(particle_type); ++i){
        HpsParticle *part= GetParticle(particle_type, i);
        if(part){
          vector<double> vert= part->getVertexPosition();
          printf("Mass: %6.4f  Energy: %6.3f  Vertex:(%6.3f,%6.3f,%6.3f) V.chi2: %8.2f \n",
                 part->getMass(),part->getEnergy(),vert[0],vert[1],vert[2],part->getVertexFitChi2());
          TRefArray *daughter_refs = part->getParticles();
          for(int j=0; j<daughter_refs->GetEntries(); ++j){
            PrintParticle(static_cast<HpsParticle *>(daughter_refs->At(j)),j);
          }
        }
      }
    }
  }
  
  
  if(do_all || (opts.find("part")!=std::string::npos)){
    for(auto particle_type: {HpsParticle::FINAL_STATE_PARTICLE,HpsParticle::OTHER_ELECTRONS}){
      if(GetNumberOfParticles(particle_type)){
        if(particle_type == HpsParticle::FINAL_STATE_PARTICLE) cout << "Final State Particles: \n";
        else if(particle_type == HpsParticle::OTHER_ELECTRONS) cout << "Other Electrons: \n";
      }
      
      for(int i=0; i<GetNumberOfParticles(particle_type); ++i){
        HpsParticle *part= GetParticle(particle_type, i);
        PrintParticle(part,i);
      }
    }
  }
  
  if(do_all || (opts.find("track")!=std::string::npos)){
    cout << "Tracks: \n";
    for(int nt=0;nt< GetNumberOfTracks();++nt){
      SvtTrack *trk = GetTrack(nt);
      vector<double> pos=trk->getPositionAtEcal();
      vector<double> mom=trk->getMomentum(b_field);
      printf("%2d: Ch:%2d Type: %2d X2: %6.2f P:(%5.2f,%5.2f,%5.2f) L:(%5.2f,%5.2f,%5.2f) Time: %5.3f \n",nt,trk->getCharge(),trk->getType(),trk->getChi2(),mom[0],mom[1],mom[2],pos[0],pos[1],pos[2],trk->getTrackTime());
      printf("    CovMat: ");
      vector<float> covmat = trk->getCovariantMatrix();
      for(int i=0;i<15;++i) printf("%5.2f ",covmat[i]);
      printf("\n");
    }
    cout << "GBL Tracks: \n";
    for(int ng=0;ng<GetNumberOfGblTracks();++ng){
      GblTrack *gtrk = GetGblTrack(ng);
      vector<double> pos=gtrk->getPositionAtEcal();
      vector<double> mom=gtrk->getMomentum(b_field);
      printf("%2d: Ch:%2d Type: %2d X2: %6.2f P:(%5.2f,%5.2f,%5.2f) L:(%5.2f,%5.2f,%5.2f) Time: %5.3f \n",ng,gtrk->getCharge(),gtrk->getType(),gtrk->getChi2(),mom[0],mom[1],mom[2],pos[0],pos[1],pos[2],gtrk->getTrackTime());
    }
  }

  if(do_all || (opts.find("ecal")!=std::string::npos)){
    cout << "Clusters: \n";
    for(int nc=0;nc<GetNumberOfEcalClusters();++nc){
      EcalCluster *clus = GetEcalCluster(nc);
      vector<double> pos=clus->getPosition();
      EcalHit *seed=clus->getSeed();
      printf("%2d: Energy: %5.2f L:(%7.2f,%6.2f,%6.1f) T: %6.2f Seed: (%3d,%3d) hits:",nc+1,clus->getEnergy(),pos[0],pos[1],pos[2],clus->getClusterTime(),TrXhit(seed->getXCrystalIndex()),seed->getYCrystalIndex());
      TRefArray *ec_hits = clus->getEcalHits();
      for(int nh=0;nh<ec_hits->GetEntries();++nh){
        EcalHit *hit = (EcalHit *)ec_hits->At(nh);
        printf("(%3d,%3d) ",TrXhit(hit->getXCrystalIndex()),hit->getYCrystalIndex());
      }
      printf("\n");
    }
    cout << "Hits: \n";
    for(int nh=0;nh<GetNumberOfEcalHits();++nh){
      EcalHit *hit=GetEcalHit(nh);
      printf("%2d: E:%7.2f  T:%7.2f  (%3d,%3d) \n",nh,hit->getEnergy(),hit->getTime(),TrXhit(hit->getXCrystalIndex()),hit->getYCrystalIndex());
    }
  }

  if(do_all || (opts.find("svt")!=std::string::npos)){
    cout << "SVT -- Implement me!\n";
  }
}

HpsParticle    *BaseAna::GetFSParticle(int n){
  // Returns the n-th FINAL_STATE_PARTICLE from the current event.
  if(n >= GetNumberOfParticles(HpsParticle::FINAL_STATE_PARTICLE)){
    cout << "Asking for particle that is not in list. ";
    return(nullptr);
  }else{
    return GetParticle(HpsParticle::FINAL_STATE_PARTICLE,n);
  }
}

double BaseAna::GetAbsMomentum(HpsParticle *part){
  // Return the absolute of the momentum for the particle.
  vector<double> mom = part->getMomentum();
  return( sqrt( mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2] ));
}

double BaseAna::GetFSAbsMomentum(int n){
  HpsParticle *part = GetFSParticle(n);
  if(part) return GetAbsMomentum(part);
  else return -999;
}

void BaseAna::DrawEcal(int n){
  // Draw a representation of the ECAL with he hits for the current event.
  // Argument n=2 sets the line width.
  //
  TH2F *hit_map = EcalHitMap();
  TH2F *clusters =ClusterMap();
  
  FancyPlot(hit_map,0);
  clusters->SetLineWidth(n);
  clusters->Draw("box same");
  clusters->Draw("text same)");
  
}

TH2F *BaseAna::EcalHitMap(void){
  // Return a 2D histogram showing the distribution of ECal hits for this event.
  // Use:
  // ... // define ana
  // TH2F my_hitmap(ana->EcalHitMap());
  // my_hitmap.Draw()
  //
  TH2F *hit_map =(TH2F *)gROOT->FindObject("ecal_hit_map");
  if(hit_map == nullptr){
    hit_map =  new TH2F("ecal_hit_map","Ecal Hit Map",(ecal_nx+1)*2,-ecal_nx-0.5,ecal_nx+2-0.5,(ecal_ny+1)*2+1,-ecal_ny-1.5,ecal_ny+1.5);
  }else{
    hit_map->Reset();
  }
  
  for(int nhit=0;nhit<GetNumberOfEcalHits();++nhit){
    EcalHit *hit= GetEcalHit(nhit);
    int xhit=hit->getXCrystalIndex();
    int yhit=hit->getYCrystalIndex();
//    xhit = (xhit>=0?xhit:xhit+1);
    hit_map->Fill(TrXhit(xhit),yhit,hit->getEnergy());
  }
  return(hit_map);

}

TH2F *BaseAna::ClusterMap(int n_cl){
  // Return a 2D histogram showing the distribution of ECal clusters for this event.
  // Use:
  // ... // define ana
  // TH2F my_hitmap(ana->EcalHitMap());
  // my_hitmap.Draw()
  //

  TH2F *cluster_map =(TH2F *)gROOT->FindObject("cluster_hit_map");
  if(cluster_map == nullptr){
    cluster_map =  new TH2F("cluster_map","Ecal Hit Map",(ecal_nx+1)*2,-ecal_nx-0.5,ecal_nx+2-0.5,(ecal_ny+1)*2+1,-ecal_ny-1.5,ecal_ny+1.5);
  }else{
    cluster_map->Reset();
  }

  if(n_cl > GetNumberOfEcalClusters()) n_cl=0;
  string name("Cluster Map");
  if(n_cl>0) name.append(" ").append(to_string(n_cl));
  cluster_map->SetName(name.data());
  
  for(int nc=(n_cl==0?0:n_cl-1);nc<(n_cl!=0?n_cl:GetNumberOfEcalClusters());++nc){
    if(gDebug ) cout << "Histo for ECAL Cluster" << nc << endl;
    EcalCluster *clus= GetEcalCluster(nc);
    // EcalHit *seed = clus->getSeed();
    TRefArray *r_hits = clus->getEcalHits();
    for(int nh=0;nh<r_hits->GetEntries();++nh){
      EcalHit *hit=(EcalHit *)r_hits->At(nh);
      int xhit=hit->getXCrystalIndex();
      int yhit=hit->getYCrystalIndex();
//      xhit = (xhit>=0?trxhit:xhit+1);
      cluster_map->Fill(TrXhit(xhit),yhit,float(nc+1));
    }
  }

  cluster_map->SetMaximum(1.);
  cluster_map->SetLineColor(kRed);
  return(cluster_map);
}

void BaseAna::FancyPlot(TH2F *histo,int opt){
  //
  // Fancy plot of the Calorimeter.  For the TEST RUN. (NO LEAD GLASS)
  //
  //this defines the position of the top right region of big boxes, others will fall in by symmetry
  int ecal_x_first=1;
  int ecal_nx_first=0;
  int ecal_y_first=1;
  int ecal_ny_first=-1;
  
  int ecal_nx=23;
  int ecal_ny=5;
  
  TAxis *xax=histo->GetXaxis();
  TAxis *yax=histo->GetYaxis();
    
  TH2F *ones_lb;
  
  if (!gROOT->FindObject(TString(histo->GetName()).Append("_oneslb"))){//if this one exists all the rest probably exist too
    ones_lb=new TH2F(TString(histo->GetName()).Append("_oneslb"),"oneslb",(ecal_nx+1)*2,-ecal_nx-0.5,ecal_nx+2-0.5,(ecal_ny+1)*2+1,-ecal_ny-1.5,ecal_ny+1.5);
  }
  else {
    ones_lb= (TH2F*) gROOT->FindObject(TString(histo->GetName()).Append("_oneslb"));
    ones_lb->Clear();
  }
  
  if(!(opt& 0x4) ) histo->SetMaximum();
  if(histo->GetMaximum() < 1){
    histo->SetMaximum(1.1);
  }
  
  float SetMax=histo->GetMaximum();
  if(SetMax<1.1)SetMax=1.1;
  
  xax=ones_lb->GetXaxis();
  yax=ones_lb->GetYaxis();
  
  //this chunk of code just puts the grid in the right place
  for (int i=0;i<ecal_nx;i++){
    for (int j=0;j<ecal_ny;j++){
      
      ones_lb->SetBinContent(xax->FindBin(ecal_x_first+i),yax->FindBin(ecal_y_first+j),1);
      ones_lb->SetBinContent(xax->FindBin(ecal_x_first+i),yax->FindBin(ecal_ny_first-j),1);
      if(j==0 && i>0 && i<10){
        ;
      }else{
        ones_lb->SetBinContent(xax->FindBin(ecal_nx_first-i),yax->FindBin(ecal_ny_first-j),1);
        ones_lb->SetBinContent(xax->FindBin(ecal_nx_first-i),yax->FindBin(ecal_y_first+j),1);
      }
    }
  }
  
  ones_lb->Scale(SetMax);//scale them so the boxes are big enough
  
  //draw stuff
  
  if (opt & 0x2){
    
    if(opt&0x1){
      histo->Draw("Acol");
    }else{
      histo->Draw("Acolz");
    }
    ones_lb->Draw("boxsame");
    
    double xtal_size = 6.665;
    double xlow = ( xax->GetFirst() -1 + xax->GetXmin() );
    double xhi    ( xax->GetLast() - xax->GetFirst() + 1 + xax->GetXmin());
    double ylow = ( yax->GetFirst() -1 + yax->GetXmin() );
    double yhi  = ( yax->GetLast() - yax->GetFirst() + 1 + yax->GetXmin());
    
    double xlow2 = xtal_size * (xlow-1);
    double xhi2  = xtal_size * (xhi -1);
    double ylow2 = xtal_size * ylow;
    double yhi2  = xtal_size * yhi;
    
    TGaxis *lefta = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),gPad->GetUxmin(),gPad->GetUymax(),ylow2,yhi2,50510);
    TGaxis *bota  = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymin(),xlow2,xhi2,50510);
    lefta->SetName("leftaxis");
    lefta->SetTitle("y [mm]");
    bota->SetName("bottomaxis");
    bota->SetTitle("x [mm]");
    lefta->Draw();
    bota->Draw();
    
    if(opt & 0x4){
      TGaxis *righta = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),ylow,yhi,510,"");
      TGaxis *topa   = new TGaxis(gPad->GetUxmin(),gPad->GetUymax(),gPad->GetUxmax(),gPad->GetUymax(),xlow,xhi,510,"-");
      righta->SetName("rightaxis");
      righta->SetTitle("y crystal index");
      topa->SetName("topaxis");
      topa->SetTitle("x crystal index");
      //      righta->Draw();
      topa->Draw();
    }
  }else{
    xax->SetTitle("x crystal index");
    yax->SetTitle("y crystal index");
    if(opt&0x1){
      histo->Draw("col");
    }else{
      histo->Draw("colz");
    }
    ones_lb->Draw("boxsame");
  }
}


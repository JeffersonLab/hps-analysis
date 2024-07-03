//
// Created by Maurik Holtrop on 7/17/20.
//

#include "Moller.h"

Moller::Moller(std::string files) {
    if (!files.empty()) {
        ch.Add(files.c_str());
    } else {
        return;
    }
}



void Moller::Setup() {
    ///
    /// Setup a standard analysis.
    ///
    /// The analysis here is the same as the 2016 Bumphunt as documented in the analysis note.
    ///
    std::cout << "Setting up version " << __MOLLER_VERSION__ << "\n";
     // Step 1 - Select only events with 2 and only 2 electrons.
}

void Moller::Process() {
    ///
    /// Trigger Processing of all the filter steps and print a summary.
    ///
    if( filters.size() == 0){
        std::cout << "Please run Setup() first.\n";
        return;
    }
    std::cout << "Start processing.\n";
    TStopwatch timer;
    timer.Start();
    Print();
    timer.Stop();
    std::cout << "Process time was " << timer.RealTime() << " s  CPU: " << timer.CpuTime() << " s \n";
}

RNode Moller::Pass(RNode in, bool exact) {
   // Cut to select events that have two electrons. If exact = true (default) then cut on two and only
   // two electrons.
   // Parameters:  in     - Input RNode (RDateFrame)
   //              exact  - If true select exactly 2 electrons.

   return in.Filter("int n; for(int i=0;i<part_pdg.size();++i){if(part_pdg[i]==11) n++;} return n==2;");
}


RNode Moller::Cut_2_electrons(RNode in, bool exact) {
   // Cut to select events that have two electrons. If exact = true (default) then cut on two and only
   // two electrons.
   // Parameters:  in     - Input RNode (RDateFrame)
   //              exact  - If true select exactly 2 electrons.

   auto cut_two_electrons = [exact](RVec<int> part_pdg)->bool {
      int n=0;
      for(int i=0;i<part_pdg.size();++i){
         if(part_pdg[i]==11) n++;
      }
      return (exact && n==2) || (!exact && n>=2);
   };

   return in.Filter(cut_two_electrons,{"part_pdg"});
}

void Moller::Print(Option_t *opt) {

}
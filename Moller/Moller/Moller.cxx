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
     // Step 1 - Select only events with the pair1 trigger == true.
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

void Moller::Print(Option_t *opt) {

}
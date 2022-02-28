# HPS Analysis

A collection of tools for analysis of the HPS data sets, which do not already have a repository location elsewhere.

## MiniDst

[MiniDst](https://github.com/JeffersonLab/hps-analysis/tree/master/MiniDst)

* Sub directory that contains the code to create simple "flat" NTuple type ROOT files from LCIO or 2015/2016 DST files.

The MiniDst code creates a simple tuple style Data Summary Tape (DST) for HPS data from LCIO files. LCIO files
are the files produced by the Java based event reconstruction code. 

The advantages of a tuple style DST are that the data can be easily expanded, or reduced, in event size without 
breaking any existing code. It also makes the data easier to explore, with TTree::Draw() and the TBrowser. 
This type of data storage is also much easier to work with if you want to use the very powerful 
[RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) analysis
framework from ROOT.

Detailed instructions on this code are in the
[Wiki for hps-analysis](https://github.com/JeffersonLab/hps-analysis/wiki)

You can find some examples of how to use this code in the 
[examples](https://github.com/JeffersonLab/hps-analysis/tree/master/MiniDst/Examples) 
directory.

## Tools with other repo locations:
* [hps-java](https://github.com/JeffersonLab/hps-java)
    * The reconstruction code for HPS EVIO data to LCIO reconscructed files.
* [hpstr](https://github.com/JeffersonLab/hpstr)
    * The "Hipster" collection of tools to convert LCIO files to ROOT files and analyze them.
    * hpstr writes C++ class based ROOT files and also contains tools to work with them.
    * Contains converters for LCIO to ROOT.
    * Contains analysis codes for hpstr ROOT files.
    * contains post-root file analysis: bumphunter.
* [EvioTool](https://github.com/JeffersonLab/EvioTool)
    * Code for analyzing EVIO files directly from the ROOT prompt.
    * Also permits writing programs to analyze EVIO files directly.
    * Contains the HPS_Trigger_Filter program to filter EVIO files.
* [hps-trigger](https://github.com/JeffersonLab/hps-trigger)
    * Code related to analyzing the HPS trigger.
    
# HPS Analysis

A collection of tools for analysis of the HPS data sets, which do not already have a repository location elsewhere.

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
    
## Tools in this reo:
* [MiniDst](https://github.com/JeffersonLab/hps-analysis/tree/master/MiniDst)
    * Sub directory that contains code to create simple "flat" NTuple type ROOT files from LCIO or 2015/2016 DST files.
 
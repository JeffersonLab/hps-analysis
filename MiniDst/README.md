# MiniDst

This code creates a simple MiniDst for HPS data from LCIO files.

## Mini-, Micro-, Nano-DSTs
The idea is to create easy to use ROOT files that have a very simple structure. 
The basic structure of the ROOT files is a set of basic values (int, double), 
columns of values (vector<int> or vector<double>) or, for a few items, 
columns of columns (vector<vector<int> >,vector<vector<double> >). 

This design makes the ROOT files easy to read with simple TTree or RDataFrame commands, yet still
allows for a full, sophisticated analysis. Having such a simple structure makes is very easy to
trim the data set to only those variables you actually need in your analysis, saving space and time.
Code that uses such trimmed (micro-, nano-) DSTs are still compatible with ROOT files that 
have a larger subset of data in them. Code that requires a larger set of variables can make use of
the friend tree mechanism to add additional information without the need to re-write the base files.

A quick example of how to make a "pico-DST" from a larger one (eg mini-DST):

```c++
using namespace ROOT;
TChain ch("MiniDST");
ch.Add("minidst_*.root"); // Load all the available minidst root files.
RDataFrame df(ch);        // Open the TChain in an RDataFrame
df.Snapshot("NanoDST","nano_dst.root", // Write a select set of items to a nano dst.
        {"v0_mass", "v0_px", "v0_py", "v0_pz", "v0_vertex_x", "v0_vertex_y", "v0_vertex_z"})
```
Note that the output goes into a single file, which can be quite large if you had a lot of input files!

More details on RDataFrame are found in the ROOT manual pages: [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)

### Event Structure:

The event structure details are found in the file 
[MiniDst.h](https://github.com/JeffersonLab/hps-analysis/blob/master/MiniDst/MiniDst/MiniDst.h)
where the TTree Branches are set and given a name. A quick way to get the available columns is to ask the RDataFrame:

```c++
using namespace ROOT;
TChain ch("MiniDST");
ch.Add("minidst_*.root"); // Load all the available minidst root files.
RDataFrame df(ch);        // Open the TChain in an RDataFrame
auto column_list = df.GetColumnNames(); // Get a list of all the column names.
```

There are only a few items in the list that are vectors of vectors. These are, for instance, the indexes from a track to
the hits on the track, where you can have N tracks in the event with M(n) hits each, requiring a 2-dimensional array to
store the integers.
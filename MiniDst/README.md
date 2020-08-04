# MiniDst

This code creates a simple MiniDst for HPS data from LCIO files.

## Mini - Micro - Nano - Pico
The idea is to create easy to use ROOT files that have a very simple structure. 
The basic structure of the ROOT files is a set of basic values (int, double), 
columns of values (vector<int> or vector<double>) or, for a few items, 
columns of columns (vector<vector<int> >,vector<vector<double> >). 

This design makes the ROOT files easy to read with simple TTree or RDataFrame commands, yet still
allows for a full, sophisticated analysis. Having such a simple structure makes is very easy to
trim the data set to only those variables you actually need in your analysis, saving space and time.
Code that uses such trimmed (micro-, nano- or pico-) DSTs are still compatible with ROOT files that 
have a larger subset of data in them. Code that requires a larger set of variables can make use of
the friend tree mechanism to add additional information without the need to re-write the base files.

### Event Structure:

The event structure details are found in the file 
[MiniDst.h](https://github.com/JeffersonLab/hps-analysis/blob/master/MiniDst/MiniDst/MiniDst.h)
where the TTree Branches are set and given a name. 


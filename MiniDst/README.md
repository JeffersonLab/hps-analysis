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

## How to build MiniDst

#### Required Software:

* cmake version 3
    * Available for all Linux and Mac versions. On some systems you need to specify "cmake3" instead of just "cmake".
    Check your version with "cmake --version".
* c++17 compiler: gcc 7+, clang
    * I recommend you use at least gcc 9.3 (on MacOS, clang 11), which has far better support for lambda functions.
    * Ubuntu 20.04 LTS comes with gcc 9.3 pre-installed.
* ROOT version 6.18 or better.
    * The RDataFrame component of ROOT is rather new and still seeing active development and bug fixes. I would 
    recommend using 6.22, or the master branch from git.
    * You need to compile ROOT with the C++17 option: -Dcxx17=ON
* LCIO
    * Since we are converting LCIO to ROOT we need to be able to read the LCIO. 
    * The original repository of LCIO is at their GitHub repository: [LCIO](https://github.com/iLCSoft/LCIO.git) and
    contains detailed the installation instructions. Unfortunately, as of August 2020, the master and v02-14-x versions 
    of LCIO have very serious memory leaks for simple event loops, making this branch of the code useless. 
    * HPS has been using the older v02.07.05 version of LCIO, because that one is still Java compatible. For this code, 
    you need a c++17 capable version of LCIO. You can find Omar's version 2.7.5 of LCIO
    which is made c++17 compatible: [Omar's LCIO](https://github.com/LDMX-Software/lcio).
    A fork of that version is also MacOS compatible, and _attempts_ to clean _all_ the issues with c++17 compilation
    is also available at [Maurik's LCIO](https://github.com/mholtrop/LCIO.git). You should use the default branch:
    v02-07-05-cxx17.
    
    
#### At JLAB

At JLab, there is a 9.2.0 version of gcc available using "module", however there is no pre-installed ROOT 
version that is C++17 compatible. This situation is currently resolved by building a "master" version (6.23.99) 
of ROOT with 9.2.0 and C++17 compatibility. To use this:

```bash
module use /group/clas12/packages/local/etc/modulefiles
module use /apps/modulefiles
module load cmake/3.5.1
module load gcc/9.2.0
source /home/holtrop/root/bin/thisroot.sh
``` 

You will find pre-compiled libraries for LCIO and MiniDst in /home/holtrop/lib

### Building MiniDst

The building of MiniDST follows the standard cmake scheme. Below the recipe for building at JLab. 

* Setup your environment, as detailed above.
* Checkout the code from GitHub: 
```
git clone https://github.com/JeffersonLab/hps-analysis.git
```
* Go there and create a build directory for MiniDST:
```
cd hps-analysis/MiniDst
mkdir build
cd build
``` 
* Run cmake and compile the code:
```
cmake -DCMAKE_INSTALL_PREFIX=${HOME} -DLCIO_DIR=/home/holtrop ..
make -j 8
make install  # Install the code in your home directory's bin and lib.
```


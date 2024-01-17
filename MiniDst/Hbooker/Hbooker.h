//
// Created by Maurik Holtrop on 7/24/20.
//
// Class and code to write out HBook files.
// (Yes, old stuff.)
//
#ifndef HBOOKER_H
#define HBOOKER_H
#include <iostream>
#include <variant>
#include <chrono>
#include "TObject.h"
#include "MiniDst.h"

#include <stdlib.h>
#ifdef __CINT__
#define g77Fortran
#endif

#include "cfortran/cfortran.h"
#include "cfortran/packlib.h"

#define __Hbooker__Version__ "1.0.0"
using namespace std;

#define PAWC_SIZE 500000

typedef struct { float PAW[PAWC_SIZE]; } PAWC_DEF;
#define PAWC COMMON_BLOCK(PAWC,pawc)

typedef struct {
   int run_number;
   int event_number;
   int trigger;
   double rf_time1;
   double rf_time2;
} DATAC_DEF;
#define DATAC COMMON_BLOCK(DATAC,datac)


class Hbooker : public MiniDst{

   string hb_output_file_name;
   int hbook_id=10;
   int record_size=8191;
   int istat;

public:
   static string _version_() { return (__Hbooker__Version__); };
   Hbooker();
   explicit Hbooker(string fileout);
   void Start() override;
   void Process() override;
   void End() override;
   void Clear() override;
};

#endif

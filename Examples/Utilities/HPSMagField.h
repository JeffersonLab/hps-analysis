//
// Created by Maurik Holtrop on 2/28/23.
//

#ifndef HPS_ANALYSIS_HPSMAGFIELD_H
#define HPS_ANALYSIS_HPSMAGFIELD_H

#include <stack>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "TEveTrackPropagator.h"   // Contains the TEveMagField class.

class HPSMagField : public TEveMagField{

// protected:
public:
   std::vector<double> Bx,By,Bz;
   int nx=0,ny=0,nz=0;
   double x_min = 1e16, x_max = -1e16;
   double y_min = 1e16, y_max = -1e16;
   double z_min = 1e16, z_max = -1e16;
   double field_max = 0;

public:
   HPSMagField() = default;
   explicit HPSMagField(const std::string &field_map, double scale=1.);
   ~HPSMagField() override = default;
   [[nodiscard]] Double_t  GetMaxFieldMagD() const override { return field_max; }
   [[nodiscard]] TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const override;
   [[nodiscard]] double GetFieldBx(Double_t x, Double_t y, Double_t z) const;
   [[nodiscard]] double GetFieldBy(Double_t x, Double_t y, Double_t z) const;
   [[nodiscard]] double GetFieldBz(Double_t x, Double_t y, Double_t z) const;
   [[nodiscard]] int GetIndex(double x, double y, double z) const;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
ClassDef(HPSMagField,1);
#pragma clang diagnostic pop

};


#endif //HPS_ANALYSIS_HPSMAGFIELD_H

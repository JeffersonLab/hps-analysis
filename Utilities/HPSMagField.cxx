//
// Created by Maurik Holtrop on 2/28/23.
//

#include "HPSMagField.h"

ClassImp(HPSMagField);

HPSMagField::HPSMagField(const std::string &field_map, double scale) :TEveMagField() {
   // Read the field map file and fill the class vectors.

   std::ifstream in(field_map);
   if(!in.is_open()){
      std::cerr << "Could not open field map file: " << field_map << std::endl;
      exit(1);
   }
   std::istringstream iss;
   std::string line;
   std::getline(in,line);
   std::getline(in,line);
   iss.str(line);
   if(!(iss >> nx >> ny >> nz)){
      std::cerr << "Could not read the number of points in the field map." << std::endl;
      return;
   }

   int count = 0;
   while(std::getline(in,line)){
      int comment_number;
      std::string comment;
      iss.clear();
      iss.str(line);
      iss >> comment_number >> comment;
//      if(!(iss >> comment_number >> comment)) {
//         std::cout << "No longer reading comments" << std::endl;
//         std::cout << "Last line: " << line << std::endl;
//         break; } // Error
      std::cout << line << " == Comment[" << comment_number << "] " << comment << std::endl;
      if(comment_number==0 || ++count >= 7) break;
   }
   while(std::getline(in,line)){
      double x,y,z,Bx,By,Bz;
      iss.clear();
      iss.str(line);
      if(!(iss >> x >> y >> z >> Bx >> By >> Bz)) {
         std::cout << "ERROR reading the fieldmap. Courageously continuing." << std::endl;
         break; } // Error
      if( std::abs(Bx)>field_max) field_max= abs(Bx*scale*1000);
      if( std::abs(By)>field_max) field_max= abs(By*scale*1000);
      if( std::abs(Bz)>field_max) field_max= abs(Bz*scale*1000);
      this->Bx.push_back(Bx*scale*1000);
      this->By.push_back(By*scale*1000);
      this->Bz.push_back(Bz*scale*1000);
      if(x<x_min) x_min=x;
      if(x>x_max) x_max=x;
      if(y<y_min) y_min=y;
      if(y>y_max) y_max=y;
      if(z<z_min) z_min=z;
      if(z>z_max) z_max=z;
   }
   std::cout << "Field map loaded   : " << field_map << std::endl;
   std::cout << "Field map dimension: " << nx << " x " << ny << " x " << nz << " points." << std::endl;
   std::cout << "Field map range    : " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << std::endl;
}

int HPSMagField::GetIndex(double x, double y, double z) const {
   int ix = (x-x_min)/(x_max-x_min)*(nx-1);
   if(ix<0) ix=0;
   if(ix>=nx) ix=nx-1;
   int iy = (y-y_min)/(y_max-y_min)*(ny-1);
   if(iy<0) iy=0;
   if(iy>=ny) iy=ny-1;
   int iz = (z-z_min)/(z_max-z_min)*(nz-1);
   if(iz<0) iz=0;
   if(iz>=nz) iz=nz-1;
   return ix + iy*nx + iz*nx*ny;
}
TEveVectorD HPSMagField::GetFieldD(Double_t x, Double_t y, Double_t z) const{
   // Return the field at the given point.
   int index = GetIndex(x,y,z);
   return TEveVectorD(Bx[index],By[index],Bz[index]);
}

double HPSMagField::GetFieldBx(Double_t x, Double_t y, Double_t z) const{
   // Return the field at the given point.
   int index = GetIndex(x,y,z);
   return Bx[index];
}

double HPSMagField::GetFieldBy(Double_t x, Double_t y, Double_t z) const{
   // Return the field at the given point.
   int index = GetIndex(x,y,z);
   return By[index];
}

double HPSMagField::GetFieldBz(Double_t x, Double_t y, Double_t z) const{
   // Return the field at the given point.
   int index = GetIndex(x,y,z);
   return Bz[index];
}

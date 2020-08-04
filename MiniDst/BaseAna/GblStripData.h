/**
 *
 * @author: 	Per Hansson Adrian <phansson@slac.stanford.edu>
 * @section institution
 * 				SLAC
 * @version:    v 0.1
 * @date:       February 3, 2014
 */

#ifndef _GBL_STRIP_DATA_H_
#define _GBL_STRIP_DATA_H_

//--- C++ ---//
#include <iostream>

//--- ROOT ---//
#include <TObject.h>
#include <TClonesArray.h>
#include <TVector3.h>

class GblStripData : public TObject {

 private:
  int m_id;
  double m_path3D;
  double m_path;
  TVector3 m_u;
  TVector3 m_v;
  TVector3 m_w;
  TVector3 m_tdir_global;
  double m_tphi;
  TVector3 m_tpos;
  double m_phi;
  double m_lambda;
  double m_umeas;
  double m_umeas_err;
  double m_ms_angle;

 public:
  GblStripData();
  virtual ~GblStripData();

  // Setters
  void SetId(const int& id) { m_id = id; }
  void SetPath(const double& val) { m_path = val; }
  void SetPath3D(const double& val) { m_path3D = val; }
  void SetU(const double& x, const double& y, const double& z) { 
    m_u.Clear();
    m_u.SetXYZ(x,y,z);
  }
  void SetV(const double& x, const double& y, const double& z) { 
    m_v.Clear();
    m_v.SetXYZ(x,y,z);
  }
  void SetW(const double& x, const double& y, const double& z) { 
    m_w.Clear();
    m_w.SetXYZ(x,y,z);
  }
  void SetGlobalTrackDir(const double& x, const double& y, const double& z) { 
    m_tdir_global.Clear();
    m_tdir_global.SetXYZ(x,y,z);
  }
  void SetPhi(const double& phi) {
    m_phi = phi;
  }
  void SetLambda(const double& lambda) {
    m_lambda = lambda;
  }
  void SetUmeas(const double& umeas) {
    m_umeas = umeas;
  }
  void SetTrackPos(const double& x, const double& y, const double& z) { 
    m_tpos.Clear();
    m_tpos.SetXYZ(x,y,z);
  }
  void SetUmeasErr(const double& umeas_err) {
    m_umeas_err = umeas_err;
  }
  void SetMSAngle(const double& angle) {
    m_ms_angle = angle;
  }

  // Getters
  int GetId() const { return m_id; }
  TVector3 GetU() const { return m_u; }
  TVector3 GetV() const { return m_v; }
  TVector3 GetW() const { return m_w; }
  TVector3 GetGlobalTrackDir() const { return m_tdir_global; }  
  double GetPath() const { return m_path;}
  double GetPath3D() const { return m_path3D;}
  double GetPhi() const { return m_phi;}
  double GetLambda() const { return m_lambda;}
  double GetUmeas() const { return m_umeas;}
  TVector3 GetTrackPos() const { return m_tpos; }  
  double GetUmeasErr() const { return m_umeas_err;}
  double GetMSAngle() const { return m_ms_angle;}

  // General
  std::string toString() const;
  

  ClassDef(GblStripData,1) // Strip information needed by GBL

};



#endif

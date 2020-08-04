/**
 *
 * @author: 	phansson@slac.stanford.edu
 * @section institution
 * 				SLAC
 * @version:    v 0.1
 * @date:       February 3, 2014
 */

//-- DST --//
#include "GblStripData.h"

//-- C++ --//
#include <sstream>

ClassImp(GblStripData)

GblStripData::GblStripData()
	: 	TObject()
{}


GblStripData::~GblStripData()
{
}

std::string GblStripData::toString() const {
    std::ostringstream oss;
    oss << "GblStripData:\n";
    oss << "\t" << "Id        " << GetId() << "\n";
    oss << "\t" << "path      " << GetPath() << "\n";
    oss << "\t" << "path3D    " << GetPath3D() << "\n";
    oss << "\t" << "u         " << "(" << GetU().x() << "," <<  GetU().y() << "," <<  GetU().z() << ")" << "\n";
    oss << "\t" << "v         " << "(" << GetV().x() << "," <<  GetV().y() << "," <<  GetV().z() << ")"  << "\n";
    oss << "\t" << "w         " << "(" << GetW().x() << "," <<  GetW().y() << "," <<  GetW().z() << ")"  << "\n";
    oss << "\t" << "umeas     " << GetUmeas() << " +- " << GetUmeasErr() << "\n";
    oss << "\t" << "tpos      " << "(" << GetTrackPos().x() << "," << GetTrackPos().y() << "," << GetTrackPos().z() << ")"  << "\n";
    oss << "\t" << "tdir (gl) " << "(" << GetGlobalTrackDir().x() << "," << GetGlobalTrackDir().y() << "," << GetGlobalTrackDir().z() << ")"  << "\n";
    oss << "\t" << "track phi " << GetPhi() << "\n";
    oss << "\t" << "track lambda " << GetLambda() << "\n";
    oss << "\t" << "MS angle  " << GetMSAngle() << "\n";
    return oss.str();
  }

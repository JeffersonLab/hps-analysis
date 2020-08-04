/**
 * @section purpose: Stores track information needed by GBL
 * @author: 	Per Hansson Adrian <phansson@slac.stanford.edu>
 * @date:       February 3, 2014
 */

#include "GblTrackData.h"

ClassImp(GblTrackData)

GblTrackData::GblTrackData() 
	: TObject(), seed_track(NULL), m_gbl_strip_hits(new TRefArray()),
	  m_prjPerToCl(3,3), n_gbl_strip_hits(0), m_kappa(0), m_theta(0), m_phi(0),
      m_d0(0), m_z0(0) 
{}

GblTrackData::~GblTrackData() {
	delete m_gbl_strip_hits;
}

void GblTrackData::Clear(Option_t* /* option */) {
	TObject::Clear();
	m_gbl_strip_hits->Delete();
	n_gbl_strip_hits = 0;
}

void GblTrackData::addStrip(GblStripData* strip) {
	++n_gbl_strip_hits;
	m_gbl_strip_hits->Add(strip);
}

void GblTrackData::setPrjPerToCl(const unsigned int row, 
								 const unsigned int col,
								 const double val) {
	m_prjPerToCl(row, col) = val;
}

double GblTrackData::getKappa() const {
	// Note: Omega and Kappa both refer to the track
	//		 curvature.
	return (static_cast<SvtTrack*>(seed_track.GetObject())->getOmega());
}

double GblTrackData::getTheta() const {
	double tan_lambda = static_cast<SvtTrack*>(seed_track.GetObject())->getTanLambda();
	return TMath::PiOver2() - atan(tan_lambda);
}

double GblTrackData::getPhi() const {
	return (static_cast<SvtTrack*>(seed_track.GetObject())->getPhi0());
}

double GblTrackData::getD0() const {
	return (static_cast<SvtTrack*>(seed_track.GetObject())->getD0());
}

double GblTrackData::getZ0() const {
	return (static_cast<SvtTrack*>(seed_track.GetObject())->getZ0());
}

std::string GblTrackData::toString() const {
	std::ostringstream oss;
	oss << "GblTrackData:\n";
	oss << "\t" << "kappa " << getKappa() << "\n";
	oss << "\t" << "theta " << getTheta() << "\n";
	oss << "\t" << "phi   " << getPhi() << "\n";
	oss << "\t" << "d0    " << getD0() << "\n";
	oss << "\t" << "z0    " << getZ0() << "\n";
	oss << "\t" << n_gbl_strip_hits << " GBL strips:\n";
	TIter iter(m_gbl_strip_hits->MakeIterator());
	const GblStripData* strip = NULL;
	while( (strip = static_cast<GblStripData*>(iter())) != NULL) {
		oss << strip->toString() << "\n";
	}
	return oss.str();
}

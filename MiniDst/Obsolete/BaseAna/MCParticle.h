/**
 * @file MCParticle.h
 * @brief Class which implements an MC particle that stores information about
 *        tracks from the simulation.
 * @author: Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _MC_PARTICLE_H_
#define _MC_PARTICLE_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <map>
#include <string>
#include <vector>

//----------//
//   ROOT   //
//----------//
#include "TObject.h"
#include "TRefArray.h"

namespace hpsdst {

class MCParticle : public TObject {

public:

    /** Constructor */
    MCParticle();

    /** Destructor */
    virtual ~MCParticle();

    /**
     * Clear the data in this object.
     */
    void Clear(Option_t *option = "");

    /**
     * Print out information of this object.
     */
    void Print(Option_t *option = "") const;

    /**
     * Get the energy of the particle [GeV].
     *
     * @return The energy of the particle.
     */
    double getEnergy() const { return energy_; }

    /**
     * Get the PDG code of the particle.
     *
     * @return The PDG code of the particle.
     */
    int getPdgID() const { return pdg_id_; }

    /**
     * Get the generator status of the particle.
     *
     * @return The generator status.
     */
    int getGenStatus() const { return gen_status_; }

    /**
     * Get the global time of the particle's creation [ns].
     *
     * @return The global time of the particle's creation.
     */
    double getTime() const { return time_; }

    /**
     * Get the XYZ vertex of the particle's creation [mm].
     *
     * @return The vertex of the particle.
     */
    std::vector<double> getVertex() const { return {x_, y_, z_}; }


    /**
     * Get the endpoint of the particle where it was destroyed
     * or left the detector [mm].
     *
     * @return The endpoint of the particle
     */
    std::vector<double> getEndPoint() const { return {end_x_, end_y_, end_z_}; }

    /**
     * Get the XYZ momentum of the particle [MeV].
     *
     * @return The momentum of the particle.
     */
    std::vector<double> getMomentum() const { return {px_, py_, pz_}; }

    /**
     * Get the mass of the particle [GeV].
     *
     * @return The mass of the particle.
     */
    double getMass() const { return mass_; }

    /**
     * The charge of the particle (units of electron charge).
     *
     * @return The charge of the particle.
     */
    double getCharge() const { return charge_; }

    /**
     * Get the number of daughter particles.
     */
    int getDaughterCount() const {
        return daughters_->GetEntriesFast();
    }

    /**
     * Get a daughter particle by index.
     *
     * @param daughter_index The index of the daughter particle.
     */
    MCParticle *getDaughter(int daughter_index) {
        return static_cast<MCParticle *>(daughters_->At(daughter_index));
    }

    /**
     * Get the number of parent particles.
     *
     * @return The number of parent particles.
     */
    int getParentCount() const {
        return parents_->GetEntriesFast();
    }

    /**
     * Get a parent particle by index.
     *
     * @param parent_index The index of the parent particle.
     */
    MCParticle *getParent(int parent_index) {
        return static_cast<MCParticle *>(parents_->At(parent_index));
    }

    /**
     * Set the energy of the particle [MeV].
     *
     * @param energy The energy of the particle.
     */
    void setEnergy(const double energy) { this->energy_ = energy; }

    /**
     * Set the PDG code of the hit.
     *
     * @param Pdg_id The PDG code of the hit.
     */
    void setPdgID(const int Pdg_id) { pdg_id_ = Pdg_id; }

    /**
     * Set the generator status of the hit.
     *
     * @param Get_status The generator status of the hit.
     */
    void setGenStatus(const int Get_status) { gen_status_ = Get_status; }

    /**
     * Set the global time of the particle's creation [ns].
     *
     * @param Time The global time of the particle's creation.
     */
    void setTime(const double Time) { time_ = Time; }

    /**
     * Set the vertex of the particle [mm].
     *
     * @param X The vertex X position.
     * @param Y The vertex Y position.
     * @param Z The vertex Z position.
     */
    void setVertex(const double X, const double Y, const double Z) {
        x_ = X;
        y_ = Y;
        z_ = Z;
    }

    /**
     * Set the end point of the particle [mm].
     *
     * @param End_x The X end point.
     * @param End_y The Y end point.
     * @param End_z The Z end point.
     */
    void setEndPoint(const double End_x, const double End_y, const double End_z) {
        end_x_ = End_x;
        end_y_ = End_y;
        end_z_ = End_z;
    }

    /**
     * Set the momentum of the particle [GeV].
     *
     * @param Px The X momentum.
     * @param Py The Y momentum.
     * @param Pz The Z momentum.
     */
    void setMomentum(const double Px, const double Py, const double Pz) {
        px_ = Px;
        py_ = Py;
        pz_ = Pz;
    }

    /**
     * Set the mass of the particle [GeV].
     *
     * @param Mass The mass of the particle.
     */
    void setMass(const double Mass) { mass_ = Mass; }

    /**
     * Set the charge of the particle.
     *
     * @param Charge The charge of the particle.
     */
    void setCharge(const double Charge) { charge_ = Charge; }

    /**
     * Add a daughter particle.
     *
     * @param daughter The daughter particle.
     */
    void addDaughter(MCParticle *daughter) { daughters_->Add(daughter); }

    /**
     * Add a parent particle.
     *
     * @param parent The parent particle.
     */
    void addParent(MCParticle *parent) { parents_->Add(parent); }

public:

    /** The energy of the particle. */
    double energy_{0};

    /** The PDG code of the particle. */
    int pdg_id_{0};

    /** The generator status. */
    int gen_status_{-1};

    /** The global creation time. */
    double time_{0};

    /** The X vertex. */
    double x_{0};

    /** The Y vertex. */
    double y_{0};

    /** The Z vertex. */
    double z_{0};

    /** The X end point. */
    double end_x_{0};

    /** The Y end point. */
    double end_y_{0};

    /** The Z end point. */
    double end_z_{0};

    /** The X momentum.*/
    double px_{0};

    /** The Y momentum. */
    double py_{0};

    /** The Z momentum. */
    double pz_{0};

    /** The particle's mass. */
    double mass_{0};

    /** The particle's charge. */
    double charge_{0};

    /** The list of daughter particles. */
    TRefArray *daughters_{new TRefArray};

    /** The list of parent particles. */
    TRefArray *parents_{new TRefArray};

    /**
     * ROOT class definition.
     */
    ClassDef(MCParticle, 2);

}; // MCParticle

}

using namespace hpsdst;

#endif // _MC_PARTICLE_H_

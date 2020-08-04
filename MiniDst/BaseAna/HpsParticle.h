/**
 *
 * @file HpsParticle.h
 * @brief Class used to describe an HPS particle.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date March 29, 2013
 *
 */

#ifndef __HPS_PARTICLE_H__
#define __HPS_PARTICLE_H__

//------------//
//--- ROOT ---//
//------------//
#include <TObject.h>
#include <TClonesArray.h>
#include <TRef.h>
#include <TRefArray.h>
#include <TMatrixD.h>

/** Forward declarations of HPS event classes */
class SvtTrack;
class EcalCluster;

class HpsParticle : public TObject { 
    
public:
    
    /** Enum constants used to denote the different particle types */
    enum ParticleType {
        FINAL_STATE_PARTICLE = 0,
        UC_V0_CANDIDATE      = 1,
        BSC_V0_CANDIDATE     = 2,
        TC_V0_CANDIDATE      = 3,
        UC_MOLLER_CANDIDATE  = 4,
        BSC_MOLLER_CANDIDATE = 5,
        TC_MOLLER_CANDIDATE  = 6,
        OTHER_ELECTRONS      = 7,
        UC_VC_CANDIDATE      = 8
    };
    
    /** Default Constructor. */
    HpsParticle();
    
    /**
     * Copy Constructor.
     *
     * @param particle_obj HpsParticle that will be cloned
     */
    HpsParticle(const HpsParticle &particle_obj);
    
    /** Destructor. */
    ~HpsParticle();
    
    /**
     * Copy assignment operator.
     *
     * @param particle_obj HpsParticle object that will be copied
     * @return A copy of the given HpsParticle object
     */
    HpsParticle &operator=(const HpsParticle &particle_obj);
    
    /** Reset this HpsParticle object */
    void Clear(Option_t *option="");
    
    /**
     * Add a reference to an SvtTrack object.  This will be used to
     * associate a particle with the SVT track that composes it.
     *
     * @param svt_track SVT track whose reference will be added
     */
    void addTrack(SvtTrack* svt_track);
    
    /**
     * Add a reference to an EcalCluster object.  This will be used to
     * associated a particle with the Ecal cluster that composes it.
     *
     * @param ecal_cluster Ecal cluster whose reference will be added
     */
    void addCluster(EcalCluster* ecal_cluster);
    
    /**
     * Add a reference to an HpsParticle object.  This will be used to
     * add daughter particles to this particle.
     *
     * @param particle Daughter particle composing this particle
     */
    void addParticle(HpsParticle* particle);
    
    /**
     * Set the charge of the particle.
     *
     * @param Charge Particle charge
     */
    void setCharge(const int Charge) { charge = Charge; };
    
    /**
     * Set the 1/2 covariant matrix for the track
     *
     * @parm An array of 15 doubles, which is the 1/2 covariant matrix.
     */
    void setCovMatrix(const float *covariant_matrix){
        memcpy(&covmatrix,covariant_matrix,sizeof(covmatrix));
    }
    
    /**
     * Set the overall goodness of the PID for this particle.
     *
     * @param Goodness_pid The goodness of the PID.
     */
    void setGoodnessOfPID(const double Goodness_pid) { goodness_pid = Goodness_pid; };
    
    /**
     * Set the type of this particle.
     *
     * @param Type The type of this particle
     */
    void setType(const int Type) { type = Type; };
    
    /**
     * Set the PDG ID of this particle.
     *
     * @param Pdg The PDG ID of this particle
     */
    void setPDG(const int Pdg) { pdg = Pdg; };
    
    /**
     * Set the energy of the particle in GeV.
     *
     * @param Energy The energy of this particle
     */
    void setEnergy(const double Energy) { energy = Energy; };
    
    /**
     * Set the invariant mass of the particle in GeV.
     *
     * @param Mass The invariant mass of the particle
     */
    void setMass(const double Mass) { mass = Mass; };
    
    /**
     * Set the momentum of the particle in GeV.
     *
     * @param momentum An array containing the three momentum components
     *                 of the particle.
     */
    void setMomentum(const double* momentum);
    
    /**
     * Set the corrected momentum of the paritcle in GeV.
     *
     * @param momentum An array containing the three momentum components
     *                 of the particle.
     */
    void setCorrMomentum(const double* momentum);
    
    /**
     * Set the vertex position of the particle.
     *
     * @param vtx_pos An array containing the three vertex position
     *                components of the particle
     */
    void setVertexPosition(const float* vtx_pos);
    
    /**
     * Set the chi^2 of the vertex fit.
     *
     * @param Vtx_fit_chi2 The chi^2 of the vertex fit
     */
    void setVertexFitChi2(const double Vtx_fit_chi2) { vtx_fit_chi2 = Vtx_fit_chi2; };
    
    /**
     * Get the SVT tracks associated with this particle.
     *
     * @return An array of references to the tracks associated with this
     *         particle
     */
    TRefArray* getTracks() const;
    
    /**
     * Get the Ecal clusters associated with this particle.
     *
     * @return An array of references to the Ecal clusters associated
     *         with this particle
     */
    TRefArray* getClusters() const;
    
    /**
     * Get the daughter particles composing this particle.
     *
     * @return An array of references to the daughter particles associated
     *         with this particle
     */
    TRefArray* getParticles() const;
    
    /**
     * Get the charge of the particle.
     *
     * @return The particle charge
     */
    int getCharge() const { return charge; };
    
    /**
     * Get the 1/2 covariant matrix as vector
     *
     * @return 1/2 covariant matrix in a vector<float> of 15 entries.
     */
    std::vector<float> getCovariantMatrix() const{
        std::vector<float> retvec(sizeof(covmatrix));
        retvec.assign(covmatrix,covmatrix+sizeof(covmatrix));    // Copy data into vector.
        return(retvec);                           // C++11 returns matrix by move.
    };
    
    /**
     * Get the [rown,coln] element of the covariant matrix.
     *
     * @return covariant matrix [row,col] value (float)
     */
    float getCovariantMatrix(int i,int j) const{
        if(i<0 || j<0 || i>3 || j>3) return(0);
        // int k=0; for(int q=0;q<=j;++q) k+=q;
        if(j<i){
            int tmp=i;
            i = j;
            j = tmp;
        }
        int k = i + j + (j-1)*(j>1)+(j-2)*(j>2);
        return(covmatrix[k]);
    };

    /**
     * Get the 1/2 covariant matrix as TMatrixD
     * Note that this is MUCH less efficient than getCovariantMatrix or getCovariantMatrix(i,j)
     *
     * @return 1/2 covariant matrix in a 5x5 TMatrix.
     */
    TMatrixD getCovariantTMatrix() const{
        double data[16];
        int k=0;
        for(int i=0;i<4;i++)
            for(int j=0;j<i;j++){
                data[i*4+j] = covmatrix[k];
                data[j*4+i] = covmatrix[k];
                k++;
            }
        TMatrixD retmat(4,4,data);
        return(retmat);        // This copies the matrix twice! (ouch)
    };

    /**
     * Get the overall goodness of the PID for this particle.
     *
     * @return The goodness of the PID.
     */
    double getGoodnessOfPID() const { return goodness_pid; };
    
    /**
     * Get the type of this particle.
     *
     * @return The type of this particle.
     */
    int getType() const { return type; };
    
    /**
     * Get the particle ID.
     *
     * @return The particle ID
     */
    int getPDG() const { return pdg; };
    
    /**
     * Get the energy of the particle in GeV.
     *
     * @return The particle energy in GeV
     */
    double getEnergy() const { return energy; };
    
    /**
     * Get the invariant mass of the particle in GeV.
     *
     * @return The invariant mass of the particle in GeV
     */
    double getMass() const { return mass; };
    
    /**
     * Get the momentum of the particle in GeV.
     *
     * @return The momentum of the particle
     */
    std::vector<double> getMomentum() const;
    
    /**
     * Get the corrected momentum of the paritcle in GeV.
     *
     * @return The corrected momentum of the particle.
     */
    std::vector<double> getCorrMomentum() const;
    
    /**
     * Get the vertex position of the particle.
     *
     * @return The vertex position of the particle
     */
    std::vector<double> getVertexPosition() const;
    
    /**
     * Get the chi^2 of the vertex fit.
     *
     * @return The chi^2 of the vertex fit
     */
    double getVertexFitChi2() const { return vtx_fit_chi2; };
    
    ClassDef(HpsParticle, 3);
    
public:
    
    /** An array of references to tracks associated with this particle */
    TRefArray* svt_tracks;
    
    /**
     * An array of references to Ecal clusters associated with this
     * particle
     */
    TRefArray* ecal_clusters;
    
    /**
     *  An array of references to daughter particles associated with this
     *  particle
     */
    TRefArray* particles;
    
    /** The number of daughters associated with this particle */
    int n_daughters;
    
    /** The charge of this particle */
    int charge;
    
    /** The type of this particle */
    int type;
    
    /** The PDG ID of this particle */
    int pdg;
    
    /** The goodness of PID of this particle. */
    double goodness_pid;
    
    /** The x component of the momentum of this particle in GeV */
    double px;
    
    /** The x component of the corrected momentum of this particle in GeV */
    double px_corr;
    
    /** The y component of the momentum of this particle in GeV */
    double py;
    
    /** The y component of the corrected momentum of this particle in GeV */
    double py_corr;
    
    /** The z component of the momentum of this particle in GeV */
    double pz;
    
    /** The z component of the corrected momentum of this particle in GeV */
    double pz_corr;
    
    /** The x component of the vertex of this particle in mm*/
    double vtx_x;
    
    /** The y component of the vertex of this particle in mm */
    double vtx_y;
    
    /** The z component of the vertex of this particle in mm */
    double vtx_z;
    
    /** The chi^2 of the vertex fit */
    double vtx_fit_chi2;
    
    /** The 1/2 Covariant Matrix. This is the lower 1/2. **/
    float covmatrix[10];
    
    /** The energy of the particle in GeV */
    double energy;
    
    /** The invariant mass of the particle in GeV */
    double mass;
    
};  // HpsParticle

#endif // _HPS_PARTICLE_H_

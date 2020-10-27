//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_BODY_H
#define HDB5_IO_BODY_H

#include "MathUtils/MathUtils.h"
#include <vector>
#include <unordered_map>
#include <memory>

#include "Mask.h"
#include "Mesh.h"
#include "PoleResidue.h"

namespace HDB5_io {

  // Forward declaration
  class HydrodynamicDataBase;

  class WaveDrift;

  /**
  * \class Body
  * \brief Class for storing the added mass, radiation damping, excitation components, mesh, hydrostatic stiffness, etc.
  * All components are stored, even null DOF components, resulting in 6 DOF vectors and 6x6 matrices.
  */
  class Body {

   public:

    typedef std::unordered_map<unsigned int, std::shared_ptr<mathutils::LookupTable1D<double, mathutils::Vector6d<double>>>> HDBinterpolator;

    /// Constructor for a body
    /// \param id index of the body in the HDB
    /// \param name name of the body
    /// \param hdb pointer to the HDB
    Body(unsigned int id, const std::string &name, HydrodynamicDataBase *hdb);

//    void SetName(std::string name) {m_name = name;}

    /// Return the name of the BEM body database
    /// \return Name of the BEM body
    std::string GetName() const { return m_name; }

    /// Return the identifier of the BEM Body
    /// \return Identifier
    unsigned int GetID() const { return m_id; }

    //
    // Setters
    //

    /// Define the position of body stored in BEM body database
    /// \param position Position of the body
    void SetPosition(const mathutils::Vector3d<double> &position);

    /// Define the mask on the force components
    /// \param mask Mask on the force components
    void SetForceMask(mathutils::Vector6d<bool> mask);

    /// Set the complex matrix of the diffraction coefficient
    /// \param iangle Corresponding wave direction
    /// \param diffractionMatrix Complex matrix of the diffraction coefficient
    void SetDiffraction(unsigned int iangle, const Eigen::MatrixXcd &diffractionMatrix);

    /// Set the complex vector of the diffraction coefficient
    /// \param iangle Corresponding wave direction
    /// \param iw Corresponding wave frequency
    /// \param diffractionVector Complex vector of the diffraction coefficient
    void SetDiffraction(unsigned int iangle, unsigned int iw, const Eigen::VectorXcd &diffractionVector);

    /// Set the complex matrix of the Froude-Krylov coefficient
    /// \param iangle Corresponding wave direction
    /// \param froudeKrylovMatrix Complex matrix of the Froude-Krylov coefficient
    void SetFroudeKrylov(unsigned int iangle, const Eigen::MatrixXcd &froudeKrylovMatrix);

    /// Set the complex vector of the Froude-Krylov coefficient
    /// \param iangle Corresponding wave direction
    /// \param iw Corresponding wave frequency
    /// \param froudeKrylovVector Complex vector of the Froude-Krylov coefficient
    void SetFroudeKrylov(unsigned int iangle, unsigned int iw, const Eigen::VectorXcd &froudeKrylovMatrix);

    /// Set the complex matrix of the wave excitation coefficient
    /// \param iangle Corresponding wave direction
    /// \param excitationMatrix Complex matrix of the wave excitation coefficients
    void SetExcitation(unsigned int iangle, const Eigen::MatrixXcd &excitationMatrix);

    /// Set the complex vector of the wave excitation coefficient
    /// \param iangle Corresponding wave direction
    /// \param iw Corresponding wave frequency
    /// \param excitationVector Complex vector of the wave excitation coefficients
    void SetExcitation(unsigned int iangle, unsigned int iw, const Eigen::VectorXcd &excitationVector);

    /// Compute the excitation loads from the diffraction loads and the Froude-Krylov loads.
    void ComputeExcitation();

    /// Set the infinite added mass of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param CMInf Infinite added mass matrix
    void SetInfiniteAddedMass(Body *BodyMotion, const mathutils::Matrix66<double> &CMInf);

    /// Set the zero frequency added mass of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param CMInf zero frequency added mass matrix
    void SetZeroFreqAddedMass(Body *BodyMotion, const mathutils::Matrix66<double> &CMZero);

    /// Set the radiation mask of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param mask radiation mask of the BEM body with respect to the motion of another BEM body
//    void SetRadiationMask(Body *BodyMotion, const mathutils::Matrix66<bool> &mask);
    void SetRadiationMask(Body *BodyMotion, const mathutils::Matrix66<bool> &mask);

    /// Set the complex matrix of the response amplitude operator
    /// \param iangle Corresponding wave direction
    /// \param RAO Complex matrix of the response amplitude operator
    void SetRAO(unsigned int iangle, const Eigen::MatrixXcd &RAO);

    /// Set the complex vector of the response amplitude operator
    /// \param iangle Corresponding wave direction
    /// \param iw Corresponding wave frequency
    /// \param RAO Complex vector of the response amplitude operator
    void SetRAO(unsigned int iangle, unsigned int iw, const Eigen::VectorXcd &RAO);

    enum interpolatedData {
      IRF_K, IRF_KU, ADDED_MASS, RADIATION_DAMPING
    };

    void SetHDBInterpolator(interpolatedData type, Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listData);

    /// Set the hydrostatic stiffness Matrix
    /// \param hydrostaticStiffnessMatrix Hydrostatic stiffness matrix
    void SetStiffnessMatrix(const mathutils::Matrix33<double> &hydrostaticStiffnessMatrix);

    /// Set the hydrostatic stiffness Matrix
    /// \param hydrostaticStiffnessMatrix Hydrostatic stiffness matrix
    void SetStiffnessMatrix(const mathutils::Matrix66<double> &hydrostaticStiffnessMatrix);

    void SetInertia(const mathutils::Matrix66<double> &inertiaMatrix);

    void SetMooring(const mathutils::Matrix66<double> &mooringMatrix);

    void SetLinearDamping(const mathutils::Matrix66<double> &linearDampingMatrix);

    /// Load the mesh, from vertices and connectivity
    /// \param vertices vertices container
    /// \param faces connectivity of all faces
    void LoadMesh(const std::vector<mathutils::Vector3d<double>> &vertices, const std::vector<Eigen::VectorXi> &faces);

    /// Add a vector of modal coefficients for the 6 degrees of freedom of this body, generated by one degree of freedom
    /// of the BodyMotion body
    /// \param BodyMotion body at the origin of the perturbation
    /// \param modalCoefficients modal coefficients
    void AddModalCoefficients(Body *BodyMotion, std::vector<PoleResidue> modalCoefficients);

    //
    // Getters
    //

    /// Check if the RAO are calculated
    /// \return true if RAO were calculated
    bool HasRAO() const;

    bool HasModal(Body *BodyMotion) const;

    bool HasIRF() const;

    bool HasInertia() const;

    bool HasHydrostatic() const;

    bool HasMooring() const;

    bool HasDamping() const;

    bool HasZeroFreqAddedMass(Body *BodyMotion) const;

    /// Return the position of the body as stored in the HDB
    /// \return position of the body
    mathutils::Vector3d<double> GetPosition() const;

    /// Return the mask value applied on a specific motion mode
    /// \param iforce Index of force
    /// \return Mask on the force mode
    Mask GetForceMask() const;

    /// Get the mesh of the body
    /// \return mesh of the body
    Mesh *GetMesh() const;

#ifdef H5_USE_VTK

    /// Visualize the mesh of the body
    void VisualizeMesh() const;

#endif

    /// Get the hydrostatic stiffness of this body
    /// \return 3x3 matrix containing the hydrostatic components (heave, roll, pitch)
    mathutils::Matrix33<double> GetHydrostaticStiffnessMatrix() const;

    mathutils::Matrix66<double> GetInertiaMatrix() const;

    mathutils::Matrix66<double> GetMooringMatrix() const;

    mathutils::Matrix66<double> GetDampingMatrix() const;

    /// Get the diffraction components for this body
    /// \param iangle index of the angle
    /// \return matrix containing the diffraction component for the given angle
    Eigen::MatrixXcd GetDiffraction(unsigned int iangle) const;

    /// Get the diffraction components for this body
    /// \param iangle index of the angle
    /// \param iforce index of the force dof
    /// \return matrix containing the diffraction component for the given angle
    Eigen::VectorXcd GetDiffraction(unsigned int iangle, unsigned int iforce) const;

    /// Get the Froude-Krylov components for this body
    /// \param iangle index of the angle
    /// \return matrix containing the Froude-Krylov component for the given angle
    Eigen::MatrixXcd GetFroudeKrylov(unsigned int iangle) const;

    /// Get the Froude-Krylov components for this body
    /// \param iangle index of the angle
    /// \param iforce index of the force dof
    /// \return matrix containing the Froude-Krylov component for the given angle
    Eigen::VectorXcd GetFroudeKrylov(unsigned int iangle, unsigned int iforce) const;

    /// Get the excitation components for this body
    /// \param iangle index of the angle
    /// \return matrix containing the excitation component for the given angle
    Eigen::MatrixXcd GetExcitation(unsigned int iangle) const;

    /// Get the excitation components for this body
    /// \param iangle index of the angle
    /// \param iforce index of the force dof
    /// \return matrix containing the excitation component for the given angle
    Eigen::VectorXcd GetExcitation(unsigned int iangle, unsigned int iforce) const;

    /// Get the infinite added mass, resulting from a motion of body BodyMotion
    /// \param BodyMotion body which motion create added mass on this body
    /// \return 6x6 matrix added mass
    mathutils::Matrix66<double> GetInfiniteAddedMass(Body *BodyMotion) const;

    /// Get the zero frequency added mass, resulting from a motion of body BodyMotion
    /// \param BodyMotion body which motion create added mass on this body
    /// \return 6x6 matrix added mass
    mathutils::Matrix66<double> GetZeroFreqAddedMass(Body *BodyMotion) const;

    /// Get the radiation mask, between this body and BodyMotion body
    /// \param BodyMotion body
    /// \return radiation mask
    mathutils::Matrix66<bool> GetRadiationMask(Body *BodyMotion) const;

    /// Get the response amplitude operator for this body
    /// \param iangle index of the angle
    /// \return matrix containing the response amplitude operator for the given angle
    Eigen::MatrixXcd GetRAO(unsigned int iangle) const;

    /// Get the response amplitude operator for this body
    /// \param iangle index of the angle
    /// \param iforce index of the force dof
    /// \return matrix containing the response amplitude operator for the given angle and dof
    Eigen::VectorXcd GetRAO(unsigned int iangle, unsigned int iforce) const;

    /// Get the infinite added mass, resulting from a motion of this body
    /// \return 6x6 matrix added mass
    mathutils::Matrix66<double> GetSelfInfiniteAddedMass();

    /// Get the interpolator corresponding to the data type given
    /// \param type type of the data interpolated (IRF_K, IRF_KU, ADDED_MASS, RADIATION_DAMPING)
    /// \return interpolator
    HDBinterpolator *GetHDBInterpolator(interpolatedData type);

    /// Get the interpolated data for the following parameters
    /// \param type type of the data interpolated (IRF_K, IRF_KU, ADDED_MASS, RADIATION_DAMPING)
    /// \param BodyMotion body at the origin of the perturbation
    /// \param idof index of the degree of freedom considered
    /// \param frequencies set of frequencies for which the data are interpolated
    /// \return interpolated data in matrix form (6 x nfreq)
    Eigen::MatrixXd
    GetHDBInterpolatedData(interpolatedData type, Body *BodyMotion, unsigned int idof,
                           mathutils::VectorN<double> frequencies);

    /// Get the modal coefficients (poles and residues) for the 6DOFs
    /// \param BodyMotion body at the origin of the perturbation
    /// \param idof index of the degree of freedom at the origin of the perturbation
    /// \return vector of modal coefficients
    std::vector<PoleResidue> GetModalCoefficients(Body *BodyMotion, int idof);

    /// Get the modal coefficients (poles and residues) for the iforce dof, generated by the idof dof of BodyMotion body
    /// \param BodyMotion body at the origin of the perturbation
    /// \param idof index of the degree of freedom at the origin of the perturbation
    /// \param iforce index of the degree of freedom considered
    /// \return modal coefficients
    PoleResidue GetModalCoefficients(Body *BodyMotion, int idof, int iforce);

   protected:

    HydrodynamicDataBase *m_HDB;                   ///< HDB containing this data container
    unsigned int m_id;                             ///< ID of the BEM Body
    std::string m_name;                            ///< Name of the body
    mathutils::Vector3d<double> m_position;        ///< Position of the body COG

    Mask m_forceMask;                              ///< Mask applied on the force

    std::shared_ptr<Mesh> m_mesh;                  ///< mesh of the body

    // TODO :replace these std::vector with 2D-interpolator, or at least unordered_map ?

    /// Excitation loads for all wave directions (std::vector), all dof (rows of Eigen::MatrixXcd) and all wave frequencies (cols of Eigen::MatrixXcd).
    std::vector<Eigen::MatrixXcd> m_excitation;

    /// Froude-Krylov loads for all wave directions (std::vector), all dof (rows of Eigen::MatrixXcd) and all wave frequencies (cols of Eigen::MatrixXcd).
    std::vector<Eigen::MatrixXcd> m_froudeKrylov;

    /// Diffraction loads for all wave directions (std::vector), all dof (rows of Eigen::MatrixXcd) and all wave frequencies (cols of Eigen::MatrixXcd).
    std::vector<Eigen::MatrixXcd> m_diffraction;

    bool m_isRAO = false;

    /// Response Amplitude Operators for all wave directions (std::vector), all dof (rows of Eigen::MatrixXcd) and all wave frequencies (cols of Eigen::MatrixXcd).
    std::vector<Eigen::MatrixXcd> m_RAO;

    bool m_isHydrostatic = false;
    mathutils::Matrix33<double> m_hydrostaticStiffnessMatrix;   ///< Hydrostatic matrix
    bool m_isInertia = false;
    mathutils::Matrix66<double> m_inertia;         ///< Inertia matrix
    bool m_isMooring = false;
    mathutils::Matrix66<double> m_mooringStiffnessMatrix;         ///< Mooring stiffness matrix
    bool m_isDamping = false;
    mathutils::Matrix66<double> m_linearDampingMatrix;         ///< Linear damping matrix (not radiation).

    std::unordered_map<Body *, mathutils::Matrix66<bool>> m_radiationMask;          ///< Radiation mask
    std::unordered_map<Body *, mathutils::Matrix66<double>> m_infiniteAddedMass;    ///< Infinite added mass for each body
    std::unordered_map<Body *, mathutils::Matrix66<double>> m_zeroFreqAddedMass;    ///< Zero frequency added mass for each body

    std::unordered_map<Body *, std::vector<std::vector<PoleResidue>>> m_modalCoefficients;       ///<

    std::shared_ptr<HDBinterpolator> m_interpK;                     ///< Impulse response function interpolator
    std::shared_ptr<HDBinterpolator> m_interpKu;                    ///< Impulse response function speed dependent interpolator
    bool m_isIRF = false;

    std::shared_ptr<HDBinterpolator> m_addedMass;                   ///< added mass interpolator
    std::shared_ptr<HDBinterpolator> m_radiationDamping;            ///< radiation damping interpolator
//    std::vector<std::vector<mathutils::Interp1dLinear<double, std::complex<double>>>> m_waveDirInterpolators;   ///<

    /// Allocate the excitation containers
    /// \param nFrequencies number of frequencies
    /// \param nDirections number of wave directions
    void AllocateAll();

  };

} // namespace HDB5_io

#endif //HDB5_IO_BODY_H

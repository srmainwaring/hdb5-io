//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_BODY_H
#define HDB5_IO_BODY_H

#include "MathUtils/Matrix.h"
#include "MathUtils/Matrix33.h"
#include "MathUtils/Matrix66.h"
#include "MathUtils/Vector3d.h"
#include <vector>
#include <unordered_map>
#include <memory>

//#include "Mode.h"
#include "Mask.h"
#include "Mesh.h"

namespace HDB5_io {

  // Forward declaration
  class HydrodynamicDataBase;

  class Discretization1D;

  /**
  * \class Body
  * \brief Class for storing the added mass, radiation damping, excitation components, mesh, hydrostatic stiffness, etc.
  * All components are stored, even null DOF components, resulting in 6 DOF vectors and 6x6 matrices.
  */
  class Body {

   public:

    /// Constructor for a body
    /// \param id index of the body in the HDB
    /// \param name name of the body
    /// \param hdb pointer to the HDB
    Body(unsigned int id, const std::string &name, HydrodynamicDataBase *hdb);

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
    void SetForceMask(mathutils::Vector6d<int> mask);

    /// Define the mask on the DOF of the body
    /// \param mask Mask on the DOF of the body
    void SetMotionMask(mathutils::Vector6d<int> mask);

    /// Set the complex matrix of the diffraction coefficient
    /// \param iangle Corresponding wave direction
    /// \param diffractionMatrix Complex matrix of the diffraction coefficient
    void SetDiffraction(unsigned int iangle, const Eigen::MatrixXcd &diffractionMatrix);

    /// Set the complex matrix of the froude-krylov coefficient
    /// \param iangle Corresponding wave direction
    /// \param froudeKrylovMatrix Complex matrix of the diffraction coefficient
    void SetFroudeKrylov(unsigned int iangle, const Eigen::MatrixXcd &froudeKrylovMatrix);

    /// Set the complex matrix of the wave excitation coefficient
    /// \param iangle Corresponding wave direction
    /// \param excitationMatrix Complex matrix of the wave excitation coefficients
    void SetExcitation(unsigned int iangle, const Eigen::MatrixXcd &excitationMatrix);

    /// Compute the excitation loads from the diffraction loads and the Froude-Krylov loads.
    void ComputeExcitation();

    /// Set the infinite added mass of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param CMInf Infinite added mass matrix
    void SetInfiniteAddedMass(Body *BodyMotion, const mathutils::Matrix66<double> &CMInf);

    /// Set the radiation mask of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param mask radiation mask of the BEM body with respect to the motion of another BEM body
//    void SetRadiationMask(Body *BodyMotion, const mathutils::Matrix66<bool> &mask);
    void SetRadiationMask(Body *BodyMotion, const mathutils::Matrix66<int> &mask);

    /// Set the impulse response function of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotionBEM body to which the motion is considered
    /// \param listIRF List of impulse response function (size nforce x ntime) for each DOF
    void SetImpulseResponseFunctionK(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listIRF);

    /// Set the impulse response function (steady-speed dependent) of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param listIRF List of impulse response function (size nforce x ntime) for each DOF
    void SetImpulseResponseFunctionKu(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listIRF);

    /// Set the impulse response function of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotionBEM body to which the motion is considered
    /// \param listIRF List of impulse response function (size nforce x ntime) for each DOF
    void SetAddedMass(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listAddedMass);

    /// Set the impulse response function (steady-speed dependent) of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param listIRF List of impulse response function (size nforce x ntime) for each DOF
    void SetRadiationDamping(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listRadiationDamping);

    /// Set the hydrostatic stiffness Matrix
    /// \param hydrostaticStiffnessMatrix Hydrostatic stiffness matrix
    void SetStiffnessMatrix(const mathutils::Matrix33<double> &hydrostaticStiffnessMatrix);

    /// Set the hydrostatic stiffness Matrix
    /// \param hydrostaticStiffnessMatrix Hydrostatic stiffness matrix
    void SetStiffnessMatrix(const mathutils::Matrix66<double> &hydrostaticStiffnessMatrix);

    /// Load the mesh, from vertices and connectivity
    /// \param vertices vertices container
    /// \param faces connectivity of all faces
    void LoadMesh(const std::vector<mathutils::Vector3d<double>> &vertices, const std::vector<Eigen::VectorXi> &faces);

//    /// Set the wave drift coefficient database
//    void SetWaveDrift();

    //
    // Getters
    //

    /// Return the position of the body as stored in the HDB
    /// \return position of the body
    mathutils::Vector3d<double> GetPosition() const;

    /// Return the mask value applied on a specific motion mode
    /// \param imotion Index of motion
    /// \return Mask on the motion mode
    Mask GetMotionMask() const;

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

    /// Get the radiation mask, between this body and BodyMotion body
    /// \param BodyMotion body
    /// \return radiation mask
    mathutils::Matrix66<bool> GetRadiationMask(Body *BodyMotion) const;

    /// Get the infinite added mass, resulting from a motion of this body
    /// \return 6x6 matrix added mass
    mathutils::Matrix66<double> GetSelfInfiniteAddedMass();

    /// Get the interpolator for the impulse response function (IRF)
    /// \param BodyMotion body which motion is at the origin of the IRF
    /// \param idof index of the dof considered
    /// \return interpolator of the IRF
    mathutils::Interp1d<double, mathutils::Vector6d<double>> *GetIRFInterpolatorK(Body *BodyMotion, unsigned int idof);

    /// Get the interpolator for the impulse response function (IRF) with the advance speed corrections
    /// \param BodyMotion body which motion is at the origin of the IRF
    /// \param idof index of the dof considered
    /// \return interpolator of the IRF
    mathutils::Interp1d<double, mathutils::Vector6d<double>> *GetIRFInterpolatorKu(Body *BodyMotion, unsigned int idof);

    /// Get the added mass interpolator
    /// \param BodyMotion body which motion create added mass on this body
    /// \param idof index of the dof considered
    /// \return interpolator of the added mass
    mathutils::Interp1d<double, mathutils::Vector6d<double>> *
    GetAddedMassInterpolator(Body *BodyMotion, unsigned int idof);

    /// Get the radiation damping interpolator
    /// \param BodyMotion body which motion create radiation damping on this body
    /// \param idof index of the dof considered
    /// \return interpolator of the radiation damping
    mathutils::Interp1d<double, mathutils::Vector6d<double>> *
    GetRadiationDampingInterpolator(Body *BodyMotion, unsigned int idof);

    /// Get the components, for the given interpolator, corresponding to the frequencies
    /// \param interpolator interpolator containing the data
    /// \param frequencies frequencies to interpolate
    /// \return matrix component
    Eigen::MatrixXd
    GetMatrixComponentFromIterator(mathutils::Interp1d<double, mathutils::Vector6d<double>> *interpolator,
                                   Discretization1D frequencies);

//    std::shared_ptr<FrWaveDriftPolarData> GetWaveDrift() const;


   private:
    HydrodynamicDataBase *m_HDB;                   ///< HDB from which BEM data are extracted
    unsigned int m_id;                             ///< ID of the BEM Body
    std::string m_name;                            ///< Name of the body
    mathutils::Vector3d<double> m_position;        ///< Position of the body from HDB

    Mask m_forceMask;                              ///< Mask applied on the force
    Mask m_motionMask;                             ///< Mask applied on the DOF

    std::shared_ptr<Mesh> m_mesh;                  ///< mesh of the body

    std::vector<Eigen::MatrixXcd> m_excitation;    ///< Complex coefficient of the excitation force
    std::vector<Eigen::MatrixXcd> m_froudeKrylov;  ///< Complex coefficient of the froude-krylov force
    std::vector<Eigen::MatrixXcd> m_diffraction;   ///< Complex coefficient of the diffraction force

    mathutils::Matrix33<double> m_hydrostaticStiffnessMatrix;   ///< Hydrostatic matrix

    std::unordered_map<Body *, mathutils::Matrix66<bool>> m_radiationMask;          ///< Radiation mask
    std::unordered_map<Body *, mathutils::Matrix66<double>> m_infiniteAddedMass;    ///< Infinite added mass for each body

    typedef std::unordered_map<Body *, std::vector<std::shared_ptr<mathutils::Interp1d<double, mathutils::Vector6d<double>>> >> HDBinterpolator;
    HDBinterpolator m_interpK;                     ///< Impulse response function interpolator
    HDBinterpolator m_interpKu;                    ///< Impulse response function speed dependent interpolator

    HDBinterpolator m_addedMass;                   ///< added mass interpolator
    HDBinterpolator m_radiationDamping;            ///< radiation damping interpolator
//    std::vector<std::vector<mathutils::Interp1dLinear<double, std::complex<double>>>> m_waveDirInterpolators;   ///<

//    std::shared_ptr<FrWaveDriftPolarData> m_waveDrift;  ///< List of wave drift coefficients

    /// Allocate the excitation containers
    /// \param nFrequencies number of frequencies
    /// \param nDirections number of wave directions
    void AllocateAll(unsigned int nFrequencies, unsigned int nDirections);

  };

} // namespace HDB5_io

#endif //HDB5_IO_BODY_H
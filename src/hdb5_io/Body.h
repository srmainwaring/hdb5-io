//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_BODY_H
#define HDB5_IO_BODY_H

#include "MathUtils/MathUtils.h"
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

    /// Set the complex matrix of the response amplitude operator
    /// \param iangle Corresponding wave direction
    /// \param RAO Complex matrix of the response amplitude operator
    void SetRAO(unsigned int iangle, const Eigen::MatrixXcd &RAO);


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

    /// Load the mesh, from vertices and connectivity
    /// \param vertices vertices container
    /// \param faces connectivity of all faces
    void LoadMesh(const std::vector<mathutils::Vector3d<double>> &vertices, const std::vector<Eigen::VectorXi> &faces);

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

    bool HasRAO() const {
      return m_isRAO;
    }

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

    HDBinterpolator *GetHDBInterpolator(interpolatedData type);

    Eigen::MatrixXd
    GetHDBInterpolatedData(interpolatedData type, Body *BodyMotion, unsigned int idof, Discretization1D frequencies);


    std::shared_ptr<WaveDrift> GetWaveDrift() const;


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

    bool m_isRAO = false;
    std::vector<Eigen::MatrixXcd> m_RAO;           ///< response amplitude operators

    mathutils::Matrix33<double> m_hydrostaticStiffnessMatrix;   ///< Hydrostatic matrix

    std::unordered_map<Body *, mathutils::Matrix66<bool>> m_radiationMask;          ///< Radiation mask
    std::unordered_map<Body *, mathutils::Matrix66<double>> m_infiniteAddedMass;    ///< Infinite added mass for each body

    std::shared_ptr<HDBinterpolator> m_interpK;                     ///< Impulse response function interpolator
    std::shared_ptr<HDBinterpolator> m_interpKu;                    ///< Impulse response function speed dependent interpolator

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

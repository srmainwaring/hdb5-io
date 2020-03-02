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

namespace HDB5_io {

  class Body {

   public:

    Body(unsigned int id, const std::string& name) : m_id(id), m_name(name) {};

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
    void SetPosition(const mathutils::Vector3d<double>& position) { m_position = position; }

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
    void SetInfiniteAddedMass(Body *BodyMotion, const Eigen::MatrixXd &CMInf);

    /// Set the impulse response function of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotionBEM body to which the motion is considered
    /// \param listIRF List of impulse response function (size nforce x ntime) for each DOF
    void SetImpulseResponseFunctionK(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listIRF);

    /// Set the impulse response function (steady-speed dependent) of the BEM body with respect to the motion of another BEM body
    /// \param BodyMotion BEM body to which the motion is considered
    /// \param listIRF List of impulse response function (size nforce x ntime) for each DOF
    void SetImpulseResponseFunctionKu(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listIRF);

    /// Set the hydrostatic stiffness Matrix
    /// \param hydrostaticStiffnessMatrix Hydrostatic stiffness matrix
    void SetStiffnessMatrix(const mathutils::Matrix33<double> &hydrostaticStiffnessMatrix);

    /// Set the hydrostatic stiffness Matrix
    /// \param hydrostaticStiffnessMatrix Hydrostatic stiffness matrix
    void SetStiffnessMatrix(const mathutils::Matrix66<double> &hydrostaticStiffnessMatrix);

//    /// Set the wave drift coefficient database
//    void SetWaveDrift();

    //
    // Getters
    //

    /// Return the mask value applied on a specific motion mode
    /// \param imotion Index of motion
    /// \return Mask on the motion mode
    Mask GetMotionMask() const { return m_motionMask; }

    /// Return the mask value applied on a specific motion mode
    /// \param iforce Index of force
    /// \return Mask on the force mode
    Mask GetForceMask() const { return m_forceMask; }

    Eigen::MatrixXcd GetDiffraction(unsigned int iangle) const;

    Eigen::VectorXcd GetDiffraction(unsigned int iangle, unsigned int iforce) const;

    Eigen::MatrixXcd GetFroudeKrylov(unsigned int iangle) const;

    Eigen::VectorXcd GetFroudeKrylov(unsigned int iangle, unsigned int iforce) const;

    Eigen::MatrixXcd GetExcitation(unsigned int iangle) const;

    Eigen::VectorXcd GetExcitation(unsigned int iangle, unsigned int iforce) const;

    mathutils::Matrix66<double> GetInfiniteAddedMass(Body *BodyMotion) const;

    mathutils::Matrix66<double> GetSelfInfiniteAddedMass();

//    mathutils::Interp1d<double, mathutils::Vector6d<double>> *
//    GetIRFInterpolatorK(FrBEMBody *BodyMotion, unsigned int idof);
//
//    mathutils::Interp1d<double, mathutils::Vector6d<double>> *
//    GetIRFInterpolatorKu(FrBEMBody *BodyMotion, unsigned int idof);

    mathutils::Matrix33<double> GetHydrostaticStiffnessMatrix() const { return m_hydrostaticStiffnessMatrix; }

//    std::shared_ptr<FrWaveDriftPolarData> GetWaveDrift() const;

    void AllocateAll(unsigned int Nfrequencies, unsigned int nDirections);


   private:
    unsigned int m_id;                             ///< ID of the BEM Body
    std::string m_name;                            ///< Name of the body
    mathutils::Vector3d<double> m_position;        ///< Position of the body from HDB

    Mask m_forceMask;                              ///< Mask applied on the force
    Mask m_motionMask;                             ///< Mask applied on the DOF

//    std::vector<Mode> m_forceModes;                ///< List of activated force modes
//    std::vector<Mode> m_motionModes;               ///< List of activated motion modes

//    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_excitationMask;   ///<
    std::vector<Eigen::MatrixXcd> m_excitation;    ///< Complex coefficient of the excitation force
    std::vector<Eigen::MatrixXcd> m_froudeKrylov;  ///< Complex coefficient of the froude-krylov force
    std::vector<Eigen::MatrixXcd> m_diffraction;   ///< Complex coefficient of the diffraction force

    mathutils::Matrix33<double> m_hydrostaticStiffnessMatrix;   ///< Hydrostatic matrix

    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> m_radiationMask;   ///< Radiation mask
    std::unordered_map<Body *, mathutils::Matrix66<double>> m_infiniteAddedMass;    ///< Infinite added mass for each body
    std::unordered_map<Body *, std::vector<std::shared_ptr<mathutils::Interp1d<double, mathutils::Vector6d<double>>> >> m_interpK; ///< Impulse response function interpolator
    std::unordered_map<Body *, std::vector<std::shared_ptr<mathutils::Interp1d<double, mathutils::Vector6d<double>>> >> m_interpKu; ///< Impulse response function speed dependent interpolator

//    std::vector<std::vector<mathutils::Interp1dLinear<double, std::complex<double>>>> m_waveDirInterpolators;   ///<

//    std::shared_ptr<FrWaveDriftPolarData> m_waveDrift;  ///< List of wave drift coefficients

  };

} // namespace HDB5_io

#endif //HDB5_IO_BODY_H

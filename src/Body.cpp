//
// Created by lletourn on 26/02/20.
//

#include "Body.h"

namespace HDB5_io {


  //
  // Setters
  //

  void Body::SetForceMask(mathutils::Vector6d<int> mask) {
    m_forceMask.SetMask(mask);
  }

  void Body::SetMotionMask(mathutils::Vector6d<int> mask) {
    m_motionMask.SetMask(mask);
  }

  void Body::SetDiffraction(unsigned int iangle, const Eigen::MatrixXcd &diffractionMatrix) {
//    assert(iangle < GetNbWaveDirections());
//    assert(diffractionMatrix.rows() == 6);
//    assert(diffractionMatrix.cols() == GetNbFrequencies());
    m_diffraction[iangle] = diffractionMatrix;
  }

  void Body::SetFroudeKrylov(unsigned int iangle, const Eigen::MatrixXcd &froudeKrylovMatrix) {
//    assert(iangle < GetNbWaveDirections());
//    assert(froudeKrylovMatrix.rows() == 6);
//    assert(froudeKrylovMatrix.cols() == GetNbFrequencies());
    m_froudeKrylov[iangle] = froudeKrylovMatrix;
  }

  void Body::SetExcitation(unsigned int iangle, const Eigen::MatrixXcd &excitationMatrix) {
//    assert(iangle < GetNbWaveDirections());
//    assert(excitationMatrix.rows() == 6);
//    assert(excitationMatrix.cols() == GetNbFrequencies());
    m_excitation[iangle] = excitationMatrix;
  }

  void Body::ComputeExcitation() {

    /// This subroutine computes the excitation loads from the diffraction loads and the Froude-Krylov loads.

    assert(m_diffraction.size() == m_froudeKrylov.size());
//    for (unsigned int iangle = 0; iangle < GetNbWaveDirections(); ++iangle) {
    for (unsigned int iangle = 0; iangle < m_diffraction.size(); ++iangle) {
      m_excitation[iangle] = m_diffraction[iangle] + m_froudeKrylov[iangle];
    }
  }

  void Body::SetInfiniteAddedMass(Body *BodyMotion, const mathutils::Matrix66<double> &CMInf) {
    m_infiniteAddedMass[BodyMotion] = CMInf;
  }

  void Body::SetRadiationMask(HDB5_io::Body *BodyMotion, const mathutils::Matrix66<bool> &mask){
    m_radiationMask[BodyMotion] = mask;
  }

  void Body::SetStiffnessMatrix(const mathutils::Matrix33<double> &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix;
  }

  void Body::SetStiffnessMatrix(const mathutils::Matrix66<double> &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix.block<3, 3>(2, 2);
  }

  //
  // Getters
  //

  Eigen::MatrixXcd Body::GetDiffraction(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_diffraction[iangle];
  }

  Eigen::VectorXcd Body::GetDiffraction(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
//    assert(iforce < 6);
    return m_diffraction[iangle].row(iforce);
  }

  Eigen::MatrixXcd Body::GetFroudeKrylov(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_froudeKrylov[iangle];
  }

  Eigen::VectorXcd Body::GetFroudeKrylov(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
//    assert(iforce < 6);
    return m_froudeKrylov[iangle].row(iforce);
  }

  Eigen::MatrixXcd Body::GetExcitation(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_excitation[iangle];
  }

  Eigen::VectorXcd Body::GetExcitation(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
//    assert(iforce < 6);
    return m_excitation[iangle].row(iforce);
  }

  mathutils::Matrix66<double> Body::GetInfiniteAddedMass(Body *BodyMotion) const {
    return m_infiniteAddedMass.at(BodyMotion);
  }

  mathutils::Matrix66<double> Body::GetSelfInfiniteAddedMass() {
    return m_infiniteAddedMass[this];
  }

  void Body::AllocateAll(unsigned int nFrequencies, unsigned int nDirections) {

    // This subroutine allocates the arrays for the hdb.

    auto nbForce = GetForceMask().GetNbDOF();

    // --> Allocating arrays for excitations

//    auto nbWaveDir = GetNbWaveDirections();
    m_diffraction.reserve((unsigned long) nDirections);
    m_froudeKrylov.reserve((unsigned long) nDirections);
    m_excitation.reserve((unsigned long) nDirections);

//    auto nbFreq = GetNbFrequencies();
    for (int i = 0; i < nDirections; ++i) {
      Eigen::MatrixXcd mat(nbForce, nFrequencies);
      m_diffraction.push_back(mat);
      m_froudeKrylov.push_back(mat);
      m_excitation.push_back(mat);
    }

//    /// --> Allocating arrays for radiation
//    auto nbBodies = GetNbBodies();
//    m_radiationMask.reserve(nbBodies);
//
//    auto nbTime = GetNbTimeSamples();
//    for (unsigned int ibody = 0; ibody < nbBodies; ++ibody) {
//
//      auto body = m_HDB->GetBody(ibody);
//      auto nbMotion = GetMotionMask().GetNbDOF();
//
//      Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> mask(nbForce, nbMotion);
//      mask.setConstant(true);
//      m_radiationMask.push_back(mask);
//
//    }



  }

}
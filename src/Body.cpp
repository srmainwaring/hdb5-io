//
// Created by lletourn on 26/02/20.
//

#include "Body.h"

namespace HDB5_io {


  //
  // Setters
  //

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

  void Body::SetInfiniteAddedMass(Body *BodyMotion, const Eigen::MatrixXd &CMInf) {
    assert(CMInf.rows() == 6);
    assert(CMInf.cols() == 6);
    m_infiniteAddedMass[BodyMotion] = CMInf;
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

}
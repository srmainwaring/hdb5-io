//
// Created by pierre-yves on 13/01/2021.
//

#include "Kochin.h"
#include "HydrodynamicDataBase.h"

namespace HDB5_io {

  Kochin::Kochin(HydrodynamicDataBase *hdb, const double &kochin_step) : m_HDB(hdb), m_kochin_step(kochin_step){

    // Constructor of the class.

    // Number of angular steps.
    ComputeNbKochinAngles();

    // Allocations.
    auto nDirections = m_HDB->GetWaveDirectionDiscretization().size();
    auto nFrequencies = m_HDB->GetFrequencyDiscretization().size();
    m_kochin_diffraction.reserve((unsigned long) nDirections);
    m_kochin_diffraction_derivate.reserve((unsigned long) nDirections);

    for (int i = 0; i < nDirections; ++i) {
      Eigen::MatrixXd mat(m_nb_kochin_angle, nFrequencies);
      m_kochin_diffraction.push_back(mat);
      m_kochin_diffraction_derivate.push_back(mat);
    }

  }

  double Kochin::GetKochinStep() const {
    return m_kochin_step;
  }

  void Kochin::ComputeNbKochinAngles(){

    // Compute the number of Kochin angles.

    m_nb_kochin_angle = 0;
    for (double theta = 0; theta <= MU_2PI; theta += m_kochin_step) {
      ++m_nb_kochin_angle;
    }

  }

  int Kochin::GetNbKochinAngles() {

    // Getter of the number of Kochin angles.

    return m_nb_kochin_angle;

  }

  void Kochin::SetDiffractionKochin(unsigned int iwave, unsigned int iw, const Eigen::VectorXd &diffractionKochinVector) {

    // Setter of the diffraction Kochin function for a single wave frequency, all angles and a single wave direction.

    assert(iwave < m_HDB->GetWaveDirectionDiscretization().size());
    assert(iw < m_HDB->GetFrequencyDiscretization().size());
    assert(diffractionKochinVector.rows() == m_nb_kochin_angle);
    assert(diffractionKochinVector.cols() == 1);
    m_kochin_diffraction[iwave].col(iw) = diffractionKochinVector;

  }

  void Kochin::SetDiffractionKochin(unsigned int iwave, const Eigen::MatrixXd &diffractionKochinMatrix) {

    // Setter of the diffraction Kochin function for all wave frequencies, all angles and a single wave direction.

    assert(iwave < m_HDB->GetWaveDirectionDiscretization().size());
    assert(diffractionKochinMatrix.rows() == m_nb_kochin_angle);
    assert(diffractionKochinMatrix.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_kochin_diffraction[iwave] = diffractionKochinMatrix;

  }

  Eigen::MatrixXd Kochin::GetDiffractionKochin(unsigned int iwave) {

    // Getter of the diffraction Kochin function for all single wave frequencies, all angles and a single wave direction.

    return m_kochin_diffraction[iwave];

  }

  void Kochin::SetDiffractionKochinDerivative(unsigned int iwave, unsigned int iw
                                            , const Eigen::VectorXd &diffractionKochinDerivativeVector) {

    // Setter of the angular derivative of the diffraction Kochin function for a single wave frequency, all angles and a single wave direction.

    assert(iwave < m_HDB->GetWaveDirectionDiscretization().size());
    assert(iw < m_HDB->GetFrequencyDiscretization().size());
    assert(diffractionKochinDerivativeVector.rows() == m_nb_kochin_angle);
    assert(diffractionKochinDerivativeVector.cols() == 1);
    m_kochin_diffraction_derivate[iwave].col(iw) = diffractionKochinDerivativeVector;

  }

  void Kochin::SetDiffractionKochinDerivative(unsigned int iwave, const Eigen::MatrixXd &diffractionKochinDerivativeMatrix) {

    // Setter of the angular derivative of the diffraction Kochin function for all wave frequencies, all angles and a single wave direction.

    assert(iwave < m_HDB->GetWaveDirectionDiscretization().size());
    assert(diffractionKochinDerivativeMatrix.rows() == m_nb_kochin_angle);
    assert(diffractionKochinDerivativeMatrix.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_kochin_diffraction_derivate[iwave] = diffractionKochinDerivativeMatrix;

  }

  Eigen::MatrixXd Kochin::GetDiffractionKochinDerivative(unsigned int iwave) {

    // Getter of the angular derivative of the diffraction Kochin function for all single wave frequencies, all angles and a single wave direction.

    return m_kochin_diffraction_derivate[iwave];

  }

  void Kochin::SetRadiationKochin(Body *Body, const Eigen::VectorXd &radiationKochinVector) {

    // Setter of the radiation Kochin function for a body, a single wave frequency and all angles.

    // Initialization.
    if (m_kochin_radiation.count(Body) == 0) {
      std::vector<Eigen::VectorXd> tempVect;
      tempVect.reserve(m_HDB->GetFrequencyDiscretization().size());
      m_kochin_radiation.insert(std::make_pair(Body, tempVect));
    }

    assert(radiationKochinVector.rows() == m_nb_kochin_angle);
    m_kochin_radiation.at(Body).push_back(radiationKochinVector);

  }

  void Kochin::SetRadiationKochinDerivative(Body *Body, const Eigen::VectorXd &radiationKochinDerivativeVector) {

    // Setter of the radiation Kochin function for a body, a single wave frequency and all angles.

    // Initialization.
    if (m_kochin_radiation_derivate.count(Body) == 0) {
      std::vector<Eigen::VectorXd> tempVect;
      tempVect.reserve(m_HDB->GetFrequencyDiscretization().size());
      m_kochin_radiation_derivate.insert(std::make_pair(Body, tempVect));
    }

    assert(radiationKochinDerivativeVector.rows() == m_nb_kochin_angle);
    m_kochin_radiation_derivate.at(Body).push_back(radiationKochinDerivativeVector);

  }

}
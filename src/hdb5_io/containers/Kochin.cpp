//
// Created by pierre-yves on 13/01/2021.
//

#include "Kochin.h"
#include "HydrodynamicDataBase.h"

namespace HDB5_io {

  Kochin::Kochin(HydrodynamicDataBase *hdb) : m_HDB(hdb){

    // Constructor of the class.

    // Allocations.
    auto nDirections = m_HDB->GetWaveDirectionDiscretization().size();
    m_kochin_diffraction.reserve((unsigned long) nDirections);
    m_kochin_diffraction_derivate.reserve((unsigned long) nDirections);

  }

  void Kochin::SetKochinStep(const double &kochin_step) {
    m_kochin_step = kochin_step;
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

  void Kochin::SetDiffractionKochin(unsigned int iwave, unsigned int iw, const Eigen::MatrixXcd &diffractionKochinVector) {

    // Setter of the diffraction Kochin function for a single wave frequency, all angles and all wave directions.

    assert(iwave < m_HDB->GetWaveDirectionDiscretization().size());
    assert(iw < m_HDB->GetFrequencyDiscretization().size());
    assert(diffractionKochinVector.rows() == m_nb_kochin_angle);
    assert(diffractionKochinVector.cols() == 1);
    m_kochin_diffraction[iwave].col(iw) = diffractionKochinVector;

  }

  void Kochin::SetDiffractionKochinDerivate(unsigned int iwave, unsigned int iw
                                            , const Eigen::MatrixXcd &diffractionKochinDerivateVector) {

    // Setter of the angular derivative of the diffraction Kochin function for a single wave frequency, all angles and all wave directions.

    assert(iwave < m_HDB->GetWaveDirectionDiscretization().size());
    assert(iw < m_HDB->GetFrequencyDiscretization().size());
    assert(diffractionKochinDerivateVector.rows() == m_nb_kochin_angle);
    assert(diffractionKochinDerivateVector.cols() == 1);
    m_kochin_diffraction_derivate[iwave].col(iw) = diffractionKochinDerivateVector;

  }

}
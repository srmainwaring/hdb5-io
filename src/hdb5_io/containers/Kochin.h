//
// Created by pierre-yves on 13/01/2021.
//

#ifndef HDB5_IO_KOCHIN_H
#define HDB5_IO_KOCHIN_H

#include "MathUtils/MathUtils.h"
#include <vector>
#include <unordered_map>

#include "Body.h"

namespace HDB5_io {

  // Forward declaration.
  class HydrodynamicDataBase;

  /**
  * Class for managing the different Kochin functions and their angular derivative.
  */
  class Kochin {

   public:

    /// Constructor of the class.
    Kochin(HydrodynamicDataBase *hdb);

    /// Setter for the Kochin angular step.
    void SetKochinStep(const double &kochin_step);

    /// Getter for the Kochin angular step.
    double GetKochinStep() const;

    /// Compute the number of Kochin angles.
    void ComputeNbKochinAngles();

    /// Setter of the diffraction Kochin function for a single wave frequency, all angles and all wave directions.
    void SetDiffractionKochin(unsigned int iwave, unsigned int iw, const Eigen::MatrixXcd &diffractionKochinVector);

    /// Setter of the angular derivative of the diffraction Kochin function for a single wave frequency, all angles and all wave directions.
    void SetDiffractionKochinDerivate(unsigned int iwave, unsigned int iw, const Eigen::MatrixXcd &diffractionKochinDerivateVector);
    
   private:

    /// HDB containing this data container.
    HydrodynamicDataBase *m_HDB;

    /// Kochin angular step.
    double m_kochin_step;

    /// Number of kochin angles.
    int m_nb_kochin_angle;

    /// Diffraction elementary Kochin functions for all wave directions (std::vector),
    /// all angles (rows of Eigen::MatrixXcd) and all wave frequencies (column of Eigen::MatrixXcd).
    std::vector<Eigen::MatrixXcd> m_kochin_diffraction;

    /// Angular differentiation of the diffraction elementary Kochin functions for all wave directions (std::vector),
    /// all angles (rows of Eigen::MatrixXcd) and all wave frequencies (column of Eigen::MatrixXcd).
    std::vector<Eigen::MatrixXcd> m_kochin_diffraction_derivate;

    /// Radiation elementary Kochin functions for all bodies (Body), all angles (rows of Eigen::VectorXcd)
    /// and all wave frequencies (std::vector).
    std::unordered_map<Body *, std::vector<Eigen::VectorXcd>> m_kochin_radiation;

    /// Angular differentiation of the radiation elementary Kochin functions for all bodies (Body), all angles (rows of Eigen::VectorXcd)
    //    /// and all wave frequencies (std::vector).
    std::unordered_map<Body *, std::vector<Eigen::VectorXcd>> m_kochin_radiation_derivate;

  };

}

#endif //HDB5_IO_KOCHIN_H

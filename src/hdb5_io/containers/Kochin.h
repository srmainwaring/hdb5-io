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
    Kochin(HydrodynamicDataBase *hdb, const double &kochin_step);

    /// Getter for the Kochin angular step.
    double GetKochinStep() const;

    /// Compute the number of Kochin angles.
    void ComputeNbKochinAngles();

    /// Getter of the number of Kochin angles.
    int GetNbKochinAngles();

    /// Setter of the diffraction Kochin function for a single wave frequency, all angles and a single wave direction.
    void SetDiffractionKochin(unsigned int iwave, unsigned int iw, const Eigen::VectorXcd &diffractionKochinVector);

    /// Setter of the diffraction Kochin function for all wave frequencies, all angles and a single wave direction.
    void SetDiffractionKochin(unsigned int iwave, const Eigen::MatrixXcd &diffractionKochinMatrix);

    /// Getter of the diffraction Kochin function for all wave frequencies, all angles and a single wave direction.
    Eigen::MatrixXcd GetDiffractionKochin(unsigned int iwave);

    /// Setter of the angular derivative of the diffraction Kochin function for a single wave frequency, all angles and a single wave direction.
    void SetDiffractionKochinDerivative(unsigned int iwave, unsigned int iw, const Eigen::VectorXcd &diffractionKochinDerivativeVector);

    /// Setter of the angular derivative of the diffraction Kochin function for all wave frequencies, all angles and single wave direction.
    void SetDiffractionKochinDerivative(unsigned int iwave, const Eigen::MatrixXcd &diffractionKochinDerivativeMatrix);

    /// Getter of the angular derivative of the diffraction Kochin function for all wave frequencies, all angles and a single wave direction.
    Eigen::MatrixXcd GetDiffractionKochinDerivative(unsigned int iwave);

    /// Setter of the radiation Kochin function for a single body, a single dof, all wave frequencies and all angles.
    void SetRadiationKochin(Body *Body, const Eigen::MatrixXcd &radiationKochinMatrix);

    /// Setter of the radiation Kochin function for a single body, a single dof, a single wave frequency and all angles.
    void SetRadiationKochin(Body *Body, unsigned int iw, unsigned int idof, const Eigen::MatrixXcd &radiationKochinMatrix);

    /// Getter of the radiation Kochin function for a single body, all wave frequencies, a single dof and all angles.
    Eigen::MatrixXcd GetRadiationKochin(Body *Body, unsigned int idof);

    /// Setter of the angular derivative of the radiation Kochin function for a single body, a single wave frequency, all dof and all angles.
    void SetRadiationKochinDerivative(Body *Body, const Eigen::MatrixXcd &radiationKochinDerivativeMatrix);

    /// Setter of the angular derivative of the radiation Kochin function for a single body, a single wave frequency, a single dof and all angles.
    void SetRadiationKochinDerivative(Body *Body, unsigned int iw, unsigned int idof, const Eigen::MatrixXcd &radiationKochinDerivativeMatrix);

    /// Getter of the angular derivative of the radiation Kochin function for a single body, all wave frequencies, a single dof and all angles.
    Eigen::MatrixXcd GetRadiationKochinDerivative(Body *Body, unsigned int idof);

   private:

    /// HDB containing this data container.
    HydrodynamicDataBase *m_HDB;

    /// Kochin angular step.
    double m_kochin_step;

    /// Number of kochin angles.
    int m_nb_kochin_angle;

    /// Diffraction elementary Kochin functions for all wave directions (std::vector),
    /// all angles (rows of Eigen::MatrixXd) and all wave frequencies (column of Eigen::MatrixXd).
    std::vector<Eigen::MatrixXcd> m_kochin_diffraction;

    /// Angular differentiation of the diffraction elementary Kochin functions for all wave directions (std::vector),
    /// all angles (rows of Eigen::MatrixXd) and all wave frequencies (column of Eigen::MatrixXd).
    std::vector<Eigen::MatrixXcd> m_kochin_diffraction_derivate;

    /// Radiation elementary Kochin functions for all bodies (Body), all dof (std::vector), all angles (rows of Eigen::MatrixXd)
    /// and all all wave frequencies (cols of Eigen::MatrixXd).
    std::unordered_map<Body *, std::vector<Eigen::MatrixXcd>> m_kochin_radiation;

    /// Angular differentiation of the radiation elementary Kochin functions for all bodies (Body), all dof (std::vector), all angles (rows of Eigen::MatrixXd)
    //    /// and all all wave frequencies (cols of Eigen::MatrixXd).
    std::unordered_map<Body *, std::vector<Eigen::MatrixXcd>> m_kochin_radiation_derivate;

  };

}

#endif //HDB5_IO_KOCHIN_H

//
// Created by lletourn on 16/04/20.
//

#ifndef HDB5_IO_POLERESIDUE_H
#define HDB5_IO_POLERESIDUE_H

#include <complex>
#include <vector>

#include <Eigen/Dense>

namespace HDB5_io {

  template<typename T>
  class PoleResiduePair {

   public:
    PoleResiduePair(const T &pole, const T &residue) : m_pole(pole), m_residue(residue) {}

    inline const T &pole() const {
      return m_pole;
    }

    inline const T &residue() const {
      return m_residue;
    }

   private:
    T m_pole;
    T m_residue;
  };

  using RealPoleResiduePair = PoleResiduePair<double>;
  using CCPoleResiduePair = PoleResiduePair<std::complex<double>>;

  class PoleResidue {

   public:

    PoleResidue() = default;

    unsigned int nb_real_poles() const;

    unsigned int nb_cc_poles() const;

    unsigned int order() const;

    void AddPoleResidue(const double &pole, const double &residue);

    void AddPoleResidue(const std::complex<double> &pole, const std::complex<double> &residue);

    std::vector<RealPoleResiduePair> GetRealPairs() const;

    std::vector<CCPoleResiduePair> GetComplexPairs() const;

    std::vector<double> GetRealPoles() const;

    std::vector<double> GetRealResidues() const;

    Eigen::MatrixXd GetComplexPoles() const;

    Eigen::MatrixXd GetComplexResidues() const;

   private:

    std::vector<RealPoleResiduePair> m_real_pairs;
    std::vector<CCPoleResiduePair> m_cc_pairs;

  };

} // end namespace HDB5_io
#endif //HDB5_IO_POLERESIDUE_H

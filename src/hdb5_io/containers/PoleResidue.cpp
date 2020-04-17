//
// Created by lletourn on 16/04/20.
//

#include "PoleResidue.h"

namespace HDB5_io {


  std::vector<double> PoleResidue::GetRealPoles() const {
    std::vector<double> pole;
    for (auto pair : m_real_pairs) {
      pole.emplace_back(pair.pole());
    }
    return pole;
  }

  std::vector<double> PoleResidue::GetRealResidues() const {
    std::vector<double> residue;
    for (auto pair : m_real_pairs) {
      residue.emplace_back(pair.residue());
    }
    return residue;
  }

  Eigen::MatrixXd PoleResidue::GetComplexPoles() const {
    Eigen::MatrixXd matrix(2, nb_cc_poles());
    int i=0;
    for (auto pair : m_cc_pairs) {
      auto pole = pair.pole();
      matrix.col(i) << pole.real(), pole.imag();
      i++;
    }
    return matrix;
  }

  Eigen::MatrixXd PoleResidue::GetComplexResidues() const {
    Eigen::MatrixXd matrix(2, nb_cc_poles());
    int i=0;
    for (auto pair : m_cc_pairs) {
      auto residue = pair.residue();
      matrix.col(i) << residue.real(), residue.imag();
      i++;
    }
    return matrix;
  }

  void PoleResidue::AddPoleResidue(const double &pole, const double &residue) {
    m_real_pairs.emplace_back(pole, residue);
  }

  void PoleResidue::AddPoleResidue(const std::complex<double> &pole, const std::complex<double> &residue) {
    m_cc_pairs.emplace_back(pole, residue);
  }

  std::vector<RealPoleResiduePair> PoleResidue::GetRealPairs() const {
    return m_real_pairs;
  }

  std::vector<CCPoleResiduePair> PoleResidue::GetComplexPairs() const {
    return m_cc_pairs;
  }

  unsigned int PoleResidue::nb_real_poles() const {return m_real_pairs.size();}

  unsigned int PoleResidue::nb_cc_poles() const {return m_cc_pairs.size();}

  unsigned int PoleResidue::order() const {return nb_real_poles() + 2 * nb_cc_poles();}
}
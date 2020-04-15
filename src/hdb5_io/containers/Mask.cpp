//
// Created by lletourn on 27/02/20.
//

#include "Mask.h"

namespace HDB5_io {

  void Mask::SetMask(mathutils::Vector6d<int> mask) {

    for (unsigned int i = 0; i < 6; i++) { assert(mask(i) == 0 or mask(i) == 1); }

    for (unsigned int i = 0; i < 6; i++) {
      if (mask(i) == 1) {
        m_mask(i) = true;
        m_listDOF.push_back(i);
      } else {
        m_mask(i) = false;
      }
    }

    m_nbDOF = (unsigned int) mask.sum();

    m_matrix = Eigen::MatrixXd::Zero(6, m_nbDOF);
    unsigned int j = 0;
    for (unsigned int i = 0; i < 6; i++) {
      if (mask(i)) {
        m_matrix(i, j) = 1.;
        j += 1;
      }
    }
  }

  mathutils::Vector6d<bool> Mask::GetMask() const {
    return m_mask;
  }

  mathutils::MatrixMN<double> Mask::GetMatrix() const {
    return m_matrix;
  }


} // namespace HDB5_io

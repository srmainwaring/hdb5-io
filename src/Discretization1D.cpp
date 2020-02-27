//
// Created by lletourn on 26/02/20.
//

#include "Discretization1D.h"
#include "MathUtils/VectorGeneration.h"

namespace HDB5_io {

  // ----------------------------------------------
  // Discretization1D
  // ----------------------------------------------

  std::vector<double> Discretization1D::GetVector() const {
    return mathutils::linspace<double>(m_xmin, m_xmax, m_nx);
  }

  void Discretization1D::SetStep(double delta) {
    m_nx = 1 + (unsigned int) ((m_xmax - m_xmin) / delta);
  }

  double Discretization1D::GetStep() const {
    return (m_xmax - m_xmin) / double(m_nx - 1);
  }


} // namespace HDB5_io

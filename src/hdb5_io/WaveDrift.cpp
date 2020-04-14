//
// Created by lletourn on 07/04/20.
//

#include "WaveDrift.h"

namespace HDB5_io {




  WaveDrift::WaveDrift(): m_symmetry_X(false), m_symmetry_Y(false) {
    m_data = std::make_unique<mathutils::LookupTable2d<double>>();
  }

  void WaveDrift::SetSymmetries(bool symmetry_X, bool symmetry_Y) {
    m_symmetry_X = symmetry_X;
    m_symmetry_Y = symmetry_Y;
  }

  std::vector<bool> WaveDrift::GetSymmetries() const {
    return {m_symmetry_X, m_symmetry_Y};
  }

  void WaveDrift::SetFrequencies(const std::vector<double> &frequencies) {
    m_data->SetX(frequencies);
  }

  void WaveDrift::SetWaveDirections(const std::vector<double> &angles) {
    m_data->SetY(angles);
  }

  void WaveDrift::AddData(const std::string &name, const std::vector<double> &coeffs) {
    m_data->AddData(name, coeffs);
  }

  double WaveDrift::Eval(const std::string &name, double frequency, double angle) const {
    return m_data->Eval(name, frequency, angle);
  }
}
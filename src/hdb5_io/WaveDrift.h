//
// Created by lletourn on 07/04/20.
//

#ifndef HDB5_IO_WAVEDRIFT_H
#define HDB5_IO_WAVEDRIFT_H

#include <memory>
#include "Discretization1D.h"
#include "MathUtils/LookupTable2D.h"
#include "MathUtils/Matrix.h"

namespace HDB5_io {

  class WaveDrift {

   public:

    WaveDrift();

    void SetSymmetries(bool symmetry_X, bool symmetry_Y);

    std::vector<bool> GetSymmetries() const;

    void SetFrequencies(const std::vector<double>& frequencies);

    void SetWaveDirections(const std::vector<double>& angles);

    /// Adding new data with wave drift coefficients
    /// \param name Name of the data
    /// \param coeffs >Wave drift coefficient values
    void AddData(const std::string &name, const std::vector<double>& coeffs);

    double Eval(const std::string &name, double frequency, double angle) const;

//    double Eval(const std::string &name, std::vector<double> frequencies, std::vector<double> angles) const {
//      return m_data->Eval(name, frequencies, angles);
//    }

   private:

    bool m_symmetry_X;
    bool m_symmetry_Y;

    std::unique_ptr<mathutils::LookupTable2d<double>> m_data;

  };

}

#endif //HDB5_IO_WAVEDRIFT_H
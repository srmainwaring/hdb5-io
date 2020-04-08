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

    void SetSymmetries(bool symmetry_X, bool symmetry_Y) {
      m_symmetry_X = symmetry_X;
      m_symmetry_Y = symmetry_Y;
    }

//    void SetSurge(std::vector<double> frequencies, std::vector<double> waveDirections, mathutils::MatrixMN<double> data);

   private:

    bool m_symmetry_X;
    bool m_symmetry_Y;

    std::unique_ptr<mathutils::LookupTable2d<double>> m_data;

  };

}

#endif //HDB5_IO_WAVEDRIFT_H

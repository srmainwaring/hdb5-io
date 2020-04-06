//
// Created by lletourn on 06/04/20.
//

#ifndef HDB5_IO_MESH_H
#define HDB5_IO_MESH_H

#include "meshoui/mesh.h"
#include "MathUtils/Matrix.h"

namespace HDB5_io {

  class Mesh : public meshoui::Mesh {

   public:

    Eigen::MatrixXd GetVertices() const;

    Eigen::MatrixXi GetFaces();

  };

}

#endif //HDB5_IO_MESH_H

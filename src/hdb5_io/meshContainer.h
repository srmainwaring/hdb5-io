//
// Created by lletourn on 03/04/20.
//

#ifndef SANDBOX_VTK_MESHCONTAINER_H
#define SANDBOX_VTK_MESHCONTAINER_H

#include <memory>
#include "meshoui/mesh.h"
#include "meshoui/vtkmesh.h"

class meshContainer {
 public:

  explicit meshContainer(const std::string &filename) {
    m_mesh = std::make_shared<meshoui::Mesh>();
    m_mesh->Load(filename);
  }


  meshContainer() {
    m_mesh = std::make_shared<meshoui::Mesh>();
    // Build a box
    std::vector<meshoui::Vector3d> vertices;
    vertices.emplace_back(-1, -1,  1);
    vertices.emplace_back(1, -1,  1);
    vertices.emplace_back(1,  1,  1);
    vertices.emplace_back(-1,  1,  1);
    vertices.emplace_back(-1, -1, -1);
    vertices.emplace_back(1, -1, -1);
    vertices.emplace_back(1,  1, -1);
    vertices.emplace_back(-1,  1, -1);
    std::vector<Eigen::VectorXi> faces;
    faces.emplace_back(Eigen::Vector3i(0,1,2));
    faces.emplace_back(Eigen::Vector3i(2,3,0));
    faces.emplace_back(Eigen::Vector3i(0,4,1));
    faces.emplace_back(Eigen::Vector3i(1,4,5));
    faces.emplace_back(Eigen::Vector3i(1,5,2));
    faces.emplace_back(Eigen::Vector3i(2,5,6));
    faces.emplace_back(Eigen::Vector3i(2,6,3));
    faces.emplace_back(Eigen::Vector3i(3,6,7));
    faces.emplace_back(Eigen::Vector3i(3,7,0));
    faces.emplace_back(Eigen::Vector3i(0,7,4));
    faces.emplace_back(Eigen::Vector3i(6,5,4));
    faces.emplace_back(Eigen::Vector3i(7,6,4));
    // Load it
    m_mesh->Load(vertices, faces);
  }

  void Visualize() const {
    // Creation of a VTK mesh.
    meshoui::VTKMesh vtkmesh = meshoui::VTKMesh(*m_mesh);

    // Visualization.
    vtkmesh.Visualize();
  }

 private:
  std::shared_ptr<meshoui::Mesh> m_mesh;

};


#endif //SANDBOX_VTK_MESHCONTAINER_H

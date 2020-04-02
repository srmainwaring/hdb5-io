//
// Created by lletourn on 27/02/20.
//

#include "HDB5_io.h"
#include "meshoui/meshoui.h"

using namespace meshoui;

using namespace HDB5_io;

int main() {

  HydrodynamicDataBase HDB;

  HDB.Import_HDF5("/home/lletourn/Documents/DEV/hdb5-io/test.hdb5");

  auto mesh = HDB.GetBody(0)->GetMesh();

  // VTKMesh.
  VTKMesh vtkmesh = VTKMesh(*mesh);

  // Visualization.
  vtkmesh.Visualize();


  HDB.Export_HDF5("out.hdb5");

  return 0;
}
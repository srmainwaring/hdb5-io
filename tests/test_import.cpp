//
// Created by lletourn on 27/02/20.
//

#include "hdb5_io/HDB5_io.h"

using namespace HDB5_io;

int main() {

//  auto HDB = std::make_shared<HydrodynamicDataBase>();
//
//  HDB->Import_HDF5("/home/lletourn/Documents/DEV/hdb5-io/Sphere_1484_faces_Kochin.hdb5");

  auto HDB = import_HDB("/home/lletourn/Documents/DEV/hdb5-io/Sphere_1484_faces_Kochin.hdb5");

//#ifdef H5_USE_VTK
//  HDB.GetBody(0)->VisualizeMesh();
//#endif

  HDB->Export_HDF5("out.hdb5");

  return 0;
}
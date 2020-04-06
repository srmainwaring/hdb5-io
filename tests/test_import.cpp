//
// Created by lletourn on 27/02/20.
//

#include "hdb5_io/HDB5_io.h"
//#include "meshoui/meshoui.h"

//using namespace meshoui;

using namespace HDB5_io;

int main() {

  HydrodynamicDataBase HDB;

  HDB.Import_HDF5("/home/lletourn/Documents/DEV/hdb5-io/test.hdb5");

#ifdef H5_USE_VTK
  HDB.GetBody(0)->VisualizeMesh();
#endif

  HDB.Export_HDF5("out.hdb5");

  return 0;
}
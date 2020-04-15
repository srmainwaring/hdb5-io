//
// Created by lletourn on 27/02/20.
//

#include "hdb5_io/HDB5_io.h"

using namespace HDB5_io;

int main() {

  auto HDB = import_HDB("/home/lletourn/Documents/DEV/hdb5-io/Sphere_1484_faces_Kochin.hdb5");

#ifdef H5_USE_VTK
  HDB->GetBody(0)->VisualizeMesh();
#endif

  export_HDB("out.hdb5", HDB.get());

  return 0;
}
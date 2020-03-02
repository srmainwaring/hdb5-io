//
// Created by lletourn on 27/02/20.
//

#include "HDB5_io.h"

using namespace HDB5_io;

int main() {

  HydrodynamicDataBase HDB;

  HDB.Import_HDF5("/home/lletourn/Documents/DEV/hdb5-io/test.hdb5");

  return 0;
}
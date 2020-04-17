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

  Eigen::MatrixXi ordre(6,6);
  ordre <<  4, 5, 6, 6, 5, 6,
            5, 4, 6, 5, 5, 4,
            4, 5, 6, 6, 5, 6,
            5, 4, 6, 5, 5, 4,
            4, 5, 6, 6, 5, 6,
            5, 4, 6, 5, 5, 4;

  for (int idof = 0; idof<6; idof++) {

    std::vector<PoleResidue> modalCoeff;
    for (int iforce = 0; iforce<6; iforce ++) {
      int n_pole = ordre(idof, iforce);
      Eigen::VectorXd poles(n_pole), residues(n_pole);
      poles.setRandom();
      residues.setRandom();

      Eigen::VectorXcd cplxPoles(n_pole), cplxResidues(n_pole);
      cplxPoles.setRandom();
      cplxResidues.setRandom();

      PoleResidue pair;
      for (int i=0; i < n_pole; i++) {
//        std::cout<<poles(i)<<"; "<<residues(i)<<std::endl;
        pair.AddPoleResidue(poles(i), residues(i));
        pair.AddPoleResidue(cplxPoles(i), cplxResidues(i));
      }
      modalCoeff.emplace_back(pair);

    }
    HDB->GetBody(0)->SetModalCoefficients(HDB->GetBody(0), modalCoeff);
  }

  export_HDB("out.hdb5", HDB.get());

  auto HDB_2 = import_HDB("out.hdb5");

  export_HDB("out_2.hdb5", HDB_2.get());

  return 0;
}
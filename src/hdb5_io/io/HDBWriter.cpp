//
// Created by lletourn on 15/04/20.
//

#include <Eigen/Dense>

#include "HDBWriter.h"
#include <highfive/H5Easy.hpp>

#include "../HydrodynamicDataBase.h"

#include "MathUtils/VectorN.h"
#include "MathUtils/Vector3d.h"

namespace HDB5_io {


  void HDBWriter::Write(const std::string &filename) const {

    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    WriteHDBBasics(file);

    if (m_hdb->GetWaveDrift()) {
      WriteWaveDrift(file);
    }

    auto bodies = file.createGroup("Bodies");

    for (unsigned int i = 0; i < m_hdb->GetNbBodies(); i++) {
      auto bodyGroup = bodies.createGroup("Body_" + std::to_string(i));

      auto body = m_hdb->GetBody(i);

      auto bodyPath = "Bodies/Body_" + std::to_string(i);

      WriteBodyBasics(file, "Bodies/Body_" + std::to_string(i), body);

      WriteMesh(file, bodyPath, body);

      WriteExcitation(excitationType::Diffraction, file, bodyPath + "/Excitation/Diffraction", body);
      WriteExcitation(excitationType::Froude_Krylov, file, bodyPath + "/Excitation/FroudeKrylov", body);

      WriteRadiation(file, bodyPath + "/Radiation", body);

      if (body->HasRAO())
        WriteRAO(file, bodyPath + "/RAO", body);
    }

  }


  void HDBWriter::WriteHDBBasics(HighFive::File &file) const {

    auto creationDate = m_hdb->GetCreationDate();
    HighFive::DataSet dataSet = file.createDataSet<std::string>("CreationDate", HighFive::DataSpace::From(creationDate));
    dataSet.write(creationDate);
    dataSet.createAttribute<std::string>("Description", "Date of the creation of this database.");

    auto version = 3.;
    dataSet = file.createDataSet<double>("Version", HighFive::DataSpace::From(version));
    dataSet.write(version);
    dataSet.createAttribute<std::string>("Description", "Version of the hdb5 output file.");

    auto solver = m_hdb->GetSolver();
    dataSet = file.createDataSet<std::string>("Solver", HighFive::DataSpace::From(solver));
    dataSet.write(solver);
    dataSet.createAttribute<std::string>("Description",
                                         "Hydrodynamic solver used for computing the hydrodynamic database.");

    auto nbody = m_hdb->GetNbBodies();
    dataSet = file.createDataSet<int>("NbBody", HighFive::DataSpace::From(nbody));
    dataSet.write(nbody);
    dataSet.createAttribute<std::string>("Description", "Number of hydrodynamic bodies.");

    auto normalizationLength = m_hdb->GetNormalizationLength();
    dataSet = file.createDataSet<double>("NormalizationLength", HighFive::DataSpace::From(normalizationLength));
    dataSet.write(normalizationLength);
    dataSet.createAttribute<std::string>("Description", "Normalization length.");
    dataSet.createAttribute<std::string>("Unit", "m");

    auto gravityAcceleration = m_hdb->GetGravityAcceleration();
    dataSet = file.createDataSet<double>("GravityAcc", HighFive::DataSpace::From(gravityAcceleration));
    dataSet.write(gravityAcceleration);
    dataSet.createAttribute<std::string>("Description", "Gravity acceleration.");
    dataSet.createAttribute<std::string>("Unit", "m/s**2");

    auto waterDensity = m_hdb->GetWaterDensity();
    dataSet = file.createDataSet<double>("WaterDensity", HighFive::DataSpace::From(waterDensity));
    dataSet.write(waterDensity);
    dataSet.createAttribute<std::string>("Description", "Water Density.");
    dataSet.createAttribute<std::string>("Unit", "kg/m**3");

    auto waterDepth = m_hdb->GetWaterDepth();
    dataSet = file.createDataSet<double>("WaterDepth", HighFive::DataSpace::From(waterDepth));
    dataSet.write(waterDepth);
    dataSet.createAttribute<std::string>("Description",
                                         "Water depth: 0 for infinite depth and positive for finite depth.");
    dataSet.createAttribute<std::string>("Unit", "m");
  }

  void HDBWriter::WriteDiscretizations(HighFive::File &file) const {

    auto discretizations = file.createGroup("Discretizations");
    H5Easy::dump(file, "Discretizations/Frequency",
                 static_cast<Eigen::Matrix<double, Eigen::Dynamic, 1>> (m_hdb->GetFrequencyDiscretization()));
    H5Easy::dump(file, "Discretizations/WaveDirection",
                 static_cast<Eigen::Matrix<double, Eigen::Dynamic, 1>> (m_hdb->GetWaveDirectionDiscretization()));
    H5Easy::dump(file, "Discretizations/Time",
                 static_cast<Eigen::Matrix<double, Eigen::Dynamic, 1>> (m_hdb->GetTimeDiscretization()));

  }

  void HDBWriter::WriteBodyBasics(HighFive::File &file, const std::string &path, Body *body) const {

    auto bodyGroup = file.getGroup(path);

    auto dataSet = bodyGroup.createDataSet<unsigned int>("ID", HighFive::DataSpace::From(body->GetID()));
    dataSet.write(body->GetID());
    dataSet.createAttribute<std::string>("Description", "Body index");

    dataSet = bodyGroup.createDataSet<std::string>("BodyName", HighFive::DataSpace::From(body->GetName()));
    dataSet.write(body->GetName());
    dataSet.createAttribute<std::string>("Description", "Body name");

    H5Easy::dump(file, path + "/BodyPosition", static_cast<Eigen::Vector3d> (body->GetPosition()));
    bodyGroup.getDataSet("BodyPosition").createAttribute<std::string>("Description",
                                                                      "Center of gravity of the body in the absolute frame");

    bodyGroup.createGroup("Mask");
//    H5Easy::dump(file, path + "/Mask/ForceMask", static_cast<Eigen::Matrix<bool, 6, 1>> (body->GetForceMask().GetMask()));
//    H5Easy::dump(file, path + "/Mask/MotionMask", static_cast<Eigen::Matrix<bool, 6, 1>> (body->GetMotionMask().GetMask()));

  }

  void HDBWriter::WriteMesh(HighFive::File &file, const std::string &path, Body *body) const {

    auto bodyGroup = file.getGroup(path);

    bodyGroup.createGroup("Mesh");
    H5Easy::dump(file, path + "/Mesh/NbFaces", body->GetMesh()->n_faces());
    H5Easy::dump(file, path + "/Mesh/Faces", body->GetMesh()->GetFaces());
    H5Easy::dump(file, path + "/Mesh/NbVertices", body->GetMesh()->n_vertices());
    H5Easy::dump(file, path + "/Mesh/Vertices", body->GetMesh()->GetVertices());

  }

  void HDBWriter::WriteExcitation(HDBWriter::excitationType type, HighFive::File &HDF5_file, const std::string &path,
                                  Body *body) const {

    auto angles = m_hdb->GetWaveDirectionDiscretization();

    for (unsigned int iwaveDir = 0; iwaveDir < angles.size(); ++iwaveDir) {

      auto anglePath = path + "/Angle_" + std::to_string(iwaveDir);

      auto angle = angles(iwaveDir);

      Eigen::MatrixXcd coeff;
      switch (type) {
        case Diffraction : {
          coeff = body->GetDiffraction(iwaveDir);
          break;
        }
        case Froude_Krylov : {
          coeff = body->GetFroudeKrylov(iwaveDir);
          break;
        }
      }
      auto angleGroup = HDF5_file.createGroup(anglePath);
      H5Easy::dump(HDF5_file, anglePath + "/Angle", angle);
      angleGroup.getDataSet("Angle").createAttribute<std::string>("Description", "Wave direction.");
      angleGroup.getDataSet("Angle").createAttribute<std::string>("Unit", "deg");

      H5Easy::dump(HDF5_file, anglePath + "/RealCoeffs", static_cast<Eigen::MatrixXd>(coeff.real()));
      angleGroup.getDataSet("RealCoeffs").createAttribute<std::string>("Description",
                                                                       "Real part of the Froude-Krylov loads on body " +
                                                                       std::to_string(body->GetID()) +
                                                                       " for a wave direction of " + std::to_string(
                                                                           angle) +
                                                                       " deg.");
      angleGroup.getDataSet("RealCoeffs").createAttribute<std::string>("Unit", "N/m");

      H5Easy::dump(HDF5_file, anglePath + "/ImagCoeffs", static_cast<Eigen::MatrixXd>(coeff.imag()));
      angleGroup.getDataSet("ImagCoeffs").createAttribute<std::string>("Description",
                                                                       "Imaginary part of the Froude-Krylov loads on body " +
                                                                       std::to_string(body->GetID()) +
                                                                       " for a wave direction of " + std::to_string(
                                                                           angle) +
                                                                       " deg.");
      angleGroup.getDataSet("ImagCoeffs").createAttribute<std::string>("Unit", "N/m");

    }

  }

  void HDBWriter::WriteRadiation(HighFive::File &HDF5_file, const std::string &path, Body *body) const {

    auto frequencies = m_hdb->GetFrequencyDiscretization();
    auto time = m_hdb->GetTimeDiscretization();

    for (unsigned int ibodyMotion = 0; ibodyMotion < m_hdb->GetNbBodies(); ++ibodyMotion) {

      auto bodyMotion = m_hdb->GetBody(ibodyMotion);
      auto bodyMotionPath = path + "/BodyMotion_" + std::to_string(ibodyMotion);
      auto bodyMotionGroup = HDF5_file.createGroup(bodyMotionPath);
      bodyMotionGroup.createAttribute("Description", "Hydrodynamic coefficients for motion of body " +
                                                     std::to_string(bodyMotion->GetID()) +
                                                     " that radiates waves and  generate force on body " +
                                                     std::to_string(body->GetID()) + ".");

      // Writing the infinite added mass matrix for the body.
      H5Easy::dump(HDF5_file, bodyMotionPath + "/InfiniteAddedMass",
                   static_cast<Eigen::MatrixXd>(body->GetInfiniteAddedMass(bodyMotion)));
      auto InfiniteAddedMass = bodyMotionGroup.getDataSet("InfiniteAddedMass");
      InfiniteAddedMass.createAttribute("Description",
                                        "Infinite added mass matrix that modifies the apparent mass of body " +
                                        std::to_string(bodyMotion->GetID()) +
                                        " from acceleration of body  " +
                                        std::to_string(body->GetID()) + ".");

      // Writing the radiation mask matrix for the body.
      H5Easy::dump(HDF5_file, bodyMotionPath + "/RadiationMask",
                   static_cast<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>>(body->GetRadiationMask(
                       bodyMotion)));
      auto RadiationMask = bodyMotionGroup.getDataSet("RadiationMask");
      RadiationMask.createAttribute("Description", "Radiation mask of body " +
                                                   std::to_string(bodyMotion->GetID()) +
                                                   " from acceleration of body  " +
                                                   std::to_string(body->GetID()) + ".");

      // Writing the added mass matrix for the body.
      auto AddedMassGroup = bodyMotionGroup.createGroup("AddedMass");
      AddedMassGroup.createAttribute("Description", "Added mass coefficients for acceleration of body " +
                                                    std::to_string(bodyMotion->GetID()) +
                                                    " that radiates waves and  generate force on body " +
                                                    std::to_string(body->GetID()) + ".");
      for (unsigned int i = 0; i < 6; i++) {
        H5Easy::dump(HDF5_file, bodyMotionPath + "/AddedMass/DOF_" + std::to_string(i),
                     body->GetHDBInterpolatedData(Body::interpolatedData::ADDED_MASS, bodyMotion, i, frequencies));
        auto DOF = AddedMassGroup.getDataSet("DOF_" + std::to_string(i));
        DOF.createAttribute("Description", "Added mass coefficients for an acceleration of body " +
                                           std::to_string(bodyMotion->GetID()) +
                                           " and force on body " +
                                           std::to_string(body->GetID()) + ".");
        //TODO : complete units, depending on DOF, bodyMotion, etc.
        DOF.createAttribute("Unit", "");
      }

      // Writing the radiation damping matrix for the body.
      auto RadiationDampingGroup = bodyMotionGroup.createGroup("RadiationDamping");
      RadiationDampingGroup.createAttribute("Description", "Damping coefficients for velocity of body " +
                                                           std::to_string(bodyMotion->GetID()) +
                                                           " that radiates waves and generates forces on body " +
                                                           std::to_string(body->GetID()) + ".");
      for (unsigned int i = 0; i < 6; i++) {
        H5Easy::dump(HDF5_file, bodyMotionPath + "/RadiationDamping/DOF_" + std::to_string(i),
                     body->GetHDBInterpolatedData(Body::interpolatedData::RADIATION_DAMPING, bodyMotion, i,
                                                  frequencies));
        auto DOF = RadiationDampingGroup.getDataSet("DOF_" + std::to_string(i));
        DOF.createAttribute("Description", "Wave damping coefficients for an acceleration of body " +
                                           std::to_string(bodyMotion->GetID()) +
                                           " and force on body " +
                                           std::to_string(body->GetID()) + ".");
        //TODO : complete units, depending on DOF, bodyMotion, etc.
        DOF.createAttribute("Unit", "");
      }

      // Writing the impulse response function K for the body.
      auto KGroup = bodyMotionGroup.createGroup("ImpulseResponseFunctionK");
      KGroup.createAttribute("Description", "Impulse response functions K due to the velocity of body " +
                                            std::to_string(bodyMotion->GetID()) +
                                            " that radiates waves and generates forces on body " +
                                            std::to_string(body->GetID()) + ".");
      for (unsigned int i = 0; i < 6; i++) {
        H5Easy::dump(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionK/DOF_" + std::to_string(i),
                     body->GetHDBInterpolatedData(Body::interpolatedData::IRF_K, bodyMotion, i, time));
        auto DOF = KGroup.getDataSet("DOF_" + std::to_string(i));
        DOF.createAttribute("Description", "Impulse response functions K");
      }

      // Writing the impulse response function KU for the body.
      auto KUGroup = bodyMotionGroup.createGroup("ImpulseResponseFunctionKU");
      KUGroup.createAttribute("Description", "Impulse response functions KU due to the velocity of body " +
                                             std::to_string(bodyMotion->GetID()) +
                                             " that radiates waves and generates forces on body " +
                                             std::to_string(body->GetID()) + ".");
      for (unsigned int i = 0; i < 6; i++) {
        H5Easy::dump(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionKU/DOF_" + std::to_string(i),
                     body->GetHDBInterpolatedData(Body::interpolatedData::IRF_KU, bodyMotion, i, time));
        auto DOF = KUGroup.getDataSet("DOF_" + std::to_string(i));
        DOF.createAttribute("Description", "Impulse response functions KU");
      }


    }

  }

  void HDBWriter::WriteRAO(HighFive::File &HDF5_file, const std::string &path, Body *body) const {

    auto angles = m_hdb->GetWaveDirectionDiscretization();

    for (unsigned int iwaveDir = 0; iwaveDir < angles.size(); ++iwaveDir) {

      auto anglePath = path + "/Angle_" + std::to_string(iwaveDir);

      auto angle = angles(iwaveDir);

      Eigen::MatrixXcd coeff = body->GetRAO(iwaveDir);

      auto angleGroup = HDF5_file.createGroup(anglePath);
      H5Easy::dump(HDF5_file, anglePath + "/Angle", angle);
      angleGroup.getDataSet("Angle").createAttribute<std::string>("Description", "Wave direction.");
      angleGroup.getDataSet("Angle").createAttribute<std::string>("Unit", "deg");

      H5Easy::dump(HDF5_file, anglePath + "/Amplitude", static_cast<Eigen::MatrixXd>(coeff.array().abs()));
      angleGroup.getDataSet("Amplitude").createAttribute<std::string>("Description",
                                                                      "Amplitude of the RAO of body " +
                                                                      std::to_string(body->GetID()) +
                                                                      " for a wave direction of " + std::to_string(
                                                                          angle) +
                                                                      " deg.");
      angleGroup.getDataSet("Amplitude").createAttribute<std::string>("Unit", "");

      std::cout<<"phase : "<<coeff(0,0)<<std::endl;
      H5Easy::dump(HDF5_file, anglePath + "/Phase", static_cast<Eigen::MatrixXd>(coeff.array().arg()));
      angleGroup.getDataSet("Phase").createAttribute<std::string>("Description",
                                                                  "Phase of the RAO of body " +
                                                                  std::to_string(body->GetID()) +
                                                                  " for a wave direction of " + std::to_string(
                                                                      angle) +
                                                                  " deg.");
      angleGroup.getDataSet("Phase").createAttribute<std::string>("Unit", "rad");

//      H5Easy::dump(HDF5_file, anglePath + "/RealCoeffs", static_cast<Eigen::MatrixXd>(coeff.real()));
//      angleGroup.getDataSet("RealCoeffs").createAttribute<std::string>("Description",
//                                                                       "Real part of the Froude-Krylov loads on body " +
//                                                                       std::to_string(body->GetID()) +
//                                                                       " for a wave direction of " + std::to_string(
//                                                                           m_waveDirectionDiscretization.GetVector()[iwaveDir]) +
//                                                                       " deg.");
//      angleGroup.getDataSet("RealCoeffs").createAttribute<std::string>("Unit", "N/m");
//
//      H5Easy::dump(HDF5_file, anglePath + "/ImagCoeffs", static_cast<Eigen::MatrixXd>(coeff.imag()));
//      angleGroup.getDataSet("ImagCoeffs").createAttribute<std::string>("Description",
//                                                                       "Imaginary part of the Froude-Krylov loads on body " +
//                                                                       std::to_string(body->GetID()) +
//                                                                       " for a wave direction of " + std::to_string(
//                                                                           m_waveDirectionDiscretization.GetVector()[iwaveDir]) +
//                                                                       " deg.");
//      angleGroup.getDataSet("ImagCoeffs").createAttribute<std::string>("Unit", "N/m");


    }

  }

  void HDBWriter::WriteWaveDrift(HighFive::File &HDF5_file) const {

    auto waveDriftGroup = HDF5_file.createGroup("WaveDrift");

    // Write symmetries data
    auto symmetries = m_hdb->GetWaveDrift()->GetSymmetries();

    auto sym_x = static_cast<unsigned int>(symmetries[0]);
    auto sym_y = static_cast<unsigned int>(symmetries[1]);

    HighFive::DataSet dataSet = waveDriftGroup.createDataSet<unsigned int>("sym_x", HighFive::DataSpace::From(sym_x));
    dataSet.write(sym_x);
    dataSet.createAttribute<std::string>("Description", "Symmetry along x");

    dataSet = waveDriftGroup.createDataSet<unsigned int>("sym_y", HighFive::DataSpace::From(sym_y));
    dataSet.write(sym_y);
    dataSet.createAttribute<std::string>("Description", "Symmetry along y");

    // Write wave drift data
    std::vector<std::string> components = {"surge", "sway", "yaw"};

    auto angles = m_hdb->GetWaveDirectionDiscretization();
    auto frequencies = m_hdb->GetFrequencyDiscretization();
    auto waveDrift = m_hdb->GetWaveDrift();

    for (auto & component: components) {

      auto componentGroup = waveDriftGroup.createGroup(component);

      for (unsigned int iangle = 0; iangle < angles.size(); iangle++) {
        auto angle = angles(iangle);
        auto headingGroup = componentGroup.createGroup("angle_" + std::to_string(iangle));
        Eigen::VectorXd data(frequencies.size());
        for (unsigned int i=0; i<frequencies.size(); i++) {
          data(i) = waveDrift->Eval(component, frequencies(i), angle);
        }
        H5Easy::dump(HDF5_file, "WaveDrift/"+component+"/angle_"+std::to_string(iangle)+"/data", data);
        dataSet = headingGroup.getDataSet("data");
        dataSet.createAttribute<std::string>("Description", "Wave Drift force coefficients");

        dataSet = headingGroup.createDataSet<double>("angle", HighFive::DataSpace::From(angle));
        dataSet.write(angle);
        dataSet.createAttribute<std::string>("Description", "Wave direction angle");
        // TODO : angle unit gestion to add
        dataSet.createAttribute<std::string>("Unit", "rad");

      }
    }

  }

  void export_HDB(const std::string &filename, HydrodynamicDataBase *hdb) {

    auto hdbWriter = std::make_shared<HDBWriter>(hdb);
    hdbWriter->Write(filename);

  }


}
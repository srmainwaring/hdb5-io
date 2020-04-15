//
// Created by lletourn on 14/04/20.
//

#include "HDBReader.h"

#include "../HydrodynamicDataBase.h"
#include "../Body.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

namespace HDB5_io {

  void HDBReader::Read(const std::string &filename) {
    HighFive::File file(filename, HighFive::File::ReadOnly);

    ReadHDBBasics(file);

    ReadDiscretizations(file);

    std::vector<Body*> bodies;

    for (int i = 0; i < m_hdb->GetNbBodies(); i++) {
      bodies.push_back(ReadBodyBasics(file, "Bodies/Body_" + std::to_string(i)));
    }


    for (auto &body : bodies) {

      ReadMesh(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Mesh", body);

      ReadExcitation(Diffraction, file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Excitation/Diffraction",
                     body);
      ReadExcitation(Froude_Krylov, file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Excitation/FroudeKrylov",
                     body);

      body->ComputeExcitation();

      ReadRadiation(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Radiation", body);

      if (file.getGroup("Bodies/Body_" + std::to_string(body->GetID())).exist("RAO")) {
        ReadRAO(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/RAO", body);
      }
    }

    if (file.exist("WaveDrift")) {
      ReadWaveDrift(file);
    }
  }

  void HDBReader::ReadHDBBasics(const HighFive::File &HDF5_file) {

    m_hdb->SetCreationDate(H5Easy::load<std::string>(HDF5_file, "CreationDate"));
    m_hdb->SetSolver(H5Easy::load<std::string>(HDF5_file, "Solver"));
    m_hdb->SetNormalizationLength(H5Easy::load<double>(HDF5_file, "NormalizationLength"));
    m_hdb->SetGravityAcceleration(H5Easy::load<double>(HDF5_file, "GravityAcc"));
    m_hdb->SetWaterDensity(H5Easy::load<double>(HDF5_file, "WaterDensity"));
    m_hdb->SetWaterDepth(H5Easy::load<double>(HDF5_file, "WaterDepth"));
    m_hdb->SetNbBodies(H5Easy::load<int>(HDF5_file, "NbBody"));

  }

  Body* HDBReader::ReadBodyBasics(const HighFive::File &file, const std::string &path) {

//    std::string name;
//    file.getGroup(path).getDataSet("BodyName").read(name);
    auto name = H5Easy::load<std::string>(file, path + "/BodyName");
    auto id = H5Easy::load<unsigned int>(file, path + "/ID");

    auto body = m_hdb->NewBody(id, name);

    body->SetPosition(H5Easy::load<Eigen::Vector3d>(file, path + "/BodyPosition"));

    body->SetForceMask(H5Easy::load<Eigen::Matrix<int, 6, 1>>(file, path + "/Mask/ForceMask"));
    body->SetMotionMask(H5Easy::load<Eigen::Matrix<int, 6, 1>>(file, path + "/Mask/MotionMask"));

    if (file.exist(path + "Hydrostatic")) {
      mathutils::Matrix66<double> stiffnessMatrix ;
      stiffnessMatrix = H5Easy::load<Eigen::Matrix<double, 6, 6>>(file, path + "/Hydrostatic/StiffnessMatrix");
      body->SetStiffnessMatrix(stiffnessMatrix);
    }
    if (file.exist(path + "Inertia")) {
      mathutils::Matrix66<double> inertiaMatrix ;
      inertiaMatrix = H5Easy::load<Eigen::Matrix<double, 6, 6>>(file, path + "/Inertia/InertiaMatrix");
      body->SetStiffnessMatrix(inertiaMatrix);
    }

    return body;

  }

  void
  HDBReader::ReadExcitation(HDBReader::excitationType type, const HighFive::File &HDF5_file, const std::string &path,
                            Body *body) {
    auto forceMask = body->GetForceMask();

    for (unsigned int iwaveDir = 0; iwaveDir < m_hdb->GetWaveDirectionDiscretization().size(); ++iwaveDir) {

      auto WaveDirPath = path + "/Angle_" + std::to_string(iwaveDir);

//      auto angle = H5Easy::load<double>(HDF5_file, WaveDirPath + "/Angle");
//      assert(abs(m_waveDirectionDiscretization.GetVector()[iwaveDir] - angle) < 1E-5);

      auto realCoeffs = H5Easy::load<Eigen::MatrixXd>(HDF5_file, WaveDirPath + "/RealCoeffs");
      auto imagCoeffs = H5Easy::load<Eigen::MatrixXd>(HDF5_file, WaveDirPath + "/ImagCoeffs");
      auto Dcoeffs = realCoeffs + MU_JJ * imagCoeffs;

      Eigen::MatrixXcd ExcitationCoeffs;
      if (imagCoeffs.rows() != 6) {
        assert(imagCoeffs.rows() == forceMask.GetNbDOF());
        ExcitationCoeffs.setZero();
        for (int i = 0; i < forceMask.GetNbDOF(); i++) {
          ExcitationCoeffs.row(forceMask.GetListDOF()[i]) = Dcoeffs.row(i);
        }
//        // Condense the matrix by removing the lines corresponding to the masked DOFs
//        ExcitationCoeffs = Eigen::VectorXi::Map(forceMask.GetListDOF().data(), forceMask.GetNbDOF()).replicate(1,Dcoeffs.cols()).unaryExpr(Dcoeffs);
      } else {
        ExcitationCoeffs = Dcoeffs;
      }

      switch (type) {
        case Diffraction : {
          body->SetDiffraction(iwaveDir, ExcitationCoeffs);
          break;
        }
        case Froude_Krylov : {
          body->SetFroudeKrylov(iwaveDir, ExcitationCoeffs);
          break;
        }
      }

    }

  }

  void HDBReader::ReadRadiation(const HighFive::File &HDF5_file, const std::string &path, Body *body) {

    for (unsigned int ibodyMotion = 0; ibodyMotion < m_hdb->GetNbBodies(); ++ibodyMotion) {

      auto bodyMotion = m_hdb->GetBody(ibodyMotion);
      auto bodyMotionPath = path + "/BodyMotion_" + std::to_string(ibodyMotion);

      // Reading the infinite added mass matrix for the body.
      auto infiniteAddedMass = H5Easy::load<Eigen::MatrixXd>(HDF5_file, bodyMotionPath + "/InfiniteAddedMass");
      body->SetInfiniteAddedMass(bodyMotion, infiniteAddedMass);

      // Reading the radiation mask matrix for the body.
      auto radiationMask = H5Easy::load<Eigen::MatrixXi>(HDF5_file, bodyMotionPath + "/RadiationMask");
      body->SetRadiationMask(bodyMotion, radiationMask);

      // Reading the impulse response functions.
      auto impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionK",
                                                      radiationMask);
      body->SetHDBInterpolator(Body::interpolatedData::IRF_K, bodyMotion, impulseResponseFunctionsK);

      impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionKU",
                                                 radiationMask);
      body->SetHDBInterpolator(Body::interpolatedData::IRF_KU, bodyMotion, impulseResponseFunctionsK);

      // Reading the added mass and radiation damping coefficients
      auto addedMass = ReadComponents(HDF5_file, bodyMotionPath + "/AddedMass", radiationMask);
      body->SetHDBInterpolator(Body::interpolatedData::ADDED_MASS, bodyMotion, addedMass);

      auto radiationDamping = ReadComponents(HDF5_file, bodyMotionPath + "/RadiationDamping", radiationMask);
      body->SetHDBInterpolator(Body::interpolatedData::RADIATION_DAMPING, bodyMotion, radiationDamping);

    }


  }

  void HDBReader::ReadRAO(const HighFive::File &HDF5_file, const std::string &path, Body *body) {
    auto forceMask = body->GetForceMask();

    for (unsigned int iwaveDir = 0; iwaveDir < m_hdb->GetWaveDirectionDiscretization().size(); ++iwaveDir) {

      auto WaveDirPath = path + "/Angle_" + std::to_string(iwaveDir);

//      auto angle = H5Easy::load<double>(HDF5_file, WaveDirPath + "/Angle");
//      assert(abs(m_waveDirectionDiscretization.GetVector()[iwaveDir] - angle) < 1E-5);

      auto amplitude = H5Easy::load<Eigen::MatrixXd>(HDF5_file, WaveDirPath + "/Amplitude");
      auto phase = H5Easy::load<Eigen::MatrixXd>(HDF5_file, WaveDirPath + "/Phase");
//      auto Dcoeffs = amplitude;
      Eigen::MatrixXcd Dcoeffs = amplitude.array() * Eigen::exp( MU_JJ * phase.array() * MU_PI / 180.);

      Eigen::MatrixXcd raoCoeffs;
      if (amplitude.rows() != 6) {
        assert(amplitude.rows() == forceMask.GetNbDOF());
        raoCoeffs.setZero();
        for (int i = 0; i < forceMask.GetNbDOF(); i++) {
          raoCoeffs.row(forceMask.GetListDOF()[i]) = Dcoeffs.row(i);
        }
//        // Condense the matrix by removing the lines corresponding to the masked DOFs
//        raoCoeffs = Eigen::VectorXi::Map(forceMask.GetListDOF().data(), forceMask.GetNbDOF()).replicate(1,Dcoeffs.cols()).unaryExpr(Dcoeffs);
      } else {
        raoCoeffs = Dcoeffs;
      }

      body->SetRAO(iwaveDir, raoCoeffs);

    }


  }

  std::vector<Eigen::MatrixXd>
  HDBReader::ReadComponents(const HighFive::File &HDF5_file, const std::string &path, Eigen::MatrixXi radiationMask) {

    std::vector<Eigen::MatrixXd> impulseResponseFunctionsK;

    Mask motionMask;
    Eigen::MatrixXd IRFCoeffs;

    for (unsigned int imotion = 0; imotion < 6; ++imotion) {
      auto IRF = H5Easy::load<Eigen::MatrixXd>(HDF5_file, path + "/DOF_" + std::to_string(imotion));
      motionMask.SetMask(radiationMask.row(imotion));
      if (IRF.rows() != 6) {
        assert(IRF.rows() == motionMask.GetNbDOF());
        IRFCoeffs.setZero();
        for (int i = 0; i < motionMask.GetNbDOF(); i++) {
          IRFCoeffs.row(motionMask.GetListDOF()[i]) = IRF.row(i);
        }
//          // Condense the matrix by removing the lines corresponding to the masked DOFs
//          //TODO:: passer en fonction de MathUtils ?
//          IRFCoeffs = Eigen::VectorXi::Map(motionMask.GetListDOF().data(), motionMask.GetNbDOF()).replicate(1,IRF.cols()).unaryExpr(IRF);
      } else {
        IRFCoeffs = IRF;
      }
      impulseResponseFunctionsK.push_back(IRFCoeffs);
    }

    return impulseResponseFunctionsK;

  }

  void HDBReader::ReadMesh(HighFive::File &HDF5_file, const std::string &path, Body *body) {

    auto nbVertices = H5Easy::load<int>(HDF5_file, path + "/NbVertices");
    auto nbFaces = H5Easy::load<int>(HDF5_file, path + "/NbFaces");

    auto vertices_hdb = H5Easy::load<Eigen::MatrixXd>(HDF5_file, path + "/Vertices");
    auto faces_hdb = H5Easy::load<Eigen::MatrixXi>(HDF5_file, path + "/Faces");

    std::vector<mathutils::Vector3d<double>> vertices;
    std::vector<Eigen::VectorXi> faces;

    for (unsigned int i = 0; i < nbVertices; i++) {
      mathutils::Vector3d<double> vertex = vertices_hdb.row(i);
      vertices.emplace_back(vertex);
    }

    for (unsigned int i = 0; i < nbFaces; i++) {
      Eigen::VectorXi face = faces_hdb.row(i);
      faces.emplace_back(face);
    }

    body->LoadMesh(vertices, faces);

  }

  void HDBReader::ReadWaveDrift(HighFive::File &HDF5_file) {

    auto waveDrift = std::make_shared<WaveDrift>();

//    auto frequency = H5Easy::load<Eigen::VectorXd>(HDF5_file, "WaveDrift/freq");
    auto frequency = m_hdb->GetFrequencyDiscretization();
    auto waveDirection = m_hdb->GetWaveDirectionDiscretization();
//    assert(frequency == GetFrequencyDiscretization().GetVectorN());

    waveDrift->SetFrequencies(m_hdb->GetFrequencyDiscretization());
    waveDrift->SetWaveDirections(waveDirection);

    auto sym_X = H5Easy::load<int>(HDF5_file, "WaveDrift/sym_x");
    auto sym_Y = H5Easy::load<int>(HDF5_file, "WaveDrift/sym_y");

    waveDrift->SetSymmetries(sym_X == 1, sym_Y == 1);

    Eigen::MatrixXd surge(waveDirection.size(), frequency.size());
    Eigen::MatrixXd sway(waveDirection.size(), frequency.size());
    Eigen::MatrixXd yaw(waveDirection.size(), frequency.size());

    for (unsigned int i = 0; i < waveDirection.size(); i++) {
      auto data = H5Easy::load<Eigen::VectorXd>(HDF5_file, "WaveDrift/surge/heading_" + std::to_string(i) + "/data");
      surge.row(i) = data;
      data = H5Easy::load<Eigen::VectorXd>(HDF5_file, "WaveDrift/sway/heading_" + std::to_string(i) + "/data");
      sway.row(i) = data;
      data = H5Easy::load<Eigen::VectorXd>(HDF5_file, "WaveDrift/yaw/heading_" + std::to_string(i) + "/data");
      yaw.row(i) = data;
    }

    std::vector<double> coeff_surge(&surge(0,0), surge.data()+surge.size());
    waveDrift->AddData("surge", coeff_surge);
    std::vector<double> coeff_sway(&sway(0,0), sway.data()+sway.size());
    waveDrift->AddData("sway", coeff_sway);
    std::vector<double> coeff_yaw(&yaw(0,0), yaw.data()+yaw.size());
    waveDrift->AddData("yaw", coeff_yaw);

    m_hdb->SetWaveDrift(waveDrift);

  }


  std::shared_ptr<HydrodynamicDataBase> import_HDB(const std::string& filename){
    auto hdb = std::make_shared<HydrodynamicDataBase>();

    HighFive::File file(filename, HighFive::File::ReadOnly);

    double version = 1.0;
    if (file.exist("Version"))
      version = H5Easy::load<double>(file, "Version");

    if (version<=2) {
      auto hdb_reader = std::make_shared<HDBReader_v2>(hdb.get());
      hdb_reader->Read(filename);
    }
    else {
      auto hdb_reader = std::make_shared<HDBReader_v3>(hdb.get());
      hdb_reader->Read(filename);
    }


    return hdb;
  }

  void HDBReader_v2::ReadDiscretizations(const HighFive::File &file) {

    double min, max;
    unsigned int nb;
    auto disc = file.getGroup("Discretizations").getGroup("Frequency");
    disc.getDataSet("MinFrequency").read(min);
    disc.getDataSet("MaxFrequency").read(max);
    disc.getDataSet("NbFrequencies").read(nb);
    m_hdb->SetFrequencyDiscretization(mathutils::VectorN<double>::LinSpaced(nb, min, max));

    disc = file.getGroup("Discretizations").getGroup("Time");
    disc.getDataSet("TimeStep").read(min);
    disc.getDataSet("FinalTime").read(max);
    disc.getDataSet("NbTimeSample").read(nb);
//    assert(abs(max/double(nb) - min) < 1E-5);
    m_hdb->SetTimeDiscretization(mathutils::VectorN<double>::LinSpaced(nb, 0., max));

    disc = file.getGroup("Discretizations").getGroup("WaveDirections");
    disc.getDataSet("MinAngle").read(min);
    disc.getDataSet("MaxAngle").read(max);
    disc.getDataSet("NbWaveDirections").read(nb);
    m_hdb->SetWaveDirectionDiscretization(mathutils::VectorN<double>::LinSpaced(nb, min, max));
  }



  void HDBReader_v3::ReadDiscretizations(const HighFive::File &file) {

    m_hdb->SetFrequencyDiscretization(H5Easy::load<Eigen::VectorXd>(file, "Discretizations/Frequency"));
    m_hdb->SetTimeDiscretization(H5Easy::load<Eigen::VectorXd>(file, "Discretizations/Time"));
    m_hdb->SetWaveDirectionDiscretization(H5Easy::load<Eigen::VectorXd>(file, "Discretizations/WaveDirection"));

  }


} // end namespace HDB5_io
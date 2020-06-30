//
// Created by lletourn on 14/04/20.
//

#include "HDBReader.h"

#include "hdb5_io/containers/HydrodynamicDataBase.h"
#include "hdb5_io/containers/Body.h"
#include "hdb5_io/containers/WaveDrift.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

namespace HDB5_io {

  void HDBReader::Read(const std::string &filename) {

    // HDB file.
    HighFive::File file(filename, HighFive::File::ReadOnly);

    // HDB basic data (version, date, gravity constant, density, etc.).
    ReadHDBBasics(file);

    // Wave frequency and wave directions.
    ReadDiscretizations(file);

    // Symmetries.
    if (file.exist("Symmetries")) {
      ReadSymmetries(file);
    }

    // Body basic data (index, name, position, mass, etc.).
    std::vector<Body *> bodies;
    for (int i = 0; i < m_hdb->GetNbBodies(); i++) {
      bodies.push_back(ReadBodyBasics(file, "Bodies/Body_" + std::to_string(i)));
    }

    // Other body data.
    for (auto &body : bodies) {

      // Body mesh.
      ReadMesh(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Mesh", body);

      // Diffraction loads.
      ReadExcitation(Diffraction, file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Excitation/Diffraction",
                     body);

      // Froude-Krylov loads.
      ReadExcitation(Froude_Krylov, file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Excitation/FroudeKrylov",
                     body);

      // Excitation loads.
      body->ComputeExcitation();

      // Added mass, damping, IRF, poles and residues.
      ReadRadiation(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Radiation", body);

      // RAOs.
      if (file.getGroup("Bodies/Body_" + std::to_string(body->GetID())).exist("RAO")) {
        ReadRAO(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/RAO", body);
      }
    }

    // Mean wave drift loads.
    if (file.exist("WaveDrift")) {
      ReadWaveDrift(file);
    }

    // Wave field parameters.
    if (file.exist("WaveField")) {
      ReadWaveField(file);
    }

    // Vector fitting parameters.
    if (file.exist("VectorFitting")) {
      ReadVectorFitting(file);
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

  Body *HDBReader::ReadBodyBasics(const HighFive::File &file, const std::string &path) {

    // Name.
    auto name = H5Easy::load<std::string>(file, path + "/BodyName");

    // Index.
    auto id = H5Easy::load<unsigned int>(file, path + "/ID");

    // New body.
    auto body = m_hdb->NewBody(id, name);

    // Position.
    body->SetPosition(H5Easy::load<Eigen::Vector3d>(file, path + "/BodyPosition"));

    // Force and Motion masks.
    body->SetForceMask(H5Easy::load<Eigen::Matrix<int, 6, 1>>(file, path + "/Mask/ForceMask"));
    body->SetMotionMask(H5Easy::load<Eigen::Matrix<int, 6, 1>>(file, path + "/Mask/MotionMask"));

    // Hydrostatic matrix.
    if (file.exist(path + "/Hydrostatic")) {
      mathutils::Matrix66<double> stiffnessMatrix;
      stiffnessMatrix = H5Easy::load<Eigen::Matrix<double, 6, 6>>(file, path + "/Hydrostatic/StiffnessMatrix");
      body->SetStiffnessMatrix(stiffnessMatrix);
    }

    // Inertia matrix.
    if (file.exist(path + "/Inertia")) {
      mathutils::Matrix66<double> inertiaMatrix;
      inertiaMatrix = H5Easy::load<Eigen::Matrix<double, 6, 6>>(file, path + "/Inertia/InertiaMatrix");
      body->SetInertia(inertiaMatrix);
    }

    // Mooring matrix.
    if (file.exist(path + "/Mooring")) {
      mathutils::Matrix66<double> mooringMatrix;
      mooringMatrix = H5Easy::load<Eigen::Matrix<double, 6, 6>>(file, path + "/Mooring/MooringMatrix");
      body->SetMooring(mooringMatrix);
    }

    // Damping matrix.
    if (file.exist(path + "/LinearDamping")) {
      mathutils::Matrix66<double> dampingMatrix;
      dampingMatrix = H5Easy::load<Eigen::Matrix<double, 6, 6>>(file, path + "/LinearDamping/DampingMatrix");
      body->SetLinearDamping(dampingMatrix);
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
      if (HDF5_file.exist(bodyMotionPath + "/ImpulseResponseFunctionK")) {
        auto impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionK",
                                                        radiationMask);
        body->SetHDBInterpolator(Body::interpolatedData::IRF_K, bodyMotion, impulseResponseFunctionsK);
      }

      if (HDF5_file.exist(bodyMotionPath + "/ImpulseResponseFunctionKU")) {
        auto impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionKU",
                                                   radiationMask);
        body->SetHDBInterpolator(Body::interpolatedData::IRF_KU, bodyMotion, impulseResponseFunctionsK);
      }

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

      auto DataSet = HDF5_file.getDataSet(WaveDirPath + "/Phase");
      std::string unit;
      DataSet.getAttribute("Unit").read<std::string>(unit);
      if (unit != "rad")
        phase = phase.array() * DEG2RAD;

      Eigen::MatrixXcd Dcoeffs = amplitude.array() * Eigen::exp(MU_JJ * phase.array());

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

  void HDBReader::ReadWaveField(HighFive::File &file) {

    m_hdb->SetWaveField();

  }

  void HDBReader::ReadVectorFitting(HighFive::File &file) {

    m_hdb->SetVF();
    m_hdb->SetVFRelaxed(H5Easy::load<int>(file, "VectorFitting/Relaxed"));
    m_hdb->SetVFMaxOrder(H5Easy::load<int>(file, "VectorFitting/MaxOrder"));
    m_hdb->SetVFTolerance(H5Easy::load<double>(file, "VectorFitting/Tolerance"));

  }

  void HDBReader::ReadSymmetries(HighFive::File &file) {

    m_hdb->SetSymmetries();
    m_hdb->SetSymBottom(H5Easy::load<int>(file, "Symmetries/Bottom"));
    m_hdb->SetSymXOZ(H5Easy::load<int>(file, "Symmetries/xOz"));
    m_hdb->SetSymYOZ(H5Easy::load<int>(file, "Symmetries/yOz"));

  }

  std::shared_ptr<HydrodynamicDataBase> import_HDB(const std::string &filename) {
    auto hdb = std::make_shared<HydrodynamicDataBase>();

    HighFive::File file(filename, HighFive::File::ReadOnly);

    double version = 1.0;
    if (file.exist("Version"))
      version = H5Easy::load<double>(file, "Version");

    if (version <= 2) {
      auto hdb_reader = std::make_shared<HDBReader_v2>(hdb.get());
      hdb_reader->Read(filename);
    } else {
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

    // The wave directions are written in degrees in the hdb5 file. The conversion degrees to radians is performed
    // in the method SetWaveDirectionDiscretization.
    disc = file.getGroup("Discretizations").getGroup("WaveDirections");
    disc.getDataSet("MinAngle").read(min); // In degrees.
    disc.getDataSet("MaxAngle").read(max); // In degrees.
    disc.getDataSet("NbWaveDirections").read(nb);

    // The wave directions are read in degrees from the hdb5 file. The conversion degrees to radians is performed below.
    m_hdb->SetWaveDirectionDiscretization(mathutils::VectorN<double>::LinSpaced(nb, min, max) * MU_PI_180);
  }

  Eigen::VectorXd
  HDBReader_v2::ReadWaveDriftComponents(HighFive::File &HDF5_file, const std::string &path, unsigned int i) {
    return H5Easy::load<Eigen::VectorXd>(HDF5_file, path + "/heading_" + std::to_string(i) + "/data");
  }

  void HDBReader_v2::ReadWaveDrift(HighFive::File &HDF5_file) {

    auto waveDrift = std::make_shared<WaveDrift>();

    auto frequency = H5Easy::load<Eigen::VectorXd>(HDF5_file, "WaveDrift/freq");
//    auto frequency = m_hdb->GetFrequencyDiscretization();
    auto waveDirection = m_hdb->GetWaveDirectionDiscretization();
//    assert(frequency == GetFrequencyDiscretization().GetVectorN());

    waveDrift->SetFrequencies(frequency);
    waveDrift->SetWaveDirections(waveDirection);

    auto sym_X = H5Easy::load<int>(HDF5_file, "WaveDrift/sym_x");
    auto sym_Y = H5Easy::load<int>(HDF5_file, "WaveDrift/sym_y");

    waveDrift->SetSymmetries(sym_X == 1, sym_Y == 1);

    Eigen::MatrixXd surge(waveDirection.size(), frequency.size());
    Eigen::MatrixXd sway(waveDirection.size(), frequency.size());
    Eigen::MatrixXd yaw(waveDirection.size(), frequency.size());

    for (unsigned int i = 0; i < waveDirection.size(); i++) {
      if (HDF5_file.exist("WaveDrift/surge")) {
        auto data_surge = ReadWaveDriftComponents(HDF5_file, "WaveDrift/surge", i);
        surge.row(i) = data_surge;
      }
      if (HDF5_file.exist("WaveDrift/sway")) {
        auto data_sway = ReadWaveDriftComponents(HDF5_file, "WaveDrift/sway", i);
        sway.row(i) = data_sway;
      }
      if (HDF5_file.exist("WaveDrift/yaw")) {
        auto data_yaw = ReadWaveDriftComponents(HDF5_file, "WaveDrift/yaw", i);
        yaw.row(i) = data_yaw;
      }
    }

    std::vector<double> coeff_surge(&surge(0, 0), surge.data() + surge.size());
    waveDrift->AddData("surge", coeff_surge);
    std::vector<double> coeff_sway(&sway(0, 0), sway.data() + sway.size());
    waveDrift->AddData("sway", coeff_sway);
    std::vector<double> coeff_yaw(&yaw(0, 0), yaw.data() + yaw.size());
    waveDrift->AddData("yaw", coeff_yaw);

    m_hdb->SetWaveDrift(waveDrift);

  }

  Eigen::VectorXd
  HDBReader_v3::ReadWaveDriftComponents(HighFive::File &HDF5_file, const std::string &path, unsigned int i) {
    return H5Easy::load<Eigen::VectorXd>(HDF5_file, path + "/angle_" + std::to_string(i) + "/data");
  }

  void HDBReader_v3::ReadWaveDrift(HighFive::File &HDF5_file) {

    auto waveDrift = std::make_shared<WaveDrift>();

//    auto frequency = H5Easy::load<Eigen::VectorXd>(HDF5_file, "WaveDrift/freq");
    auto frequency = m_hdb->GetFrequencyDiscretization();
    auto waveDirection = m_hdb->GetWaveDirectionDiscretization();
//    assert(frequency == GetFrequencyDiscretization().GetVectorN());

    waveDrift->SetFrequencies(m_hdb->GetFrequencyDiscretization());
    waveDrift->SetWaveDirections(waveDirection);

    auto kochin_step = H5Easy::load<double>(HDF5_file, "WaveDrift/KochinStep"); // In degree.
    auto sym_X = H5Easy::load<int>(HDF5_file, "WaveDrift/sym_x");
    auto sym_Y = H5Easy::load<int>(HDF5_file, "WaveDrift/sym_y");

    waveDrift->SetSymmetries(sym_X == 1, sym_Y == 1);
    waveDrift->SetKochinStep(kochin_step * MU_PI_180); // Conversion in radians.

    Eigen::MatrixXd surge(waveDirection.size(), frequency.size());
    Eigen::MatrixXd sway(waveDirection.size(), frequency.size());
    Eigen::MatrixXd yaw(waveDirection.size(), frequency.size());

    for (unsigned int i = 0; i < waveDirection.size(); i++) {
      auto data_surge = ReadWaveDriftComponents(HDF5_file, "WaveDrift/surge", i);
      surge.row(i) = data_surge;
      auto data_sway = ReadWaveDriftComponents(HDF5_file, "WaveDrift/sway", i);
      sway.row(i) = data_sway;
      auto data_yaw = ReadWaveDriftComponents(HDF5_file, "WaveDrift/yaw", i);
      yaw.row(i) = data_yaw;
    }

    std::vector<double> coeff_surge(&surge(0, 0), surge.data() + surge.size());
    waveDrift->AddData("surge", coeff_surge);
    std::vector<double> coeff_sway(&sway(0, 0), sway.data() + sway.size());
    waveDrift->AddData("sway", coeff_sway);
    std::vector<double> coeff_yaw(&yaw(0, 0), yaw.data() + yaw.size());
    waveDrift->AddData("yaw", coeff_yaw);

    m_hdb->SetWaveDrift(waveDrift);

  }

  void HDBReader_v3::ReadDiscretizations(const HighFive::File &file) {

    m_hdb->SetFrequencyDiscretization(H5Easy::load<Eigen::VectorXd>(file, "Discretizations/Frequency"));
    m_hdb->SetTimeDiscretization(H5Easy::load<Eigen::VectorXd>(file, "Discretizations/Time"));

    // The wave directions are read in degrees from the hdb5 file. The conversion degrees to radians is performed below.
    m_hdb->SetWaveDirectionDiscretization(H5Easy::load<Eigen::VectorXd>(file, "Discretizations/WaveDirection") * MU_PI_180);

  }

  void HDBReader_v3::ReadRadiation(const HighFive::File &file, const std::string &path, Body *body) {
    HDBReader::ReadRadiation(file, path, body);

    for (unsigned int ibodyMotion = 0; ibodyMotion < m_hdb->GetNbBodies(); ++ibodyMotion) {

      auto bodyMotion = m_hdb->GetBody(ibodyMotion);
      auto bodyMotionPath = path + "/BodyMotion_" + std::to_string(ibodyMotion);

      if (file.exist(bodyMotionPath + "/ZeroFreqAddedMass")) {
        auto zeroFreqAddedMass = H5Easy::load<Eigen::MatrixXd>(file, bodyMotionPath + "/ZeroFreqAddedMass");
        body->SetZeroFreqAddedMass(bodyMotion, zeroFreqAddedMass);
      }

      if (file.exist(bodyMotionPath + "/Modal")) {

        for (unsigned int idof = 0; idof < 6; idof++) {

          std::vector<PoleResidue> modalCoeff;
          for (unsigned int iforce = 0; iforce < 6; iforce++) {

            auto forcePath = bodyMotionPath + "/Modal/DOF_" + std::to_string(idof) + "/FORCE_" + std::to_string(iforce);

            // Real poles and residues.
            int nPoles_real = 0;
            Eigen::VectorXd poles, residues;
            if(file.exist(forcePath + "/RealPoles")){
              poles = H5Easy::load<Eigen::VectorXd>(file, forcePath + "/RealPoles");
              nPoles_real = poles.size();
              residues = H5Easy::load<Eigen::VectorXd>(file, forcePath + "/RealResidues");
              assert(residues.size() == nPoles_real);
            }

            // Complex poles and residues.
            int nPoles_cc = 0;
            Eigen::VectorXcd cplxPoles, cplxResidues;
            if(file.exist(forcePath + "/ComplexPoles/RealCoeff")) {
              // Poles.
              auto realCoeff = H5Easy::load<Eigen::VectorXd>(file, forcePath + "/ComplexPoles/RealCoeff");
              auto imagCoeff = H5Easy::load<Eigen::VectorXd>(file, forcePath + "/ComplexPoles/ImagCoeff");
              nPoles_cc = realCoeff.size();
              assert(imagCoeff.size() == nPoles_cc);
              cplxPoles = realCoeff + MU_JJ * imagCoeff;

              // Residues.
              realCoeff = H5Easy::load<Eigen::VectorXd>(file, forcePath + "/ComplexResidues/RealCoeff");
              imagCoeff = H5Easy::load<Eigen::VectorXd>(file, forcePath + "/ComplexResidues/ImagCoeff");
              assert(realCoeff.size() == nPoles_cc && imagCoeff.size() == nPoles_cc);
              cplxResidues = realCoeff + MU_JJ * imagCoeff;
            }

            // Adding to modalCoeff.
            PoleResidue pair;
            for (int i=0; i < nPoles_real; i++) {
              pair.AddPoleResidue(poles(i), residues(i));
            }
            for (int i=0; i < nPoles_cc; i++) {
              pair.AddPoleResidue(cplxPoles(i), cplxResidues(i));
            }
            modalCoeff.emplace_back(pair);

          }
          body->AddModalCoefficients(bodyMotion, modalCoeff);

        }

      }

    }

  }


} // end namespace HDB5_io
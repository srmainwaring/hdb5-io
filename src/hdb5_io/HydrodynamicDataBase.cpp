//
// Created by lletourn on 26/02/20.
//

#include <Eigen/Dense>

#include "HydrodynamicDataBase.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

namespace HDB5_io {

  void HydrodynamicDataBase::SetCreationDate(std::string date) {
    m_creationDate = date;
  }

  std::string HydrodynamicDataBase::GetCreationDate() const {
    return m_creationDate;
  }

  void HydrodynamicDataBase::SetSolver(std::string solver) {
    m_solver = solver;
  }

  std::string HydrodynamicDataBase::GetSolver() const {
    return m_solver;
  }

  void HydrodynamicDataBase::SetVersion(double version) {
    m_version = version;
  }

  double HydrodynamicDataBase::GetVersion() const {
    return m_version;
  }

  void HydrodynamicDataBase::SetGravityAcceleration(double g) {
    m_gravityAcceleration = g;
  }

  double HydrodynamicDataBase::GetGravityAcceleration() const {
    return m_gravityAcceleration;
  }

  void HydrodynamicDataBase::SetWaterDensity(double rho) {
    m_waterDensity = rho;
  }

  double HydrodynamicDataBase::GetWaterDensity() const {
    return m_waterDensity;
  }

  void HydrodynamicDataBase::SetWaterDepth(double h) {
    m_waterDepth = h;
  }

  double HydrodynamicDataBase::GetWaterDepth() const {
    return m_waterDepth;
  }

  void HydrodynamicDataBase::SetNormalizationLength(double L) {
    m_normalizationLength = L;
  }

  double HydrodynamicDataBase::GetNormalizationLength() const {
    return m_normalizationLength;
  }

  Body *HydrodynamicDataBase::GetBody(int id) const {
    return m_bodies[id].get();
  }

  Body *HydrodynamicDataBase::NewBody(unsigned int id, const std::string &name) {
    m_bodies.push_back(std::make_unique<Body>(id, name, this));
    return (m_bodies.back()).get();
  }

  int HydrodynamicDataBase::GetNbBodies() const {
    return m_nbody;
  }

  void HydrodynamicDataBase::SetFrequencyDiscretization(double wmin, double wmax, unsigned int nw) {
    m_frequencyDiscretization = Discretization1D(wmin, wmax, nw);
  }

  void HydrodynamicDataBase::SetWaveDirectionDiscretization(double tmin, double tmax, unsigned int nt) {
    m_waveDirectionDiscretization = Discretization1D(tmin, tmax, nt);
  }

  void HydrodynamicDataBase::SetTimeDiscretization(double tmin, double tmax, unsigned int nt) {
    m_timeDiscretization = Discretization1D(tmin, tmax, nt);
  }

  Discretization1D HydrodynamicDataBase::GetFrequencyDiscretization() const {
    return m_frequencyDiscretization;
  }

  Discretization1D HydrodynamicDataBase::GetWaveDirectionDiscretization() const {
    return m_waveDirectionDiscretization;
  }

  Discretization1D HydrodynamicDataBase::GetTimeDiscretization() const {
    return m_timeDiscretization;
  }


  void HydrodynamicDataBase::Import_HDF5(const std::string &HDF5_file) {

    HighFive::File file(HDF5_file, HighFive::File::ReadOnly);

    try {
      file.getDataSet("Version").read(m_version);
    } catch (HighFive::Exception &err) {
      std::cerr << err.what() << std::endl;
      m_version = 1.0;
    }

    Import_HDF5_v3(HDF5_file);

  }

  void HydrodynamicDataBase::Import_HDF5_v3(const std::string &HDF5_file) {

    HighFive::File file(HDF5_file, HighFive::File::ReadOnly);

    file.getDataSet("CreationDate").read(m_creationDate);
    file.getDataSet("Solver").read(m_solver);
    file.getDataSet("NbBody").read(m_nbody);
    file.getDataSet("NormalizationLength").read(m_normalizationLength);
    file.getDataSet("GravityAcc").read(m_gravityAcceleration);
    file.getDataSet("WaterDensity").read(m_waterDensity);
    file.getDataSet("WaterDepth").read(m_waterDepth);

    double min, max;
    unsigned int nb;
    auto disc = file.getGroup("Discretizations").getGroup("Frequency");
    disc.getDataSet("MinFrequency").read(min);
    disc.getDataSet("MaxFrequency").read(max);
    disc.getDataSet("NbFrequencies").read(nb);
    m_frequencyDiscretization = {min, max, nb};

    disc = file.getGroup("Discretizations").getGroup("Time");
    disc.getDataSet("TimeStep").read(min);
    disc.getDataSet("FinalTime").read(max);
    disc.getDataSet("NbTimeSample").read(nb);
//    assert(abs(max/double(nb) - min) < 1E-5);
    m_timeDiscretization = {0., max, nb};

    disc = file.getGroup("Discretizations").getGroup("WaveDirections");
    disc.getDataSet("MinAngle").read(min);
    disc.getDataSet("MaxAngle").read(max);
    disc.getDataSet("NbWaveDirections").read(nb);
    m_waveDirectionDiscretization = {min, max, nb};


    for (int i = 0; i < m_nbody; i++) {
      std::string name;
      unsigned int id;

      auto hdb_body = file.getGroup("Bodies").getGroup("Body_" + std::to_string(i));

      hdb_body.getDataSet("BodyName").read(name);
      hdb_body.getDataSet("ID").read(id);

      auto body = NewBody(id, name);

//      std::vector<double> position;
//      hdb_body.getDataSet("BodyPosition").read(position);
//      body->SetPosition(Eigen::Vector3d(position.data()));

//      Eigen::Vector3d position;
//      hdb_body.getDataSet("BodyPosition").read(position);
      mathutils::Vector3d<double> position;
      position = H5Easy::load<Eigen::Vector3d>(file, "Bodies/Body_" + std::to_string(i) + "/BodyPosition");
      body->SetPosition(position);

//      std::vector<int> mask;
//      hdb_body.getGroup("Mask").getDataSet("ForceMask").read(mask);
//      body->SetForceMask(Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(mask.data(), mask.size()));
//      hdb_body.getGroup("Mask").getDataSet("MotionMask").read(mask);
//      body->SetMotionMask(Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(mask.data(), mask.size()));

      mathutils::Vector6d<int> mask;
      mask = H5Easy::load<Eigen::Matrix<int, 6, 1>>(file, "Bodies/Body_" + std::to_string(i) + "/Mask/ForceMask");
      body->SetForceMask(mask);
      mask = H5Easy::load<Eigen::Matrix<int, 6, 1>>(file, "Bodies/Body_" + std::to_string(i) + "/Mask/MotionMask");
      body->SetMotionMask(mask);

      ReadMesh(file, "Bodies/Body_" + std::to_string(i) + "/Mesh", body);
    }


    for (auto &body : m_bodies) {

      ReadExcitation(Diffraction, file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Excitation/Diffraction",
                     body.get());
      ReadExcitation(Froude_Krylov, file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Excitation/FroudeKrylov",
                     body.get());

      body->ComputeExcitation();

      ReadRadiation(file, "Bodies/Body_" + std::to_string(body->GetID()) + "/Radiation", body.get());
    }

  }

//  void HydrodynamicDataBase::ReadExcitation(HighFive::Group* group, Body* body) {
//
//    auto forceMask = body->GetForceMask();
//
//    for (unsigned int iwaveDir = 0; iwaveDir < m_waveDirectionDiscretization.GetNbSample(); ++iwaveDir) {
//
//      auto angle_group = group->getGroup("Excitation/Diffraction/Angle_" + std::to_string(iwaveDir));
//
//      double angle;
//      angle_group.getDataSet("Angle").read(angle);
//      assert(abs(m_waveDirectionDiscretization.GetVector()[iwaveDir] - angle) < 1E-5);
//
//      std::vector<std::vector<double>> coeffs;
//      angle_group.getDataSet("ImagCoeffs").read(coeffs);
//      std::cout<<coeffs[0][1]<<std::endl;
//
//      for (unsigned int idof = 0; idof < 6; idof++){}
//
//
//    }
//
//  }

  void HydrodynamicDataBase::ReadExcitation(excitationType type, const HighFive::File &HDF5_file,
                                            const std::string &path, Body *body) {
    auto forceMask = body->GetForceMask();

    for (unsigned int iwaveDir = 0; iwaveDir < m_waveDirectionDiscretization.GetNbSample(); ++iwaveDir) {

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

  void HydrodynamicDataBase::ReadRadiation(const HighFive::File &HDF5_file, const std::string &path, Body *body) {

    for (unsigned int ibodyMotion = 0; ibodyMotion < m_nbody; ++ibodyMotion) {

      auto bodyMotion = this->GetBody(ibodyMotion);
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
      body->SetImpulseResponseFunctionK(bodyMotion, impulseResponseFunctionsK);

      impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionKU",
                                                 radiationMask);
      body->SetImpulseResponseFunctionKu(bodyMotion, impulseResponseFunctionsK);

      // Reading the added mass and radiation damping coefficients
      auto addedMass = ReadComponents(HDF5_file, bodyMotionPath + "/AddedMass", radiationMask);
      body->SetAddedMass(bodyMotion, addedMass);

      auto radiationDamping = ReadComponents(HDF5_file, bodyMotionPath + "/RadiationDamping", radiationMask);
      body->SetRadiationDamping(bodyMotion, radiationDamping);

    }

  }

  std::vector<Eigen::MatrixXd> HydrodynamicDataBase::ReadComponents(const HighFive::File &HDF5_file,
                                                                    const std::string &path,
                                                                    Eigen::MatrixXi radiationMask) {

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

  void HydrodynamicDataBase::ReadMesh(HighFive::File &HDF5_file, const std::string &path, Body *body) {
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

  void HydrodynamicDataBase::Export_HDF5(const std::string &HDF5_file) {

    using namespace HighFive;

    File file(HDF5_file, File::ReadWrite | File::Create | File::Truncate);

    DataSet dataSet = file.createDataSet<std::string>("CreationDate", DataSpace::From(m_creationDate));
    dataSet.write(m_creationDate);
    dataSet.createAttribute<std::string>("Description", "Date of the creation of this database.");

    dataSet = file.createDataSet<double>("Version", DataSpace::From(m_version));
    dataSet.write(m_version);
    dataSet.createAttribute<std::string>("Description", "Version of the hdb5 output file.");

    dataSet = file.createDataSet<std::string>("Solver", DataSpace::From(m_solver));
    dataSet.write(m_solver);
    dataSet.createAttribute<std::string>("Description",
                                         "Hydrodynamic solver used for computing the hydrodynamic database.");

    dataSet = file.createDataSet<int>("NbBody", DataSpace::From(m_nbody));
    dataSet.write(m_nbody);
    dataSet.createAttribute<std::string>("Description", "Number of hydrodynamic bodies.");

    dataSet = file.createDataSet<double>("NormalizationLength", DataSpace::From(m_normalizationLength));
    dataSet.write(m_normalizationLength);
    dataSet.createAttribute<std::string>("Description", "Normalization length.");
    dataSet.createAttribute<std::string>("Unit", "m");

    dataSet = file.createDataSet<double>("GravityAcc", DataSpace::From(m_gravityAcceleration));
    dataSet.write(m_gravityAcceleration);
    dataSet.createAttribute<std::string>("Description", "Gravity acceleration.");
    dataSet.createAttribute<std::string>("Unit", "m/s**2");

    dataSet = file.createDataSet<double>("WaterDensity", DataSpace::From(m_waterDensity));
    dataSet.write(m_waterDensity);
    dataSet.createAttribute<std::string>("Description", "Water Density.");
    dataSet.createAttribute<std::string>("Unit", "kg/m**3");

    dataSet = file.createDataSet<double>("WaterDepth", DataSpace::From(m_waterDepth));
    dataSet.write(m_waterDepth);
    dataSet.createAttribute<std::string>("Description",
                                         "Water depth: 0 for infinite depth and positive for finite depth.");
    dataSet.createAttribute<std::string>("Unit", "m");

    auto discretizations = file.createGroup("Discretizations");

    auto frequency = discretizations.createGroup("Frequency");
    dataSet = frequency.createDataSet<double>("MaxFrequency", DataSpace::From(m_frequencyDiscretization.GetMax()));
    dataSet.write(m_frequencyDiscretization.GetMax());
    dataSet.createAttribute<std::string>("Description", "Maximum frequency.");
    dataSet.createAttribute<std::string>("Unit", "rad/s");

    dataSet = frequency.createDataSet<double>("MinFrequency", DataSpace::From(m_frequencyDiscretization.GetMin()));
    dataSet.write(m_frequencyDiscretization.GetMin());
    dataSet.createAttribute<std::string>("Description", "Minimum frequency.");
    dataSet.createAttribute<std::string>("Unit", "rad/s");

    dataSet = frequency.createDataSet<int>("NbFrequency", DataSpace::From(m_frequencyDiscretization.GetNbSample()));
    dataSet.write(m_frequencyDiscretization.GetNbSample());
    dataSet.createAttribute<std::string>("Description", "Number of frequencies");


    auto waveDirections = discretizations.createGroup("WaveDirections");
    dataSet = frequency.createDataSet<double>("MaxAngle", DataSpace::From(m_waveDirectionDiscretization.GetMax()));
    dataSet.write(m_waveDirectionDiscretization.GetMax());
    dataSet.createAttribute<std::string>("Description", "Maximum wave direction.");
    dataSet.createAttribute<std::string>("Unit", "deg");

    dataSet = frequency.createDataSet<double>("MinAngle", DataSpace::From(m_waveDirectionDiscretization.GetMin()));
    dataSet.write(m_waveDirectionDiscretization.GetMin());
    dataSet.createAttribute<std::string>("Description", "Minimum wave direction.");
    dataSet.createAttribute<std::string>("Unit", "deg");

    dataSet = frequency.createDataSet<int>("NbWaveDirections",
                                           DataSpace::From(m_waveDirectionDiscretization.GetNbSample()));
    dataSet.write(m_waveDirectionDiscretization.GetNbSample());
    dataSet.createAttribute<std::string>("Description", "Number of wave directions.");


    auto time = discretizations.createGroup("Time");
    dataSet = frequency.createDataSet<double>("FinalTime", DataSpace::From(m_timeDiscretization.GetMax()));
    dataSet.write(m_timeDiscretization.GetMax());
    dataSet.createAttribute<std::string>("Description",
                                         "Final time for the evaluation of the impulse response functions.");
    dataSet.createAttribute<std::string>("Unit", "s");

    dataSet = frequency.createDataSet<double>("TimeStep", DataSpace::From(m_timeDiscretization.GetStep()));
    dataSet.write(m_timeDiscretization.GetStep());
    dataSet.createAttribute<std::string>("Description", "Time step.");
    dataSet.createAttribute<std::string>("Unit", "s");

    dataSet = frequency.createDataSet<int>("NbTimeSample", DataSpace::From(m_timeDiscretization.GetNbSample()));
    dataSet.write(m_timeDiscretization.GetNbSample());
    dataSet.createAttribute<std::string>("Description", "Number of time samples.");


    auto bodies = file.createGroup("Bodies");

    for (unsigned int i = 0; i < m_nbody; i++) {
      auto bodyGroup = bodies.createGroup("Body_" + std::to_string(i));

      auto body = GetBody(i);

      dataSet = bodyGroup.createDataSet<unsigned int>("ID", DataSpace::From(i));
      dataSet.write(i);
      dataSet.createAttribute<std::string>("Description", "Body index");

      dataSet = bodyGroup.createDataSet<std::string>("BodyName", DataSpace::From(body->GetName()));
      dataSet.write(body->GetName());
      dataSet.createAttribute<std::string>("Description", "Body name");
//      dataSet = bodyGroup.createDataSet<double>("BodyPosition", DataSpace::From(static_cast<Eigen::Vector3d>(body->GetPosition())));
//      dataSet.write(body->GetPosition());
//      dataSet.createAttribute<std::string>("Description", "Center of gravity of the body in the absolute frame");

      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/BodyPosition",
                   static_cast<Eigen::Vector3d> (body->GetPosition()));
      bodyGroup.getDataSet("BodyPosition").createAttribute<std::string>("Description",
                                                                        "Center of gravity of the body in the absolute frame");

      bodyGroup.createGroup("Mask");
      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/Mask/ForceMask",
                   static_cast<Eigen::Matrix<bool, 6, 1>> (body->GetForceMask().GetMask()));
      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/Mask/MotionMask",
                   static_cast<Eigen::Matrix<bool, 6, 1>> (body->GetMotionMask().GetMask()));

      bodyGroup.createGroup("Mesh");
      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/Mesh/NbFaces",
                   body->GetMesh()->n_faces());
      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/Mesh/Faces",
                   body->GetMesh()->GetFaces());
      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/Mesh/NbVertices",
                   body->GetMesh()->n_vertices());
      H5Easy::dump(file, "Bodies/Body_" + std::to_string(i) + "/Mesh/Vertices",
                   body->GetMesh()->GetVertices());


      WriteMesh(file, "Bodies/Body_" + std::to_string(i) + "/Mesh", body);

      WriteExcitation(excitationType::Diffraction, file,
                      "Bodies/Body_" + std::to_string(i) + "/Excitation/Diffraction", body);
      WriteExcitation(excitationType::Froude_Krylov, file,
                      "Bodies/Body_" + std::to_string(i) + "/Excitation/FroudeKrylov", body);

      WriteRadiation(file, "Bodies/Body_" + std::to_string(i) + "/Radiation", body);

    }


  }

  void HydrodynamicDataBase::WriteExcitation(excitationType type, HighFive::File &HDF5_file,
                                             const std::string &path, Body *body) {

    for (unsigned int iwaveDir = 0; iwaveDir < m_waveDirectionDiscretization.GetNbSample(); ++iwaveDir) {

      auto anglePath = path + "/Angle_" + std::to_string(iwaveDir);

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
      H5Easy::dump(HDF5_file, anglePath + "/Angle", m_waveDirectionDiscretization.GetVector()[iwaveDir]);
      angleGroup.getDataSet("Angle").createAttribute<std::string>("Description", "Wave direction.");
      angleGroup.getDataSet("Angle").createAttribute<std::string>("Unit", "deg");

      H5Easy::dump(HDF5_file, anglePath + "/RealCoeffs", static_cast<Eigen::MatrixXd>(coeff.real()));
      angleGroup.getDataSet("RealCoeffs").createAttribute<std::string>("Description",
                                                                       "Real part of the Froude-Krylov loads on body " +
                                                                       std::to_string(body->GetID()) +
                                                                       " for a wave direction of " + std::to_string(
                                                                           m_waveDirectionDiscretization.GetVector()[iwaveDir]) +
                                                                       " deg.");
      angleGroup.getDataSet("RealCoeffs").createAttribute<std::string>("Unit", "N/m");

      H5Easy::dump(HDF5_file, anglePath + "/ImagCoeffs", static_cast<Eigen::MatrixXd>(coeff.imag()));
      angleGroup.getDataSet("ImagCoeffs").createAttribute<std::string>("Description",
                                                                       "Imaginary part of the Froude-Krylov loads on body " +
                                                                       std::to_string(body->GetID()) +
                                                                       " for a wave direction of " + std::to_string(
                                                                           m_waveDirectionDiscretization.GetVector()[iwaveDir]) +
                                                                       " deg.");
      angleGroup.getDataSet("ImagCoeffs").createAttribute<std::string>("Unit", "N/m");


    }

  }

  void HydrodynamicDataBase::WriteRadiation(HighFive::File &HDF5_file, const std::string &path, Body *body) {

    auto frequencies = GetFrequencyDiscretization().GetVector();

    for (unsigned int ibodyMotion = 0; ibodyMotion < m_nbody; ++ibodyMotion) {

      auto bodyMotion = this->GetBody(ibodyMotion);
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

      // Writing the impulse response functions.
//      auto impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionK", radiationMask);
//      body->SetImpulseResponseFunctionK(bodyMotion, impulseResponseFunctionsK);
//
//      impulseResponseFunctionsK = ReadComponents(HDF5_file, bodyMotionPath + "/ImpulseResponseFunctionKU", radiationMask);
//      body->SetImpulseResponseFunctionKu(bodyMotion, impulseResponseFunctionsK);


      // Writing the added mass matrix for the body.
      auto AddedMassGroup = bodyMotionGroup.createGroup("AddedMass");
      AddedMassGroup.createAttribute("Description", "Added mass coefficients for acceleration of body " +
                                                    std::to_string(bodyMotion->GetID()) +
                                                    " that radiates waves and  generate force on body " +
                                                    std::to_string(body->GetID()) + ".");
      for (unsigned int i = 0; i < 6; i++) {
        H5Easy::dump(HDF5_file, bodyMotionPath + "/AddedMass/DOF_" + std::to_string(i),
                     body->GetMatrixComponentFromIterator(body->GetAddedMassInterpolator(bodyMotion, i),
                                                          GetFrequencyDiscretization()));
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
                     body->GetMatrixComponentFromIterator(body->GetRadiationDampingInterpolator(bodyMotion, i),
                                                          GetFrequencyDiscretization()));
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
                     body->GetMatrixComponentFromIterator(body->GetIRFInterpolatorK(bodyMotion, i),
                                                          GetTimeDiscretization()));
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
                     body->GetMatrixComponentFromIterator(body->GetIRFInterpolatorKu(bodyMotion, i),
                                                          GetTimeDiscretization()));
        auto DOF = KUGroup.getDataSet("DOF_" + std::to_string(i));
        DOF.createAttribute("Description", "Impulse response functions KU");
      }


    }


  }


} // namespace HDB5_io
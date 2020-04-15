//
// Created by lletourn on 26/02/20.
//

#include <Eigen/Dense>

#include "HydrodynamicDataBase.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "WaveDrift.h"

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

  Body *HydrodynamicDataBase::NewBody(unsigned int id, const std::string &name) {
    m_bodies.push_back(std::make_shared<Body>(id, name, this));
    return (m_bodies.back()).get();
  }

  Body *HydrodynamicDataBase::GetBody(int id) const {
    return m_bodies[id].get();
  }

  void HydrodynamicDataBase::SetNbBodies(int nb) {
    m_nbody = nb;
  }

  int HydrodynamicDataBase::GetNbBodies() const {
    return m_nbody;
  }

  void HydrodynamicDataBase::Export_HDF5(const std::string &HDF5_file) {

    using namespace HighFive;

    File file(HDF5_file, File::ReadWrite | File::Create | File::Truncate);

    DataSet dataSet = file.createDataSet<std::string>("CreationDate", DataSpace::From(m_creationDate));
    dataSet.write(m_creationDate);
    dataSet.createAttribute<std::string>("Description", "Date of the creation of this database.");

    m_version = 3.;
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
    H5Easy::dump(file, "Discretizations/Frequency",
                 static_cast<Eigen::Matrix<double, Eigen::Dynamic, 1>> (GetFrequencyDiscretization()));
    H5Easy::dump(file, "Discretizations/WaveDirection",
                 static_cast<Eigen::Matrix<double, Eigen::Dynamic, 1>> (GetWaveDirectionDiscretization()));
    H5Easy::dump(file, "Discretizations/Time",
                 static_cast<Eigen::Matrix<double, Eigen::Dynamic, 1>> (GetTimeDiscretization()));

    if (m_waveDrift) {
      WriteWaveDrift(file);
    }

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

      WriteExcitation(excitationType::Diffraction, file,
                      "Bodies/Body_" + std::to_string(i) + "/Excitation/Diffraction", body);
      WriteExcitation(excitationType::Froude_Krylov, file,
                      "Bodies/Body_" + std::to_string(i) + "/Excitation/FroudeKrylov", body);

      WriteRadiation(file, "Bodies/Body_" + std::to_string(i) + "/Radiation", body);

      WriteRAO(file, "Bodies/Body_" + std::to_string(i) + "/RAO", body);

    }


  }

  void HydrodynamicDataBase::WriteExcitation(excitationType type, HighFive::File &HDF5_file,
                                             const std::string &path, Body *body) {

    for (unsigned int iwaveDir = 0; iwaveDir < GetWaveDirectionDiscretization().size(); ++iwaveDir) {

      auto anglePath = path + "/Angle_" + std::to_string(iwaveDir);

      auto angle = GetWaveDirectionDiscretization()(iwaveDir);

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

  void HydrodynamicDataBase::WriteRadiation(HighFive::File &HDF5_file, const std::string &path, Body *body) {

//    auto frequencies = GetFrequencyDiscretization().GetVector();

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
                     body->GetHDBInterpolatedData(Body::interpolatedData::ADDED_MASS, bodyMotion, i,
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
                     body->GetHDBInterpolatedData(Body::interpolatedData::RADIATION_DAMPING, bodyMotion, i,
                                                  GetFrequencyDiscretization()));
//                     body->GetMatrixComponentFromIterator(body->GetRadiationDampingInterpolator(bodyMotion, i),
//                                                          GetFrequencyDiscretization()));
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
                     body->GetHDBInterpolatedData(Body::interpolatedData::IRF_K, bodyMotion, i,
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
                     body->GetHDBInterpolatedData(Body::interpolatedData::IRF_KU, bodyMotion, i,
                                                  GetTimeDiscretization()));
        auto DOF = KUGroup.getDataSet("DOF_" + std::to_string(i));
        DOF.createAttribute("Description", "Impulse response functions KU");
      }


    }


  }

  void HydrodynamicDataBase::WriteRAO(HighFive::File &HDF5_file, const std::string &path, Body *body) {

    for (unsigned int iwaveDir = 0; iwaveDir < GetWaveDirectionDiscretization().size(); ++iwaveDir) {

      auto anglePath = path + "/Angle_" + std::to_string(iwaveDir);

      auto angle = GetWaveDirectionDiscretization()(iwaveDir);

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

  void HydrodynamicDataBase::WriteWaveDrift(HighFive::File &HDF5_file) {

    auto symmetries = m_waveDrift->GetSymmetries();

    unsigned int sym_x = static_cast<unsigned int>(symmetries[0]);
    unsigned int sym_y = static_cast<unsigned int>(symmetries[1]);

    auto waveDriftGroup = HDF5_file.createGroup("WaveDrift");

    HighFive::DataSet dataSet = waveDriftGroup.createDataSet<unsigned int>("sym_x", HighFive::DataSpace::From(sym_x));
    dataSet.write(sym_x);
    dataSet.createAttribute<std::string>("Description", "Symmetry along x");

    dataSet = waveDriftGroup.createDataSet<unsigned int>("sym_y", HighFive::DataSpace::From(sym_y));
    dataSet.write(sym_y);
    dataSet.createAttribute<std::string>("Description", "Symmetry along y");

    std::vector<std::string> components = {"surge", "sway", "yaw"};

    for (auto & component: components) {

      auto componentGroup = waveDriftGroup.createGroup(component);

      for (unsigned int iangle = 0; iangle < GetWaveDirectionDiscretization().size(); iangle++) {
        auto angle = GetWaveDirectionDiscretization()(iangle);
        auto headingGroup = componentGroup.createGroup("angle_" + std::to_string(iangle));
        Eigen::VectorXd data(GetFrequencyDiscretization().size());
        for (unsigned int i=0; i<GetFrequencyDiscretization().size(); i++) {
          data(i) = m_waveDrift->Eval(component, GetFrequencyDiscretization()(i), angle);
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


} // namespace HDB5_io
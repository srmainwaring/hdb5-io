//
// Created by lletourn on 26/02/20.
//

#include "Body.h"
#include "HydrodynamicDataBase.h"


namespace HDB5_io {

  Body::Body(unsigned int id, const std::string &name, HydrodynamicDataBase *hdb) : m_id(id), m_name(name), m_HDB(hdb) {
    m_mesh = std::make_shared<meshoui::Mesh>();
  }

  //
  // Setters
  //

  void Body::SetForceMask(mathutils::Vector6d<int> mask) {
    m_forceMask.SetMask(mask);
  }

  void Body::SetMotionMask(mathutils::Vector6d<int> mask) {
    m_motionMask.SetMask(mask);
  }

  void Body::SetDiffraction(unsigned int iangle, const Eigen::MatrixXcd &diffractionMatrix) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().GetNbSample());
    assert(diffractionMatrix.rows() == 6);
    assert(diffractionMatrix.cols() == m_HDB->GetFrequencyDiscretization().GetNbSample());
    m_diffraction[iangle] = diffractionMatrix;
  }

  void Body::SetFroudeKrylov(unsigned int iangle, const Eigen::MatrixXcd &froudeKrylovMatrix) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().GetNbSample());
    assert(froudeKrylovMatrix.rows() == 6);
    assert(froudeKrylovMatrix.cols() == m_HDB->GetFrequencyDiscretization().GetNbSample());
    m_froudeKrylov[iangle] = froudeKrylovMatrix;
  }

  void Body::SetExcitation(unsigned int iangle, const Eigen::MatrixXcd &excitationMatrix) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().GetNbSample());
    assert(excitationMatrix.rows() == 6);
    assert(excitationMatrix.cols() == m_HDB->GetFrequencyDiscretization().GetNbSample());
    m_excitation[iangle] = excitationMatrix;
  }

  void Body::ComputeExcitation() {

    /// This subroutine computes the excitation loads from the diffraction loads and the Froude-Krylov loads.

    assert(m_diffraction.size() == m_froudeKrylov.size());
//    for (unsigned int iangle = 0; iangle < GetNbWaveDirections(); ++iangle) {
    for (unsigned int iangle = 0; iangle < m_diffraction.size(); ++iangle) {
      m_excitation[iangle] = m_diffraction[iangle] + m_froudeKrylov[iangle];
    }
  }

  void Body::SetInfiniteAddedMass(Body *BodyMotion, const mathutils::Matrix66<double> &CMInf) {
    m_infiniteAddedMass[BodyMotion] = CMInf;
  }

//  void Body::SetRadiationMask(HDB5_io::Body *BodyMotion, const mathutils::Matrix66<bool> &mask){
//    m_radiationMask[BodyMotion] = mask;
//  }

  void Body::SetRadiationMask(HDB5_io::Body *BodyMotion, const mathutils::Matrix66<int> &mask) {

    assert(mask.maxCoeff() <= 1 and mask.minCoeff() >= 0);

    for (unsigned int i = 0; i < 6; i++) {
      for (unsigned int j = 0; j < 6; j++) {
        m_radiationMask[BodyMotion](i, j) = mask(i, j) == 1;
      }
    }
  }

  void Body::SetImpulseResponseFunctionK(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listIRF) {
    for (auto &IRF: listIRF) {
      assert(IRF.rows() == 6);
      assert(IRF.cols() == m_HDB->GetTimeDiscretization().GetNbSample());

      auto vtime = std::make_shared<std::vector<double>>(m_HDB->GetTimeDiscretization().GetVector());

      auto vdata = std::make_shared<std::vector<mathutils::Vector6d<double>>>();
      for (unsigned int j = 0; j < IRF.cols(); ++j) {
        vdata->push_back(IRF.col(j));
      }

      auto interp = std::make_shared<mathutils::Interp1dLinear<double, mathutils::Vector6d<double>>>();
      interp->Initialize(vtime, vdata);

      m_interpK[BodyMotion].push_back(interp);
    }

  }

  void Body::SetImpulseResponseFunctionKu(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listIRF) {
    for (auto &IRF: listIRF) {
      assert(IRF.rows() == 6);
      assert(IRF.cols() == m_HDB->GetTimeDiscretization().GetNbSample());

      auto vtime = std::make_shared<std::vector<double>>(m_HDB->GetTimeDiscretization().GetVector());

      auto vdata = std::make_shared<std::vector<mathutils::Vector6d<double>>>();
      for (unsigned int j = 0; j < IRF.cols(); ++j) {
        vdata->push_back(IRF.col(j));
      }

      auto interp = std::make_shared<mathutils::Interp1dLinear<double, mathutils::Vector6d<double>>>();
      interp->Initialize(vtime, vdata);

      m_interpKu[BodyMotion].push_back(interp);
    }

  }

  void Body::SetAddedMass(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listAddedMass) {
    for (auto &addedMass: listAddedMass) {
      assert(addedMass.rows() == 6);
      assert(addedMass.cols() == m_HDB->GetFrequencyDiscretization().GetNbSample());

      auto frequencies = std::make_shared<std::vector<double>>(m_HDB->GetFrequencyDiscretization().GetVector());

      auto vdata = std::make_shared<std::vector<mathutils::Vector6d<double>>>();
      for (unsigned int j = 0; j < addedMass.cols(); ++j) {
        vdata->push_back(addedMass.col(j));
      }

      auto interp = std::make_shared<mathutils::Interp1dLinear<double, mathutils::Vector6d<double>>>();
      interp->Initialize(frequencies, vdata);

      m_addedMass[BodyMotion].push_back(interp);
    }

  }

  void Body::SetRadiationDamping(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listRadiationDamping) {
    for (auto &radiationDamping: listRadiationDamping) {
      assert(radiationDamping.rows() == 6);
      assert(radiationDamping.cols() == m_HDB->GetFrequencyDiscretization().GetNbSample());

      auto frequencies = std::make_shared<std::vector<double>>(m_HDB->GetFrequencyDiscretization().GetVector());

      auto vdata = std::make_shared<std::vector<mathutils::Vector6d<double>>>();
      for (unsigned int j = 0; j < radiationDamping.cols(); ++j) {
        vdata->push_back(radiationDamping.col(j));
      }

      auto interp = std::make_shared<mathutils::Interp1dLinear<double, mathutils::Vector6d<double>>>();
      interp->Initialize(frequencies, vdata);

      m_radiationDamping[BodyMotion].push_back(interp);
    }

  }

  void Body::SetStiffnessMatrix(const mathutils::Matrix33<double> &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix;
  }

  void Body::SetStiffnessMatrix(const mathutils::Matrix66<double> &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix.block<3, 3>(2, 2);
  }

  void Body::LoadMesh(const std::vector<mathutils::Vector3d<double>> &vertices, const std::vector<Eigen::VectorXi>& faces) {
    m_mesh->Load(vertices, faces);
  }

  //
  // Getters
  //

  Eigen::MatrixXcd Body::GetDiffraction(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_diffraction[iangle];
  }

  Eigen::VectorXcd Body::GetDiffraction(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
//    assert(iforce < 6);
    return m_diffraction[iangle].row(iforce);
  }

  Eigen::MatrixXcd Body::GetFroudeKrylov(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_froudeKrylov[iangle];
  }

  Eigen::VectorXcd Body::GetFroudeKrylov(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
//    assert(iforce < 6);
    return m_froudeKrylov[iangle].row(iforce);
  }

  Eigen::MatrixXcd Body::GetExcitation(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_excitation[iangle];
  }

  Eigen::VectorXcd Body::GetExcitation(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
//    assert(iforce < 6);
    return m_excitation[iangle].row(iforce);
  }

  mathutils::Matrix66<double> Body::GetInfiniteAddedMass(Body *BodyMotion) const {
    return m_infiniteAddedMass.at(BodyMotion);
  }

  mathutils::Matrix66<bool> Body::GetRadiationMask(Body *BodyMotion) const {
    return m_radiationMask.at(BodyMotion);
  }

  mathutils::Matrix66<double> Body::GetSelfInfiniteAddedMass() {
    return m_infiniteAddedMass[this];
  }

  mathutils::Interp1d<double, mathutils::Vector6d<double>> *
  Body::GetIRFInterpolatorK(Body *BodyMotion, unsigned int idof) {
    assert(idof < 6);
    return m_interpK[BodyMotion][idof].get();
  };

  mathutils::Interp1d<double, mathutils::Vector6d<double>> *
  Body::GetIRFInterpolatorKu(Body *BodyMotion, unsigned int idof) {
    assert(idof < 6);
    return m_interpKu[BodyMotion][idof].get();
  };

  mathutils::Interp1d<double, mathutils::Vector6d<double>> *
  Body::GetAddedMassInterpolator(Body *BodyMotion, unsigned int idof) {
    assert(idof < 6);
    return m_addedMass[BodyMotion][idof].get();
  };

  mathutils::Interp1d<double, mathutils::Vector6d<double>> *
  Body::GetRadiationDampingInterpolator(Body *BodyMotion, unsigned int idof) {
    assert(idof < 6);
    return m_radiationDamping[BodyMotion][idof].get();
  };

  Eigen::MatrixXd Body::GetMatrixComponentFromIterator(mathutils::Interp1d<double, mathutils::Vector6d<double>>* interpolator, Discretization1D frequencies) {
    Eigen::MatrixXd addedMasses(6, frequencies.GetNbSample());
    for (unsigned int i=0; i<frequencies.GetNbSample(); i++) {
      addedMasses.col(i) = interpolator->Eval(frequencies.GetVector()[i]);
    }
    return addedMasses;
  }


  void Body::AllocateAll(unsigned int nFrequencies, unsigned int nDirections) {

    // This subroutine allocates the arrays for the hdb.

    auto nbForce = GetForceMask().GetNbDOF();

    // --> Allocating arrays for excitations

//    auto nbWaveDir = GetNbWaveDirections();
    m_diffraction.reserve((unsigned long) nDirections);
    m_froudeKrylov.reserve((unsigned long) nDirections);
    m_excitation.reserve((unsigned long) nDirections);

//    auto nbFreq = GetNbFrequencies();
    for (int i = 0; i < nDirections; ++i) {
      Eigen::MatrixXcd mat(nbForce, nFrequencies);
      m_diffraction.push_back(mat);
      m_froudeKrylov.push_back(mat);
      m_excitation.push_back(mat);
    }

    /// --> Allocating arrays for radiation
    m_radiationMask.reserve(m_HDB->GetNbBodies());


  }

  void Body::BoxMesh() {
    // Build a box
    std::vector<meshoui::Vector3d> vertices;
    vertices.emplace_back(-1, -1,  1);
    vertices.emplace_back(1, -1,  1);
    vertices.emplace_back(1,  1,  1);
    vertices.emplace_back(-1,  1,  1);
    vertices.emplace_back(-1, -1, -1);
    vertices.emplace_back(1, -1, -1);
    vertices.emplace_back(1,  1, -1);
    vertices.emplace_back(-1,  1, -1);
    std::vector<Eigen::VectorXi> faces;
    faces.emplace_back(Eigen::Vector3i(0,1,2));
    faces.emplace_back(Eigen::Vector3i(2,3,0));
    faces.emplace_back(Eigen::Vector3i(0,4,1));
    faces.emplace_back(Eigen::Vector3i(1,4,5));
    faces.emplace_back(Eigen::Vector3i(1,5,2));
    faces.emplace_back(Eigen::Vector3i(2,5,6));
    faces.emplace_back(Eigen::Vector3i(2,6,3));
    faces.emplace_back(Eigen::Vector3i(3,6,7));
    faces.emplace_back(Eigen::Vector3i(3,7,0));
    faces.emplace_back(Eigen::Vector3i(0,7,4));
    faces.emplace_back(Eigen::Vector3i(6,5,4));
    faces.emplace_back(Eigen::Vector3i(7,6,4));
    // Load it
    m_mesh->Load(vertices, faces);
  }

#ifdef H5_USE_VTK
  void Body::VisualizeMesh() const {
    // Creation of a VTK mesh.
    meshoui::VTKMesh vtkmesh = meshoui::VTKMesh(*m_mesh);

    // Visualization.
    vtkmesh.Visualize();
  }
#endif

}
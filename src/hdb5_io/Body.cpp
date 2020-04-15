//
// Created by lletourn on 26/02/20.
//

#include "Body.h"

#ifdef H5_USE_VTK

#include "meshoui/vtkmesh.h"

#endif

#include "HydrodynamicDataBase.h"


namespace HDB5_io {

  Body::Body(unsigned int id, const std::string &name, HydrodynamicDataBase *hdb) : m_id(id), m_name(name), m_HDB(hdb) {
    m_mesh = std::make_shared<Mesh>();
    m_interpK = std::make_shared<HDBinterpolator>();
    m_interpKu = std::make_shared<HDBinterpolator>();
    m_addedMass = std::make_shared<HDBinterpolator>();
    m_radiationDamping = std::make_shared<HDBinterpolator>();
    AllocateAll();
  }

  //
  // Setters
  //

  void Body::SetPosition(const mathutils::Vector3d<double> &position) {
    m_position = position;
  }

  void Body::SetForceMask(mathutils::Vector6d<int> mask) {
    m_forceMask.SetMask(mask);
  }

  void Body::SetMotionMask(mathutils::Vector6d<int> mask) {
    m_motionMask.SetMask(mask);
  }

  void Body::SetDiffraction(unsigned int iangle, const Eigen::MatrixXcd &diffractionMatrix) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().size());
    assert(diffractionMatrix.rows() == 6);
    assert(diffractionMatrix.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_diffraction[iangle] = diffractionMatrix;
  }

  void Body::SetFroudeKrylov(unsigned int iangle, const Eigen::MatrixXcd &froudeKrylovMatrix) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().size());
    assert(froudeKrylovMatrix.rows() == 6);
    assert(froudeKrylovMatrix.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_froudeKrylov[iangle] = froudeKrylovMatrix;
  }

  void Body::SetExcitation(unsigned int iangle, const Eigen::MatrixXcd &excitationMatrix) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().size());
    assert(excitationMatrix.rows() == 6);
    assert(excitationMatrix.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_excitation[iangle] = excitationMatrix;
  }

  void Body::ComputeExcitation() {

    /// This subroutine computes the excitation loads from the diffraction loads and the Froude-Krylov loads.

    assert(m_diffraction.size() == m_froudeKrylov.size());
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

  void Body::SetRAO(unsigned int iangle, const Eigen::MatrixXcd &RAO) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().size());
    assert(RAO.rows() == 6);
    assert(RAO.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_RAO[iangle] = RAO;
    m_isRAO = true;
  }

  void Body::SetHDBInterpolator(interpolatedData type, Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listData) {

    unsigned int idof = 0;

    if (GetHDBInterpolator(type)->count(idof) == 0) {

      for (auto &data: listData) {
        assert(data.rows() == 6);

        auto table = std::make_shared<mathutils::LookupTable1D<double, mathutils::Vector6d<double>>>();

        switch (type) {
          case IRF_K: {
            auto EigenTime = m_HDB->GetTimeDiscretization();
            assert(data.cols() == EigenTime.size());
            std::vector<double> time(&EigenTime(0,0), EigenTime.data()+EigenTime.size());
            table->SetX(time);
            break;
          }
          case IRF_KU: {
            auto EigenTime = m_HDB->GetTimeDiscretization();
            assert(data.cols() == EigenTime.size());
            std::vector<double> time(&EigenTime(0,0), EigenTime.data()+EigenTime.size());
            table->SetX(time);
            break;
          }
          case ADDED_MASS: {
            assert(data.cols() == m_HDB->GetFrequencyDiscretization().size());
            auto EigenFreq = m_HDB->GetFrequencyDiscretization();
            std::vector<double> freq(&EigenFreq(0, 0), EigenFreq.data() + EigenFreq.size());
            table->SetX(freq);
            break;
          }
          case RADIATION_DAMPING: {
            assert(data.cols() == m_HDB->GetFrequencyDiscretization().size());
            auto EigenFreq = m_HDB->GetFrequencyDiscretization();
            std::vector<double> freq(&EigenFreq(0, 0), EigenFreq.data() + EigenFreq.size());
            table->SetX(freq);
            break;
          }
        }

        std::vector<mathutils::Vector6d<double>> vdata;
        for (unsigned int j = 0; j < data.cols(); ++j) {
          vdata.emplace_back(data.col(j));
        }
        table->AddY(BodyMotion->GetName(), vdata);

        switch (type) {
          case IRF_K: {
            m_interpK->insert(std::make_pair(idof, table));
            break;
          }
          case IRF_KU: {
            m_interpKu->insert(std::make_pair(idof, table));
            break;
          }
          case ADDED_MASS: {
            m_addedMass->insert(std::make_pair(idof, table));
            break;
          }
          case RADIATION_DAMPING: {
            m_radiationDamping->insert(std::make_pair(idof, table));
            break;
          }
        }

        idof++;
      }
    } else {
      for (auto &data: listData) {
        assert(data.rows() == 6);
        std::vector<mathutils::Vector6d<double>> vdata;
        for (unsigned int j = 0; j < data.cols(); ++j) {
          vdata.emplace_back(data.col(j));
        }
        GetHDBInterpolator(type)->at(idof)->AddY(BodyMotion->GetName(), vdata);
        idof++;
      }
    }

  }

  void Body::SetStiffnessMatrix(const mathutils::Matrix33<double> &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix;
  }

  void Body::SetStiffnessMatrix(const mathutils::Matrix66<double> &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix.block<3, 3>(2, 2);
  }

  void
  Body::LoadMesh(const std::vector<mathutils::Vector3d<double>> &vertices, const std::vector<Eigen::VectorXi> &faces) {
    m_mesh->Load(vertices, faces);
  }

  //
  // Getters
  //

  mathutils::Vector3d<double> Body::GetPosition() const {
    return m_position;
  }

  Mask Body::GetMotionMask() const {
    return m_motionMask;
  }

  Mask Body::GetForceMask() const {
    return m_forceMask;
  }

  Mesh *Body::GetMesh() const {
    return m_mesh.get();
  }

  mathutils::Matrix33<double> Body::GetHydrostaticStiffnessMatrix() const {
    return m_hydrostaticStiffnessMatrix;
  }

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

  Eigen::MatrixXcd Body::GetRAO(const unsigned int iangle) const {
//    assert(iangle < this->GetNbWaveDirections());
    return m_RAO[iangle];
  }

  Eigen::VectorXcd Body::GetRAO(const unsigned int iangle, const unsigned iforce) const {
//    assert(iangle < this->GetNbWaveDirections());
    assert(iforce < 6);
    return m_RAO[iangle].row(iforce);
  }

  Body::HDBinterpolator *Body::GetHDBInterpolator(interpolatedData type) {
    switch (type) {
      case IRF_K:
        return m_interpK.get();
      case IRF_KU:
        return m_interpKu.get();
      case ADDED_MASS:
        return m_addedMass.get();
      case RADIATION_DAMPING:
        return m_radiationDamping.get();
    }
  }

  Eigen::MatrixXd
  Body::GetHDBInterpolatedData(interpolatedData type, Body *BodyMotion, unsigned int idof,
                               mathutils::VectorN<double> frequencies) {

    Eigen::MatrixXd data(6, frequencies.size());
    for (unsigned int i = 0; i < frequencies.size(); i++) {
      data.col(i) = GetHDBInterpolator(type)->at(idof)->Eval(BodyMotion->GetName(), frequencies(i));
    }
    return data;

  }


  void Body::AllocateAll() {

    // This subroutine allocates the arrays for the hdb.

    auto nDirections = m_HDB->GetWaveDirectionDiscretization().size();
    auto nFrequencies = m_HDB->GetFrequencyDiscretization().size();
    m_diffraction.reserve((unsigned long) nDirections);
    m_froudeKrylov.reserve((unsigned long) nDirections);
    m_excitation.reserve((unsigned long) nDirections);
    m_RAO.reserve((unsigned long) nDirections);

//    auto nbFreq = GetNbFrequencies();
    for (int i = 0; i < nDirections; ++i) {
      Eigen::MatrixXcd mat(6, nFrequencies);
      m_diffraction.push_back(mat);
      m_froudeKrylov.push_back(mat);
      m_excitation.push_back(mat);
      m_RAO.push_back(mat);
    }

    /// --> Allocating arrays for radiation
    m_radiationMask.reserve(m_HDB->GetNbBodies());

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
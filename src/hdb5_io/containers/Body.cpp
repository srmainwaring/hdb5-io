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
    AllocateAll();
  }

  //
  // Setters
  //

  void Body::SetPosition(const mathutils::Vector3d<double> &position) {
    m_position = position;
  }

  void Body::SetForceMask(mathutils::Vector6d<bool> mask) {
    m_forceMask.SetMask(mask);
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

    // This method computes the excitation loads from the diffraction loads and the Froude-Krylov loads.

    assert(m_diffraction.size() == m_froudeKrylov.size());
    for (unsigned int iangle = 0; iangle < m_diffraction.size(); ++iangle) {
      m_excitation[iangle] = m_diffraction[iangle] + m_froudeKrylov[iangle];
    }
  }

  void Body::SetInfiniteAddedMass(Body *BodyMotion, const Matrix66 &CMInf) {
    m_infiniteAddedMass[BodyMotion] = CMInf;
  }

  void Body::SetZeroFreqAddedMass(Body *BodyMotion, const Matrix66 &CMZero) {
    m_zeroFreqAddedMass[BodyMotion] = CMZero;
  }

//  void Body::SetRadiationMask(HDB5_io::Body *BodyMotion, const Matrix66b &mask){
//    m_radiationMask[BodyMotion] = mask;
//  }

  void Body::SetRadiationMask(HDB5_io::Body *BodyMotion, const Matrix66b &mask) {
    assert(mask.maxCoeff() <= 1 and mask.minCoeff() >= 0);
    m_radiationMask[BodyMotion] = mask;
  }

  void Body::SetRAO(unsigned int iangle, const Eigen::MatrixXcd &RAO) {
    assert(iangle < m_HDB->GetWaveDirectionDiscretization().size());
    assert(RAO.rows() == 6);
    assert(RAO.cols() == m_HDB->GetFrequencyDiscretization().size());
    m_RAO[iangle] = RAO;
    m_isRAO = true;
  }

  void Body::SetAddedMass(Body *BodyMotion, const std::vector<Matrix66> &listData) {
    assert(listData.size() == m_HDB->GetFrequencyDiscretization().size());
    m_addedMass.insert(std::make_pair(BodyMotion,listData));
  }

  void Body::SetRadiationDamping(Body *BodyMotion, const std::vector<Matrix66> &listData) {
    assert(listData.size() == m_HDB->GetFrequencyDiscretization().size());
    m_radiationDamping.insert(std::make_pair(BodyMotion,listData));
  }

  void Body::AddAddedMass(Body *BodyMotion, const Matrix66 &Data) {
    if (m_addedMass.count(BodyMotion) == 0) {
      m_addedMass.at(BodyMotion).reserve(m_HDB->GetFrequencyDiscretization().size());
    }
    m_addedMass.at(BodyMotion).push_back(Data);
  }

  void Body::AddRadiationDamping(Body *BodyMotion, const Matrix66 &Data) {
    if (m_radiationDamping.count(BodyMotion) == 0) {
      m_radiationDamping.at(BodyMotion).reserve(m_HDB->GetFrequencyDiscretization().size());
    }
    m_radiationDamping.at(BodyMotion).push_back(Data);
  }

  void Body::SetIRF(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listData) {

    unsigned int idof = 0;

    if (m_interpK->count(idof) == 0) {

      for (auto &data: listData) {
        assert(data.rows() == 6);

        auto table = std::make_shared<mathutils::LookupTable1D<double, mathutils::Vector6d<double>>>();

        auto EigenTime = m_HDB->GetTimeDiscretization();
        assert(data.cols() == EigenTime.size());
        std::vector<double> time(&EigenTime(0, 0), EigenTime.data() + EigenTime.size());
        table->SetX(time);

        std::vector<mathutils::Vector6d<double>> vdata;
        for (unsigned int j = 0; j < data.cols(); ++j) {
          vdata.emplace_back(data.col(j));
        }
        table->AddY(BodyMotion->GetName(), vdata);

        m_interpK->insert(std::make_pair(idof, table));

        idof++;
      }
    } else {
      for (auto &data: listData) {
        assert(data.rows() == 6);
        std::vector<mathutils::Vector6d<double>> vdata;
        for (unsigned int j = 0; j < data.cols(); ++j) {
          vdata.emplace_back(data.col(j));
        }
        m_interpK->at(idof)->AddY(BodyMotion->GetName(), vdata);
        idof++;
      }
    }

    m_isIRF = true;

  }

  void Body::SetIRF_Ku(Body *BodyMotion, const std::vector<Eigen::MatrixXd> &listData) {

    unsigned int idof = 0;

    if (m_interpKu->count(idof) == 0) {

      for (auto &data: listData) {
        assert(data.rows() == 6);

        auto table = std::make_shared<mathutils::LookupTable1D<double, mathutils::Vector6d<double>>>();

        auto EigenTime = m_HDB->GetTimeDiscretization();
        assert(data.cols() == EigenTime.size());
        std::vector<double> time(&EigenTime(0, 0), EigenTime.data() + EigenTime.size());
        table->SetX(time);

        std::vector<mathutils::Vector6d<double>> vdata;
        for (unsigned int j = 0; j < data.cols(); ++j) {
          vdata.emplace_back(data.col(j));
        }
        table->AddY(BodyMotion->GetName(), vdata);

        m_interpKu->insert(std::make_pair(idof, table));

        idof++;
      }
    } else {
      for (auto &data: listData) {
        assert(data.rows() == 6);
        std::vector<mathutils::Vector6d<double>> vdata;
        for (unsigned int j = 0; j < data.cols(); ++j) {
          vdata.emplace_back(data.col(j));
        }
        m_interpKu->at(idof)->AddY(BodyMotion->GetName(), vdata);
        idof++;
      }
    }

    m_isIRF = true;

  }

  void Body::SetStiffnessMatrix(const Matrix33 &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix;
    m_isHydrostatic = true;
  }

  void Body::SetStiffnessMatrix(const Matrix66 &hydrostaticStiffnessMatrix) {
    m_hydrostaticStiffnessMatrix = hydrostaticStiffnessMatrix.block<3, 3>(2, 2);
    m_isHydrostatic = true;
  }

  void Body::SetInertia(const Matrix66 &inertiaMatrix) {
    m_inertia = inertiaMatrix;
    m_isInertia = true;
  }

  void Body::SetMooring(const Matrix66 &mooringMatrix) {
    m_mooringStiffnessMatrix = mooringMatrix;
    m_isMooring = true;
  }

  void Body::SetLinearDamping(const Matrix66 &linearDampingMatrix) {
    m_linearDampingMatrix = linearDampingMatrix;
    m_isDamping = true;
  }

  void
  Body::LoadMesh(const std::vector<mathutils::Vector3d<double>> &vertices, const std::vector<Eigen::VectorXi> &faces) {
    m_mesh->Load(vertices, faces);
  }

  void Body::AddModalCoefficients(Body *BodyMotion, std::vector<PoleResidue> modalCoefficients) {
    if (m_modalCoefficients.count(BodyMotion)>0) {
      auto coeff = m_modalCoefficients.find(BodyMotion);
      coeff->second.emplace_back(modalCoefficients);
    }else {
      std::vector<std::vector<PoleResidue>> coefficients;
      coefficients.emplace_back(modalCoefficients);
      m_modalCoefficients[BodyMotion] = coefficients;
//      auto paired = std::make_pair(BodyMotion, coefficients);
//      m_modalCoefficients.insert(paired);
    }

  }

  //
  // Getters
  //

  bool Body::HasRAO() const {
    return m_isRAO;
  }

  bool Body::HasModal(Body *BodyMotion) const {
    return m_modalCoefficients.count(BodyMotion) > 0;
  }

  bool Body::HasIRF() const {
    return m_isIRF;
  }

  bool Body::HasInertia() const {
    return m_isInertia;
  }

  bool Body::HasHydrostatic() const {
    return m_isHydrostatic;
  }

  bool Body::HasMooring() const {
    return m_isMooring;
  }

  bool Body::HasDamping() const {
    return m_isDamping;
  }

  bool Body::HasZeroFreqAddedMass(Body *BodyMotion) const {
    return m_zeroFreqAddedMass.count(BodyMotion) > 0;
  }

  mathutils::Vector3d<double> Body::GetPosition() const {
    return m_position;
  }

  Mask Body::GetForceMask() const {
    return m_forceMask;
  }

  Mesh *Body::GetMesh() const {
    return m_mesh.get();
  }

  Matrix33 Body::GetHydrostaticStiffnessMatrix() const {
    return m_hydrostaticStiffnessMatrix;
  }

  Matrix66 Body::GetInertiaMatrix() const {
    return m_inertia;
  }

  Matrix66 Body::GetMooringMatrix() const {
    return m_mooringStiffnessMatrix;
  }

  Matrix66 Body::GetDampingMatrix() const {
    return m_linearDampingMatrix;
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

  Matrix66 Body::GetInfiniteAddedMass(Body *BodyMotion) const {
    return m_infiniteAddedMass.at(BodyMotion);
  }

  Matrix66 Body::GetZeroFreqAddedMass(HDB5_io::Body *BodyMotion) const {
    return m_zeroFreqAddedMass.at(BodyMotion);
  }

  Matrix66b Body::GetRadiationMask(Body *BodyMotion) const {
    return m_radiationMask.at(BodyMotion);
  }

  Matrix66 Body::GetSelfInfiniteAddedMass() {
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

  Eigen::MatrixXd
  Body::GetIRFInterpolatedData(Body *BodyMotion, unsigned int idof,
                               mathutils::VectorN<double> frequencies) {

    Eigen::MatrixXd data(6, frequencies.size());
    for (unsigned int i = 0; i < frequencies.size(); i++) {
      data.col(i) = m_interpK->at(idof)->Eval(BodyMotion->GetName(), frequencies(i));
    }
    return data;

  }

  Eigen::MatrixXd
  Body::GetIRF_KuInterpolatedData(Body *BodyMotion, unsigned int idof,
                               mathutils::VectorN<double> frequencies) {

    Eigen::MatrixXd data(6, frequencies.size());
    for (unsigned int i = 0; i < frequencies.size(); i++) {
      data.col(i) = m_interpKu->at(idof)->Eval(BodyMotion->GetName(), frequencies(i));
    }
    return data;

  }

  std::vector<PoleResidue> Body::GetModalCoefficients(Body *BodyMotion, int idof) {
    return m_modalCoefficients.at(BodyMotion)[idof];
  }

  PoleResidue Body::GetModalCoefficients(Body *BodyMotion, int idof, int iforce) {
    return m_modalCoefficients.at(BodyMotion)[idof][iforce];
  }


  void Body::AllocateAll() {

    // This subroutine allocates the arrays for the hdb.

    auto nDirections = m_HDB->GetWaveDirectionDiscretization().size();
    auto nFrequencies = m_HDB->GetFrequencyDiscretization().size();
    m_diffraction.reserve((unsigned long) nDirections);
    m_froudeKrylov.reserve((unsigned long) nDirections);
    m_excitation.reserve((unsigned long) nDirections);
    m_RAO.reserve((unsigned long) nDirections);

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

  Body::HDBinterpolator *Body::GetIRFInterpolator() const {
    return m_interpK.get();
  }

  Body::HDBinterpolator *Body::GetIRF_KuInterpolator() const {
    return m_interpKu.get();
  }

  Eigen::MatrixXd Body::GetAddedMass(Body *BodyMotion, unsigned int iforce) {
    Eigen::MatrixXd AM = Eigen::MatrixXd::Zero(6,m_HDB->GetFrequencyDiscretization().size());

    for (int imotion = 0; imotion < 6; ++imotion) {
      for (int iw = 0; iw < m_HDB->GetFrequencyDiscretization().size(); ++iw) {
        AM(imotion, iw) = m_addedMass.at(BodyMotion)[iw](iforce, imotion);
      }
    }

    return AM;
  }

  Eigen::MatrixXd Body::GetRadiationDamping(Body *BodyMotion, unsigned int imotion) {
    Eigen::MatrixXd RD = Eigen::MatrixXd::Zero(6, m_HDB->GetFrequencyDiscretization().size());

    for (int iforce = 0; iforce < 6; ++iforce) {
      for (int iw = 0; iw < m_HDB->GetFrequencyDiscretization().size(); ++iw) {
        RD(iforce, iw) = m_radiationDamping.at(BodyMotion)[iw](iforce, imotion);
      }
    }

    return RD;
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
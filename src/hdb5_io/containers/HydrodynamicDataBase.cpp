//
// Created by lletourn on 26/02/20.
//

#include "HydrodynamicDataBase.h"

#include "WaveDrift.h"

#include "Body.h"

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

  void HydrodynamicDataBase::SetFrequencyDiscretization(const mathutils::VectorN<double> &frequency) {
    m_frequencyDiscretization = frequency;
  }

  void HydrodynamicDataBase::SetWaveDirectionDiscretization(const mathutils::VectorN<double> &directions) {
    m_waveDirectionDiscretization = directions;
  }

  void HydrodynamicDataBase::SetTimeDiscretization(const mathutils::VectorN<double> &time) {
    m_timeDiscretization = time;
  }

  mathutils::VectorN<double> HydrodynamicDataBase::GetFrequencyDiscretization() const {
    return m_frequencyDiscretization;
  }

  mathutils::VectorN<double> HydrodynamicDataBase::GetWaveDirectionDiscretization() const {
    return m_waveDirectionDiscretization;
  }

  mathutils::VectorN<double> HydrodynamicDataBase::GetTimeDiscretization() const {
    return m_timeDiscretization;
  }

  void HydrodynamicDataBase::SetWaveDrift(const std::shared_ptr<WaveDrift> &wavedrift) {
    m_waveDrift = wavedrift;
  }

  std::shared_ptr<WaveDrift> HydrodynamicDataBase::GetWaveDrift() const {
    return m_waveDrift;
  }

  void HydrodynamicDataBase::SetVF() {
    m_isVF = true;
  }

  bool HydrodynamicDataBase::GetVF() const {
    return m_isVF;
  }

  void HydrodynamicDataBase::SetVFRelaxed(const int &relaxed) {
    m_VF_relaxed = relaxed;
  }

  int HydrodynamicDataBase::GetVFRelaxed() const {
    return m_VF_relaxed;
  }

  void HydrodynamicDataBase::SetVFMaxOrder(const int &order) {
    m_VF_max_order = order;
  }

  int HydrodynamicDataBase::GetVFMaxOrder() const {
    return m_VF_max_order;
  }

  void HydrodynamicDataBase::SetVFTolerance(const double &tolerance) {
    m_VF_tolerance = tolerance;
  }

  double HydrodynamicDataBase::GetVFTolerance() const {
    return m_VF_tolerance;
  }

  void HydrodynamicDataBase::SetWaveField() {
    m_isWaveField = true;
  }

  bool HydrodynamicDataBase::GetWaveField() const {
    return m_isSymmetries;
  }

  void HydrodynamicDataBase::SetSymmetries() {
    m_isSymmetries = true;
  }

  bool HydrodynamicDataBase::GetSymmetries() const {
    return m_isWaveField;
  }

  void HydrodynamicDataBase::SetSymBottom(const bool &sym_bottom) {
    m_sym_bottom = sym_bottom;
  }

  bool HydrodynamicDataBase::GetSymBottom() const {
    return m_sym_bottom;
  }

  void HydrodynamicDataBase::SetSymXOZ(const bool &sym_xoz) {
    m_sym_xoz = sym_xoz;
  }

  bool HydrodynamicDataBase::GetSymXOZ() const {
    return m_sym_xoz;
  }

  void HydrodynamicDataBase::SetSymYOZ(const bool &sym_yoz) {
    m_sym_yoz = sym_yoz;
  }

  bool HydrodynamicDataBase::GetSymYOZ() const {
    return m_sym_yoz;
  }

} // namespace HDB5_io
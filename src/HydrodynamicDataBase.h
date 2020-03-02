//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_HYDRODYNAMICDATABASE_H
#define HDB5_IO_HYDRODYNAMICDATABASE_H

#include <iostream>
#include <memory>
#include "Discretization1D.h"
#include "Body.h"

#include "MathUtils/Vector3d.h"
#include "highfive/H5Group.hpp"

namespace HDB5_io {

  // Forward declarations


  class HydrodynamicDataBase {

   public:

    HydrodynamicDataBase();

    void Import_HDF5(const std::string &HDF5_file);

    void Export_HDF5(const std::string &HDF5_file);


    // Accessors

    void SetCreationDate(std::string date) { m_creationDate = date; }

    std::string GetCreationDate() const { return m_creationDate; }

    void SetSolver(std::string solver) { m_solver = solver; }

    std::string GetSolver() const { return m_solver; }

    void SetVersion(double version) { m_version = version; }

    double GetVersion() const { return m_version; }

    void SetGravityAcceleration(double g) { m_gravityAcceleration = g; }

    double GetGravityAcceleration() const { return m_gravityAcceleration; }

    void SetWaterDensity(double rho) { m_waterDensity = rho; }

    double GetWaterDensity() const { return m_waterDensity; }

    void SetWaterDepth(double h) { m_waterDepth = h; }

    double GetWaterDepth() const { return m_waterDepth; }

    void SetNormalizationLength(double L) { m_normalizationLength = L; }

    double GetNormalizationLength() const { return m_normalizationLength; }

    Body *GetBody(int id) const { return m_bodies[id].get(); }

    Body *NewBody(unsigned int id, const std::string& name) {
      m_bodies.push_back(std::make_unique<Body>(id, name));
      return (m_bodies.back()).get();
    }

    void SetFrequencyDiscretization(double wmin, double wmax, unsigned int nw) { m_frequencyDiscretization = Discretization1D(wmin, wmax, nw);}

    void SetWaveDirectionDiscretization(double tmin, double tmax, unsigned int nt) { m_waveDirectionDiscretization = Discretization1D(tmin, tmax, nt);}

    void SetTimeDiscretization(double tmin, double tmax, unsigned int nt) { m_timeDiscretization = Discretization1D(tmin, tmax, nt);}

    Discretization1D GetFrequencyDiscretization() const { return m_frequencyDiscretization; }

    Discretization1D GetWaveDirectionDiscretization() const { return m_waveDirectionDiscretization; }

    Discretization1D GetTimeDiscretization() const { return m_timeDiscretization; }


   private:

    std::string m_creationDate;       ///< Creation date of the HDB
    std::string m_solver;             ///< Solver which computed the hydrodynamic data base (NEMOH/HELIOS)

    double m_version;                 ///< Version of the HDB file
    double m_gravityAcceleration;     ///< Gravity coming from the HDB
    double m_waterDensity;            ///< Water density coming from the HDB
    double m_waterDepth;              ///< Water depth coming from the HDB
    double m_normalizationLength;     ///< Normalization length coming from the HDB

    int m_nbody;                      ///< Number of bodies in interaction considered in the HDB

    std::vector<std::unique_ptr<Body>> m_bodies;      ///< List of BEM body database

    Discretization1D m_frequencyDiscretization;       ///< Wave frequency discretization
    Discretization1D m_waveDirectionDiscretization;   ///< Wave direction discretization
    Discretization1D m_timeDiscretization;            ///< Time samples


    void Import_HDF5_v3(const std::string &HDF5_file);

//    void ReadExcitation(HighFive::Group* group, Body* body);

    void ReadExcitation(const HighFive::File &HDF5_file, const std::string &path, Body* body);

    void ReadRadiation(const HighFive::File &HDF5_file, const std::string &path, Body* body);

  };

} // namespace HDB5_io

#endif //HDB5_IO_HYDRODYNAMICDATABASE_H

//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_HYDRODYNAMICDATABASE_H
#define HDB5_IO_HYDRODYNAMICDATABASE_H

#include <iostream>
#include <memory>

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorN.h"

namespace HDB5_io {

  // Forward declarations

  class WaveDrift;

  class Body;

  /**
  * \class HydrodynamicDataBase
  * \brief Class for storing a hydrodynamic database.
  * All components are stored, even null DOF components, resulting in 6 DOF vectors and 6x6 matrices.
  */
  class HydrodynamicDataBase {

   public:

    // Accessors

    void SetCreationDate(std::string date);

    std::string GetCreationDate() const;

    void SetSolver(std::string solver);

    std::string GetSolver() const;

    void SetVersion(double version);

    double GetVersion() const;

    void SetGravityAcceleration(double g);

    double GetGravityAcceleration() const;

    void SetWaterDensity(double rho);

    double GetWaterDensity() const;

    void SetWaterDepth(double h);

    double GetWaterDepth() const;

    void SetNormalizationLength(double L);

    double GetNormalizationLength() const;

    Body *NewBody(unsigned int id, const std::string &name);

    Body *GetBody(int id) const;

    void SetNbBodies(int nb);

    int GetNbBodies() const;

    void SetFrequencyDiscretization(const mathutils::VectorN<double> &frequency) {m_frequencyDiscretization = frequency;}

    void SetWaveDirectionDiscretization(const mathutils::VectorN<double> &directions) {m_waveDirectionDiscretization = directions;}

    void SetTimeDiscretization(const mathutils::VectorN<double> &time) {m_timeDiscretization = time;}

    mathutils::VectorN<double> GetFrequencyDiscretization() const {return m_frequencyDiscretization;}

    mathutils::VectorN<double> GetWaveDirectionDiscretization() const {return m_waveDirectionDiscretization;}

    mathutils::VectorN<double> GetTimeDiscretization() const {return m_timeDiscretization;}

    void SetWaveDrift(const std::shared_ptr<WaveDrift> &wavedrift) {
      m_waveDrift = wavedrift;
    }

    WaveDrift* GetWaveDrift() const {return m_waveDrift.get();}


   private:

    std::string m_creationDate;       ///< Creation date of the HDB
    std::string m_solver;             ///< Solver which computed the hydrodynamic data base (NEMOH/HELIOS)

    double m_version;                 ///< Version of the HDB file
    double m_gravityAcceleration;     ///< Gravity coming from the HDB
    double m_waterDensity;            ///< Water density coming from the HDB
    double m_waterDepth;              ///< Water depth coming from the HDB
    double m_normalizationLength;     ///< Normalization length coming from the HDB

    int m_nbody;                      ///< Number of bodies in interaction considered in the HDB
    std::vector<std::shared_ptr<Body>> m_bodies;      ///< List of BEM body database

    mathutils::VectorN<double> m_frequencyDiscretization;
    mathutils::VectorN<double> m_timeDiscretization;
    mathutils::VectorN<double> m_waveDirectionDiscretization;

    std::shared_ptr<WaveDrift> m_waveDrift;            ///< wave drift components

  };

} // namespace HDB5_io

#endif //HDB5_IO_HYDRODYNAMICDATABASE_H

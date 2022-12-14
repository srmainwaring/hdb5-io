//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_HYDRODYNAMICDATABASE_H
#define HDB5_IO_HYDRODYNAMICDATABASE_H

#include <iostream>
#include <memory>

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorN.h"

namespace hdb5_io {

  // Forward declarations

  class WaveDrift;
  class Kochin;

  class Body;

  /**
  * \class HydrodynamicDataBase
  * \brief Class for storing a hydrodynamic database. Contains basic information (version, date of creation, solver, ...)
   * and body and wave drift container class containing all related hydrodynamic data.
  */
  class HydrodynamicDataBase {

   public:

    // Accessors

    void SetCreationDate(std::string date);

    std::string GetCreationDate() const;

    void SetSolver(std::string solver);

    std::string GetSolver() const;

    void SetNormalizedVersionString(std::string commit_hash);

    std::string GetNormalizedVersionString() const;

    bool IsNormalizedVersionString() const;

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

    /// Add a body to the hydrodynamic database
    /// \param id index of the body in the database
    /// \param name name o the body, must be unique
    /// \return body created
    Body *NewBody(unsigned int id, const std::string &name);

    /// Get id-th body
    /// \param id index of the body in the database
    /// \return
    Body *GetBody(int id) const;

    void SetNbBodies(int nb);

    int GetNbBodies() const;

    void SetFrequencyDiscretization(const mathutils::VectorN<double> &frequency);

    void SetWaveDirectionDiscretization(const mathutils::VectorN<double> &directions);

    void SetTimeDiscretization(const mathutils::VectorN<double> &time);

    mathutils::VectorN<double> GetFrequencyDiscretization() const;

    mathutils::VectorN<double> GetWaveDirectionDiscretization() const;

    mathutils::VectorN<double> GetTimeDiscretization() const;

    double GetMinFrequency() const;

    double GetMaxFrequency() const;

    void SetWaveDrift(const std::shared_ptr<WaveDrift> &wavedrift);

    WaveDrift* GetWaveDrift() const;

    void SetKochin(const std::shared_ptr<Kochin> &kochin);

    Kochin* GetKochin() const;

    void SetVF();

    bool GetVF() const;

    void SetVFRelaxed(const int &relaxed);

    int GetVFRelaxed() const;

    void SetVFMaxOrder(const int &order);

    int GetVFMaxOrder() const;

    void SetVFTolerance(const double &tolerance);

    double GetVFTolerance() const;

    void SetWaveField();

    bool GetWaveField() const;

    void SetSymmetries();

    bool GetSymmetries() const;

    void SetSymBottom(const bool &sym_bottom);

    bool GetSymBottom() const;

    void SetSymXOZ(const bool &sym_xOz);

    bool GetSymXOZ() const;

    void SetSymYOZ(const bool &sym_yOz);

    bool GetSymYOZ() const;

    void SetExpertParameters();

    bool GetExpertParameters() const;

    void SetSurfaceIntegrationOrder(const int &order);

    int GetSurfaceIntegrationOrder() const;

    void SetGreenFunction(const std::string &green_function);

    std::string GetGreenFunction() const;

    void SetCrmax(const int &crmax);

    int GetCrmax() const;

//    void SetWaveReferencePoint(const double &x, const double &y);

//    void GetWaveReferencePoint(double &x, double &y) const;

    void IsXDerivative();

    bool GetIsXDerivative() const;

   protected:

    std::string m_creationDate;       ///< Creation date of the HDB
    std::string m_solver;             ///< Solver which computed the hydrodynamic data base (NEMOH/HELIOS)
    bool m_is_commit_hash = false;
    std::string m_commit_hash;        ///< Normalized commit hash.

    double m_version;                 ///< Version of the HDB file
    double m_gravityAcceleration;     ///< Gravity coming from the HDB
    double m_waterDensity;            ///< Water density coming from the HDB
    double m_waterDepth;              ///< Water depth coming from the HDB
    double m_normalizationLength;     ///< Normalization length coming from the HDB

    bool m_isVF = false;
    int m_VF_relaxed; /// Relaxed vector fitting (1) or not (0).
    int m_VF_max_order; /// Maximum order of the vector fitting.
    double m_VF_tolerance; /// Tolerance of the vector fitting.

    bool m_isWaveField = false;

    bool m_isSymmetries = false;
    bool m_sym_bottom; /// Bottom symmetry.
    bool m_sym_xoz; /// (xOz) symmetry.
    bool m_sym_yoz; /// (yOz) symmetry.

    bool m_isExpertParameters = false;
    int m_surface_integration_order; /// Surface integration order.
    std::string m_Green_function; /// Method for evaluating the Green's function.
    double m_Crmax; /// Distance criteria threshold for applying the far-field approximation for computing the influence coefficients.
    double m_wave_reference_point_x; // x-coordinate of the wave reference point (m).
    double m_wave_reference_point_y; // y-coordinate of the wave reference point (m).

    int m_nbody;                      ///< Number of bodies in interaction considered in the HDB
    // FIXME :: change the vector to unordered_map with index as key, to be sure which body to get with the GetBody method
    std::vector<std::shared_ptr<Body>> m_bodies;      ///< List of BEM body database

    mathutils::VectorN<double> m_frequencyDiscretization;
    mathutils::VectorN<double> m_timeDiscretization;
    mathutils::VectorN<double> m_waveDirectionDiscretization;

    std::shared_ptr<WaveDrift> m_waveDrift;            ///< wave drift components
    std::shared_ptr<Kochin> m_kochin;            ///< Kochin functions.

    bool m_is_x_derivative = false; /// Presence if the x-derivative of the physical quantities.

  };

} // namespace hdb5_io

#endif //HDB5_IO_HYDRODYNAMICDATABASE_H

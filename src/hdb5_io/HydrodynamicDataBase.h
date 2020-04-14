//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_HYDRODYNAMICDATABASE_H
#define HDB5_IO_HYDRODYNAMICDATABASE_H

#include <iostream>
#include <memory>
#include "Discretization1D.h"
#include "Body.h"
#include "WaveDrift.h"

#include "MathUtils/Vector3d.h"
#include "highfive/H5Group.hpp"

namespace HDB5_io {

  // Forward declarations

  class WaveDrift;

  /**
  * \class HydrodynamicDataBase
  * \brief Class for storing a hydrodynamic database. Can import and export database from HDF5 data format file.
  * All components are stored, even null DOF components, resulting in 6 DOF vectors and 6x6 matrices.
  */
  class HydrodynamicDataBase {

   public:

    /// Import a hydrodynamic database from a HDF5 format file
    /// \param HDF5_file file containing the hydrodynamic database
    void Import_HDF5(const std::string &HDF5_file);

    /// Export the hydrodynamic database to a HDF5 format file
    /// \param HDF5_file file to export the hydrodynamic database
    void Export_HDF5(const std::string &HDF5_file);


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

    std::vector<std::shared_ptr<Body>> GetBodies();

    Body *GetBody(int id) const;

    void SetNbBodies(int nb);

    int GetNbBodies() const;

    void SetFrequencyDiscretization(double wmin, double wmax, unsigned int nw);

    void SetWaveDirectionDiscretization(double tmin, double tmax, unsigned int nt);

    void SetTimeDiscretization(double tmin, double tmax, unsigned int nt);

    Discretization1D GetFrequencyDiscretization() const;

    Discretization1D GetWaveDirectionDiscretization() const;

    Discretization1D GetTimeDiscretization() const;

//    void SetFrequencyDiscretization(const mathutils::VectorN<double> &frequency) {m_frequencies = frequency;}
//
//    void SetWaveDirectionDiscretization(const mathutils::VectorN<double> &directions) {m_waveDirections = directions;}
//
//    void SetTimeDiscretization(const mathutils::VectorN<double> &time) {m_time = time;}
//
//    mathutils::VectorN<double> GetFrequencyDiscretization() const {return m_frequencies;}
//
//    mathutils::VectorN<double> GetWaveDirectionDiscretization() const {return m_waveDirections;}
//
//    mathutils::VectorN<double> GetTimeDiscretization() const {return m_time;}



    /// Set the wave drift coefficient database
    void SetWaveDrift(const std::string &name, const Eigen::MatrixXd &data);

    void SetWaveDrift(const std::shared_ptr<WaveDrift> &wavedrift) {
      m_waveDrift = wavedrift;
    }

//    std::shared_ptr<WaveDrift> GetWaveDrift() const;


   private:

    enum excitationType {
      Diffraction, Froude_Krylov
    };

    std::string m_creationDate;       ///< Creation date of the HDB
    std::string m_solver;             ///< Solver which computed the hydrodynamic data base (NEMOH/HELIOS)

    double m_version;                 ///< Version of the HDB file
    double m_gravityAcceleration;     ///< Gravity coming from the HDB
    double m_waterDensity;            ///< Water density coming from the HDB
    double m_waterDepth;              ///< Water depth coming from the HDB
    double m_normalizationLength;     ///< Normalization length coming from the HDB

    int m_nbody;                      ///< Number of bodies in interaction considered in the HDB
    std::vector<std::shared_ptr<Body>> m_bodies;      ///< List of BEM body database

    Discretization1D m_frequencyDiscretization;       ///< Wave frequency discretization
    Discretization1D m_waveDirectionDiscretization;   ///< Wave direction discretization
    Discretization1D m_timeDiscretization;            ///< Time samples

    mathutils::VectorN<double> m_frequencies;
    mathutils::VectorN<double> m_time;
    mathutils::VectorN<double> m_waveDirections;

    std::shared_ptr<WaveDrift> m_waveDrift;            ///< wave drift components

    /// Import a hydrodynamic database from a HDF5 file in version 3.xx
    /// \param HDF5_file file containing the hydrodynamic database in version 3.xx
    void Import_HDF5_v2(const std::string &HDF5_file);

    /// Read the excitation components
    /// \param type excitation type (Diffraction or Froude_Krylov
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    void ReadExcitation(excitationType type, const HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Write the excitation components
    /// \param type excitation type (Diffraction or Froude_Krylov
    /// \param HDF5_file file to export the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body containing the components
    void WriteExcitation(excitationType type, HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Read the radiation components
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    void ReadRadiation(const HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Write the radiation components
    /// \param HDF5_file file to export the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body containing the components
    void WriteRadiation(HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Read the response amplitude operators
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    void ReadRAO(const HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Write the response amplitude operators
    /// \param HDF5_file file to export the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body containing the components
    void WriteRAO(HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Read the components (added mass/damping r / impulse response functions)
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param radiationMask radiation mask
    /// \return matrix containing the components
    std::vector<Eigen::MatrixXd> ReadComponents(const HighFive::File &HDF5_file, const std::string &path,
                                                Eigen::MatrixXi radiationMask);

    /// Read the mesh contained in the HDF5 file
    /// \param HDF5_file file containing the mesh
    /// \param path path to the mesh in the file
    /// \param body body to which store the mesh
    void ReadMesh(HighFive::File &HDF5_file, const std::string &path, Body *body);

    void ReadWaveDrift(HighFive::File &HDF5_file);

    void WriteWaveDrift(HighFive::File &HDF5_file);

  };

} // namespace HDB5_io

#endif //HDB5_IO_HYDRODYNAMICDATABASE_H

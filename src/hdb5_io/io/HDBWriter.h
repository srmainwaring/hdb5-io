//
// Created by lletourn on 15/04/20.
//

#ifndef HDB5_IO_HDBWRITER_H
#define HDB5_IO_HDBWRITER_H

#include <string>

#include <highfive/H5File.hpp>

namespace HDB5_io {

  // Forward Declarations
  class HydrodynamicDataBase;
  class Body;

  class HDBWriter {

   public:

    explicit HDBWriter(HydrodynamicDataBase* hdb) : m_hdb(hdb) {}

    virtual void Write(const std::string &filename) const;

   protected:

    enum excitationType {
      Diffraction, Froude_Krylov
    };

    HydrodynamicDataBase* m_hdb;

    virtual void WriteHDBBasics(HighFive::File &HDF5_file) const;

    virtual void WriteDiscretizations(HighFive::File &HDF5_file) const;

    virtual void WriteBodyBasics(HighFive::File &HDF5_file, const std::string &path, Body *body) const;

    virtual void WriteMesh(HighFive::File &HDF5_file, const std::string &path, Body *body) const;

    /// Write the excitation components
    /// \param type excitation type (Diffraction or Froude_Krylov
    /// \param HDF5_file file to export the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body containing the components
    virtual void WriteExcitation(excitationType type, HighFive::File &HDF5_file, const std::string &path, Body *body) const;

    /// Write the radiation components
    /// \param HDF5_file file to export the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body containing the components
    virtual void WriteRadiation(HighFive::File &HDF5_file, const std::string &path, Body *body) const;

    /// Write the response amplitude operators
    /// \param HDF5_file file to export the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body containing the components
    virtual void WriteRAO(HighFive::File &HDF5_file, const std::string &path, Body *body) const;

    virtual void WriteWaveDrift(HighFive::File &HDF5_file) const;

  };

  void export_HDB(const std::string &filename, HydrodynamicDataBase *hdb);

}

#endif //HDB5_IO_HDBWRITER_H

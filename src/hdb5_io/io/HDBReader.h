//
// Created by lletourn on 14/04/20.
//

#ifndef HDB5_IO_HDBREADER_H
#define HDB5_IO_HDBREADER_H

#include <string>
#include <Eigen/Dense>

#include "highfive/H5Group.hpp"

namespace HDB5_io {

  // Forward Declaration
  class HydrodynamicDataBase;
  class Body;

  class HDBReader {

   public:

    explicit HDBReader(HydrodynamicDataBase* hdb) : m_hdb(hdb) {}

    virtual void Read(const std::string& filename);

   protected:

    enum excitationType {
      Diffraction, Froude_Krylov
    };

//    double m_version;

    HydrodynamicDataBase* m_hdb;

    virtual void ReadHDBBasics(const HighFive::File &HDF5_file);

    virtual Body* ReadBodyBasics(const HighFive::File &HDF5_file, const std::string &path);

    virtual void ReadDiscretizations(const HighFive::File &file) = 0;

    /// Read the excitation components
    /// \param type excitation type (Diffraction or Froude_Krylov
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    virtual void ReadExcitation(excitationType type, const HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Read the radiation components
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    virtual void ReadRadiation(const HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Read the response amplitude operators
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    virtual void ReadRAO(const HighFive::File &HDF5_file, const std::string &path, Body *body);

    /// Read the components (added mass/radiation damping/ impulse response functions)
    /// \param HDF5_file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param radiationMask radiation mask
    /// \return matrix containing the components
    virtual std::vector<Eigen::MatrixXd> ReadComponents(const HighFive::File &HDF5_file, const std::string &path,
                                                        Eigen::MatrixXi radiationMask);

    /// Read the mesh contained in the HDF5 file
    /// \param HDF5_file file containing the mesh
    /// \param path path to the mesh in the file
    /// \param body body to which store the mesh
    virtual void ReadMesh(HighFive::File &HDF5_file, const std::string &path, Body *body);

    virtual void ReadWaveDrift(HighFive::File &HDF5_file);

  };

  std::shared_ptr<HydrodynamicDataBase> import_HDB(const std::string& filename);




  class HDBReader_v2 : public HDBReader {

   public:

    explicit HDBReader_v2(HydrodynamicDataBase* hdb) : HDBReader(hdb) {}

   protected:

    void ReadDiscretizations(const HighFive::File &file) override;


  };



  class HDBReader_v3 : public HDBReader {

   public:

    explicit HDBReader_v3(HydrodynamicDataBase* hdb) : HDBReader(hdb) {}

   protected:

    void ReadDiscretizations(const HighFive::File &file) override;


  };

} // end namespace HDB5_io
#endif //HDB5_IO_HDBREADER_H

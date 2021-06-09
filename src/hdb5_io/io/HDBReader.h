//
// Created by lletourn on 14/04/20.
//

#ifndef HDB5_IO_HDBREADER_H
#define HDB5_IO_HDBREADER_H

#include <string>
#include <Eigen/Dense>

#include <highfive/H5Group.hpp>

namespace hdb5_io {

  // Forward Declaration
  class HydrodynamicDataBase;

  class Body;

  /**
  * \class HDBReader
  * \brief Class for reading a hydrodynamic database from a .hdb5 file. 
  */
  class HDBReader {

   public:

    /// Constructor of the HDBReader
    /// \param hdb hydrodynamic database in which store the data
    explicit HDBReader(HydrodynamicDataBase *hdb) : m_hdb(hdb) {}

    /// Read a hdb5 file
    /// \param filename name of the file to read
    virtual void Read(const std::string &filename);

   protected:

    enum excitationType {
      Diffraction, Froude_Krylov
    };

//    double m_version;

    HydrodynamicDataBase *m_hdb;  ///< hydrodynamic database to store the data

    /// Read basic information contained in the hydrodynamic database (version, date of creation, solver, etc.)
    /// \param file file containing the hydrodynamic database
    virtual void ReadHDBBasics(const HighFive::File &file);

    /// Read basic information related to the body given in the path
    /// \param file file containing the hydrodynamic database
    /// \param path path to the body data in the hdb5
    /// \return container of the body hydrodynamic data
    virtual Body *ReadBodyBasics(const HighFive::File &file, const std::string &path);

    /// Read the wave direction, frequency, and time discretizations 
    /// \param file file containing the hydrodynamic database
    virtual void ReadDiscretizations(const HighFive::File &file) = 0;

    /// Read the symmetry data
    /// \param file file containing the hydrodynamic database
    virtual void ReadSymmetries(HighFive::File &file);

    /// Read the excitation components
    /// \param type excitation type (Diffraction or Froude_Krylov
    /// \param file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    virtual void
    ReadExcitation(excitationType type, const HighFive::File &file, const std::string &path, Body *body);

    /// Read the radiation components
    /// \param file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    virtual void ReadRadiation(const HighFive::File &file, const std::string &path, Body *body);

    /// Read the response amplitude operators
    /// \param file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    virtual void ReadRAO(const HighFive::File &file, const std::string &path, Body *body);

    /// Read the components (added mass/radiation damping/ impulse response functions)
    /// \param file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param radiationMask radiation mask
    /// \return matrix containing the components
    virtual std::vector<Eigen::MatrixXd> ReadComponents(const HighFive::File &file, const std::string &path,
                                                        Eigen::Matrix<bool,6,6> radiationMask);

    /// Read the wave drift data
    /// \param file file containing the hydrodynamic database
    virtual void ReadWaveDrift(HighFive::File &file) = 0;

    /// Read the mesh contained in the HDF5 file
    /// \param file file containing the mesh
    /// \param path path to the mesh in the file
    /// \param body body to which store the mesh
    virtual void ReadMesh(HighFive::File &file, const std::string &path, Body *body);

    /// Read the wave field data
    /// \param file file containing the hydrodynamic database
    virtual void ReadWaveField(HighFive::File &file);

    /// Read the wave drift data
    /// \param file file containing the hydrodynamic database
    virtual void ReadVectorFitting(HighFive::File &file);

    /// This method reads the expert numerical parameters.
    virtual void ReadExpertNumericalParameters(HighFive::File &file);

    /// Read the wave drift components, from the path given, for the i-th body
    /// \param file file containing the hydrodynamic database
    /// \param path path in the hdb5 file
    /// \param i index of the body
    /// \return wave drift components (surge, sway or yaw) in Eigen::vectorXd format
    virtual Eigen::VectorXd
    ReadWaveDriftComponents(HighFive::File &file, const std::string &path, unsigned int i) = 0;

  };

  /// Import a hydrodynamic database from a .hdb5 file
  /// \param filename name of the file to import
  /// \return hydrodynamic database
  std::shared_ptr<HydrodynamicDataBase> import_HDB(const std::string &filename);


  /**
  * \class HDBReader_v2
  * \brief Class for reading a hydrodynamic database from a .hdb5 file, stored in version 2.xx
  */
  class HDBReader_v2 : public HDBReader {

   public:

    explicit HDBReader_v2(HydrodynamicDataBase *hdb) : HDBReader(hdb) {}

   protected:

    void ReadDiscretizations(const HighFive::File &file) override;

    /// Read basic information related to the body given in the path
    /// \param file file containing the hydrodynamic database
    /// \param path path to the body data in the hdb5
    /// \return container of the body hydrodynamic data
    Body *ReadBodyBasics(const HighFive::File &file, const std::string &path) override;

    /// Read the wave drift data
    /// \param file file containing the hydrodynamic database
    void ReadWaveDrift(HighFive::File &file) override;

    Eigen::VectorXd ReadWaveDriftComponents(HighFive::File &file, const std::string &path, unsigned int i) override;

  };


  /**
  * \class HDBReader_v3
  * \brief Class for reading a hydrodynamic database from a .hdb5 file, stored in version 3.xx
  */
  class HDBReader_v3 : public HDBReader {

   public:

    explicit HDBReader_v3(HydrodynamicDataBase *hdb) : HDBReader(hdb) {}

   protected:

    void ReadDiscretizations(const HighFive::File &file) override;

    /// Read basic information related to the body given in the path
    /// \param file file containing the hydrodynamic database
    /// \param path path to the body data in the hdb5
    /// \return container of the body hydrodynamic data
    Body *ReadBodyBasics(const HighFive::File &file, const std::string &path) override;

    Eigen::VectorXd ReadWaveDriftComponents(HighFive::File &file, const std::string &path, unsigned int i) override;

    /// Read the wave drift data
    /// \param file file containing the hydrodynamic database
    virtual void ReadWaveDrift(HighFive::File &file);

    /// Read the radiation components
    /// \param file file containing the hydrodynamic database
    /// \param path path to the components in the file
    /// \param body body to which store the components
    void ReadRadiation(const HighFive::File &file, const std::string &path, Body *body) override;

  };

} // end namespace hdb5_io
#endif //HDB5_IO_HDBREADER_H

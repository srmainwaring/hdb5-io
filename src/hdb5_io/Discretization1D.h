//
// Created by lletourn on 26/02/20.
//

#ifndef HDB5_IO_DISCRETIZATION1D_H
#define HDB5_IO_DISCRETIZATION1D_H

#include <vector>

namespace HDB5_io {


  // ----------------------------------------------------------
  // Discretization1D
  // ----------------------------------------------------------

  /**
   * \class Discretization1D
   * \brief Class for the linear discretization (frequency and angle).
   */
  class Discretization1D {
   private:
    double m_xmin = 0.;         ///< Min value of the discretization
    double m_xmax = 0.;         ///< Max value of the discretization
    unsigned int m_nx = 0;      ///< Number of discrete value

   public:

    Discretization1D() = default;

    /// Constructor with specified parameters
    /// \param xmin Min abscissa
    /// \param xmax Max abscissa
    /// \param nx Number of discrete value
    Discretization1D(double xmin, double xmax, unsigned int nx)
        : m_xmin(xmin), m_xmax(xmax), m_nx(nx) {}

    /// Set / Get the min value of the discretization vector
    double GetMin() const { return m_xmin; }

    void SetMin(double xmin) { m_xmin = xmin; }

    /// Set / Get the max value of the discretization vector
    double GetMax() const { return m_xmax; }

    void SetMax(double xmax) { m_xmax = xmax; }

    /// Set / Get the number of discrete value
    unsigned int GetNbSample() const { return m_nx; }

    void SetNbSample(unsigned int nx) { m_nx = nx; }

    /// Return the discrete vector
    /// \return
    std::vector<double> GetVector() const;

    /// Define the step the discretization vector
    /// \param delta Step value
    void SetStep(double delta);

    /// Return the step of the discretization vector
    /// \return Step value
    double GetStep() const;
  };


} // namespace HDB5_io


#endif //HDB5_IO_DISCRETIZATION1D_H

//
// Created by lletourn on 27/02/20.
//

#ifndef HDB5_IO_MASK_H
#define HDB5_IO_MASK_H

#include "vector"

#include "MathUtils/MathUtils.h"

namespace HDB5_io {

   /**
   * \class Mask
   * \brief Class for the force and motions masks
   */
  class Mask {

   private:
    unsigned int m_nbDOF;                   ///< Number of degree of freedom of the body
    std::vector<int> m_listDOF;             ///< List of degree of freedom of the body
    mathutils::Vector6d<bool> m_mask;       ///< Mask applied on the degree of freedom
    mathutils::MatrixMN<double> m_matrix;   ///<

   public:
    /// Default constructor of the DOF Mask
    Mask() {}

    /// Define the DOF mask from vector
    /// \param mask DOF Mask vector
    void SetMask(mathutils::Vector6d<int> mask);

    /// Return the DOF Mask
    /// \return DOF Mask
    mathutils::Vector6d<bool> GetMask() const;

    /// Return the mask for a specific degree of freedom
    /// \param imode Degree of freedom
    /// \return True if the mask is active, False otherwise
    bool GetMask(unsigned int imode) const { return m_mask[imode]; }

    ///
    /// \return
    mathutils::MatrixMN<double> GetMatrix() const;

    /// Return the number of activated DOF in the mask
    /// \return Number of activated DOF
    unsigned int GetNbDOF() const { return m_nbDOF; }

    /// Return the list of activated DOF
    /// \return List of activated DOF
    std::vector<int> GetListDOF() const { return m_listDOF; }
  };


} // namespace HDB5_io


#endif //HDB5_IO_MASK_H

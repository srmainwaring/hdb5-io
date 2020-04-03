//
// Created by lletourn on 27/02/20.
//

#ifndef HDB5_IO_MODE_H
#define HDB5_IO_MODE_H

#include "MathUtils/Vector3d.h"

namespace HDB5_io {
  /**
   * \class Mode
   * \brief
   */
  class Mode {

   public:
    enum TYPE {
      LINEAR,
      ANGULAR
    };

   private:
    TYPE m_type;
    mathutils::Vector3d<double> m_direction;
    mathutils::Vector3d<double> m_position;
    bool m_active = true;

   public:
    Mode() = default;

    void SetTypeLINEAR() { m_type = LINEAR; }

    void SetTypeANGULAR() { m_type = ANGULAR; }

    TYPE GetType() const { return m_type; }

    void SetDirection(mathutils::Vector3d<double> &direction) { m_direction = direction; }

    mathutils::Vector3d<double> GetDirection() const { return m_direction; }

    void SetPointPosition(mathutils::Vector3d<double> &position) { m_position = position; }

    mathutils::Vector3d<double> GetPointPosition() const { return m_position; }

    void Activate() { m_active = true; }

    void Deactivate() { m_active = false; }

    bool IsActive() const { return m_active; }

  };

//  typedef Mode ForceMode;
//  typedef Mode MotionMode;

} // namespace HDB5_io


#endif //HDB5_IO_MODE_H

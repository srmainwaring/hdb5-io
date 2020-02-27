//
// Created by lletourn on 27/02/20.
//

#ifndef HDB5_IO_MODE_H
#define HDB5_IO_MODE_H

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
    Direction m_direction;
    Position m_position;
    bool m_active = true;

   public:
    Mode() = default;

    void SetTypeLINEAR() { m_type = LINEAR; }

    void SetTypeANGULAR() { m_type = ANGULAR; }

    TYPE GetType() const { return m_type; }

    void SetDirection(Direction &direction) { m_direction = direction; }

    Direction GetDirection() const { return m_direction; }

    void SetPointPosition(Position &position) { m_position = position; }

    Position GetPointPosition() const { return m_position; }

    void Activate() { m_active = true; }

    void Deactivate() { m_active = false; }

    bool IsActive() const { return m_active; }

  };

//  typedef Mode ForceMode;
//  typedef Mode MotionMode;

} // namespace HDB5_io


#endif //HDB5_IO_MODE_H

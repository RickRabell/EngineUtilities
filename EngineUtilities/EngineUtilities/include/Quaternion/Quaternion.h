#pragma once

#include "../Vectors/Vector3.h"
#include "../Matrix/Matrix4x4.h"
#include "../Utilities/EngineMath.h"

/**
  * @file Quaternion.h
  * @brief Defines the Quaternion class for representing 3D rotations.
*/

namespace EngineUtilities {

  /**
    * @class Quaternion
    * @brief Represents a rotation in 3D space using a four-component quaternion.
  */
  class Quaternion {
  public:
    /**
      * @brief The scalar (w) and vector (x, y, z) components of the quaternion.
    */
    float m_w, m_x, m_y, m_z;

    /**
      * @brief Default constructor. Initializes to the identity quaternion (no rotation).
    */
    Quaternion() : m_w(1.0f), m_x(0.0f), m_y(0.0f), m_z(0.0f) {}

    /**
      * @brief Constructor with all 4 components.
      * @param w_ The scalar component.
      * @param x_ The x component of the vector part.
      * @param y_ The y component of the vector part.
      * @param z_ The z component of the vector part.
    */
    Quaternion(float w_, float x_, float y_, float z_) : m_w(w_), m_x(x_), m_y(y_), m_z(z_) {}

    /**
      * @brief Returns the identity quaternion.
      * @return Identity quaternion.
    */
    static Quaternion 
    Identity() {
      return Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
    }

    /**
      * @brief Creates a quaternion from an axis and an angle.
      * @param axis The axis of rotation (should be normalized).
      * @param radians The angle of rotation in radians.
      * @return A quaternion representing the rotation.
    */
    static Quaternion 
    FromAxisAngle(const Vector3& axis, float radians) {
      float halfAngle = radians * 0.5f;
      float s = EngineMath::sin(halfAngle);
      
      return Quaternion(
        EngineMath::cos(halfAngle),
        axis.m_x * s,
        axis.m_y * s,
        axis.m_z * s
      );
    }

    /**
      * @brief Calculates the squared length (magnitude) of the quaternion.
      * @return The squared length.
    */
    float 
    LengthSquared() const {
      return m_w * m_w + m_x * m_x + m_y * m_y + m_z * m_z;
    }

    /**
      * @brief Calculates the length (magnitude) of the quaternion.
      * @return The length.
    */
    float 
    Length() const {
      return EngineMath::sqrt(LengthSquared());
    }

    /**
      * @brief Normalizes this quaternion in place.
    */
    void 
    Normalize() {
      float len = Length();
      if (EngineMath::absf(len) > EngineMath::epsilon) {
        float invLen = 1.0f / len;
        m_w *= invLen;
        m_x *= invLen;
        m_y *= invLen;
        m_z *= invLen;
      }
    }

    /**
      * @brief Returns a normalized copy of this quaternion.
      * @return The normalized quaternion.
    */
    Quaternion 
    Normalized() const {
      Quaternion q = *this;
      q.Normalize();
      
      return q;
    }

    /**
      * @brief Returns the conjugate of this quaternion.
      * @return The conjugate quaternion.
    */
    Quaternion 
    Conjugate() const {
      return Quaternion(m_w, -m_x, -m_y, -m_z);
    }

    /**
      * @brief Returns the inverse of this quaternion.
      * @return The inverse quaternion.
    */
    Quaternion 
    Inverse() const {
      float lenSq = LengthSquared();
      
      if (EngineMath::absf(lenSq) < EngineMath::epsilon) {
        return Quaternion::Identity();
      }
      
      return Conjugate() * (1.0f / lenSq);
    }

    /**
      * @brief Converts the quaternion to a 4x4 rotation matrix.
      * @return The equivalent rotation matrix.
    */
    Matrix4x4 
    ToMatrix() const {
      Quaternion q = Normalized();
      float xx = q.m_x * q.m_x, yy = q.m_y * q.m_y, zz = q.m_z * q.m_z;
      float xy = q.m_x * q.m_y, xz = q.m_x * q.m_z, yz = q.m_y * q.m_z;
      float wx = q.m_w * q.m_x, wy = q.m_w * q.m_y, wz = q.m_w * q.m_z;

      return Matrix4x4(
        1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy), 0,
        2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx), 0,
        2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy), 0,
        0, 0, 0, 1
      );
    }

    /**
      * @brief Adds two quaternions component-wise.
      * @param other The quaternion to add.
      * @return The resulting quaternion.
    */
    Quaternion 
    operator+(const Quaternion& other) const {
      return Quaternion(
        m_w + other.m_w,
        m_x + other.m_x,
        m_y + other.m_y,
        m_z + other.m_z
      );
    }

    /**
      * @brief Substracts two quaternions component-wise.
      * @param other The quaternion to substract.
      * @return The resulting quaternion.
    */
    Quaternion
      operator-(const Quaternion& other) const {
      return Quaternion(
        m_w - other.m_w,
        m_x - other.m_x,
        m_y - other.m_y,
        m_z - other.m_z
      );
    }

    /**
      * @brief Multiplies (concatenates) this quaternion with another.
      * @param other The quaternion to multiply by.
      * @return The resulting quaternion.
    */
    Quaternion 
    operator*(const Quaternion& other) const {
      return Quaternion(
        m_w * other.m_w - m_x * other.m_x - m_y * other.m_y - m_z * other.m_z,
        m_w * other.m_x + m_x * other.m_w + m_y * other.m_z - m_z * other.m_y,
        m_w * other.m_y - m_x * other.m_z + m_y * other.m_w + m_z * other.m_x,
        m_w * other.m_z + m_x * other.m_y - m_y * other.m_x + m_z * other.m_w
      );
    }

    /**
      * @brief Rotates a vector by this quaternion.
      * @param v The vector to rotate.
      * @return The rotated vector.
    */
    Vector3 
    operator*(const Vector3& v) const {
      Quaternion p(0, v.m_x, v.m_y, v.m_z);
      Quaternion result = (*this) * p * this->Inverse();
      
      return Vector3(result.m_x, result.m_y, result.m_z);
    }

    /**
      * @brief Multiplies the quaternion by a scalar.
      * @param scalar The scalar value.
      * @return The scaled quaternion.
    */
    Quaternion 
    operator*(float scalar) const {
      return Quaternion(m_w * scalar, m_x * scalar, m_y * scalar, m_z * scalar);
    }

  };
}
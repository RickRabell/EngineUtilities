#pragma once

#include "../Utilities/EngineMath.h"

namespace EngineUtilities {
  /*
    *  @brief Represents a 4D vector and provides arithmetic and geometric operations.
  */
  class Vector4 {
  public:
    /*
      *  @brief X component of the vector.
    */
    float m_x;
    /*
      *  @brief Y component of the vector.
    */
    float m_y;
    /*
      *  @brief Z component of the vector.
    */
    float m_z;
    /*
      *  @brief W component of the vector.
    */
    float m_w;

    /*
      *  @brief Default constructor. Initializes all components to zero.
    */
    Vector4() { 
      m_x = 0; 
      m_y = 0; 
      m_z = 0; 
      m_w = 0; 
    }
    /*
      *  @brief Initializes all components to the given value.
      *  @param value Value to assign to all components.
    */
    Vector4(float value) { 
      m_x = value; 
      m_y = value; 
      m_z = value; 
      m_w = value; 
    }
    /*
      *  @brief Initializes the vector with the given components.
      *  @param x X component.
      *  @param y Y component.
      *  @param z Z component.
      *  @param w W component.
    */
    Vector4(float x, float y, float z, float w) { 
      m_x = x; 
      m_y = y; 
      m_z = z; 
      m_w = w; 
    }

    /*
      *  @brief Adds two vectors.
      *  @param other Vector to add.
      *  @return Resulting vector.
    */
    Vector4 
    operator+(const Vector4& other) const {
      return Vector4(m_x + other.m_x, 
                     m_y + other.m_y, 
                     m_z + other.m_z, 
                     m_w + other.m_w);
    }

    /*
      *  @brief Subtracts two vectors.
      *  @param other Vector to subtract.
      *  @return Resulting vector.
    */
    Vector4 
    operator-(const Vector4& other) const {
      return Vector4(m_x - other.m_x, 
                     m_y - other.m_y, 
                     m_z - other.m_z, 
                     m_w - other.m_w);
    }

    /*
      *  @brief Multiplies the vector by a scalar.
      *  @param scalar Scalar value.
      *  @return Resulting vector.
    */
    Vector4 
    operator*(float scalar) const {
      return Vector4(m_x * scalar, 
                     m_y * scalar, 
                     m_z * scalar, 
                     m_w * scalar);
    }

    /*
      *  @brief Divides the vector by a scalar.
      *  @param scalar Scalar value.
      *  @return Resulting vector.
    */
    Vector4 
    operator/(float scalar) const {
      return Vector4(m_x / scalar, 
                     m_y / scalar, 
                     m_z / scalar, 
                     m_w / scalar);
    }

    /*
      *  @brief Adds another vector to this vector.
      *  @param other Vector to add.
      *  @return Reference to this vector.
    */
    Vector4& 
    operator+=(const Vector4& other) {
      m_x += other.m_x;
      m_y += other.m_y;
      m_z += other.m_z;
      m_w += other.m_w;
      
      return *this;
    }

    /*
      *  @brief Subtracts another vector from this vector.
      *  @param other Vector to subtract.
      *  @return Reference to this vector.
    */
    Vector4& 
    operator-=(const Vector4& other) {
      m_x -= other.m_x;
      m_y -= other.m_y;
      m_z -= other.m_z;
      m_w -= other.m_w;
      
      return *this;
    }

    /*
      *  @brief Multiplies this vector by another vector component-wise.
      *  @param other Vector to multiply.
      *  @return Reference to this vector.
    */
    Vector4& 
    operator*=(const Vector4& other) {
      m_x *= other.m_x;
      m_y *= other.m_y;
      m_z *= other.m_z;
      m_w *= other.m_w;
      
      return *this;
    }

    /*
      *  @brief Multiplies this vector by a scalar.
      *  @param scalar Scalar value.
      *  @return Reference to this vector.
    */
    Vector4& 
    operator*=(float scalar) {
      m_x *= scalar;
      m_y *= scalar;
      m_z *= scalar;
      m_w *= scalar;
      
      return *this;
    }

    /*
      *  @brief Checks if two vectors are equal.
      *  @param other Vector to compare.
      *  @return True if equal, false otherwise.
    */
    bool 
    operator==(const Vector4& other) const {
      return (m_x == other.m_x && 
              m_y == other.m_y && 
              m_z == other.m_z && 
              m_w == other.m_w);
    }

    /*
      *  @brief Checks if two vectors are not equal.
      *  @param other Vector to compare.
      *  @return True if not equal, false otherwise.
    */
    bool 
    operator!=(const Vector4& other) const {
      return (m_x != other.m_x || 
              m_y != other.m_y || 
              m_z != other.m_z || 
              m_w != other.m_w);
    }

    /*
      *  @brief Accesses a component by index.
      *  @param index Index of the component (0=x, 1=y, 2=z, 3=w).
      *  @return Reference to the component.
      *  @throws "Index out of range for Vector4" if index is invalid.
    */
    float& 
    operator[](int index) {
      switch (index) {
      case 0: return m_x;
      case 1: return m_y;
      case 2: return m_z;
      case 3: return m_w;
      default: throw "Index out of range for Vector4";
      }
    }

    /*
      *  @brief Calculates the length of a vector.
      *  @param other Vector to calculate length for.
      *  @return Length of the vector.
    */
    static float 
    length(const Vector4& other) {
      return EngineMath::sqrt(squaredLength(other));
    }

    /*
      *  @brief Calculates the squared length of a vector.
      *  @param other Vector to calculate squared length for.
      *  @return Squared length of the vector.
    */
    static float 
    squaredLength(const Vector4& other) {
      return dotProduct(other, other);
    }

    /*
      *  @brief Calculates the dot product of two vectors.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Dot product.
    */
    static float 
    dotProduct(const Vector4& a, const Vector4& b) {
      return (a.m_x * b.m_x) + 
             (a.m_y * b.m_y) + 
             (a.m_z * b.m_z) + 
             (a.m_w * b.m_w);
    }

    /*
      *  @brief Calculates the cross product of two vectors (ignores w component).
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Cross product vector.
    */
    Vector4 
    crossProduct(const Vector4& a, const Vector4& b) const {
      return Vector4(
        a.m_y * b.m_z - a.m_z * b.m_y,
        a.m_z * b.m_x - a.m_x * b.m_z,
        a.m_x * b.m_y - a.m_y * b.m_x,
        0.0f
      );
    }

    /*
      *  @brief Returns a normalized copy of the given vector.
      *  @param other Vector to normalize.
      *  @return Normalized vector. If the length of the vector is zero, returns a zero vector.
    */
    Vector4 
    normalize(const Vector4& other) const {
      float len = length(other);
      if (len == 0.0f) {
        // Return zero vector if length is zero to avoid division by zero
        return Vector4::Zero();
      }

      return Vector4(other.m_x / len, 
                     other.m_y / len, 
                     other.m_z / len, 
                     other.m_w / len);
    }

    /*
      *  @brief Returns a normalized copy of the given vector (non-const overload).
      *  @param other Vector to normalize.
      *  @return Normalized vector. If the length of the vector is zero, returns a zero vector.
    */
    static Vector4 
    normalize(Vector4& other) {
      float len = length(other);
      if (len == 0.0f) {
        // Return zero vector if length is zero to avoid division by zero
        return Vector4::Zero();
      }

      return Vector4(other.m_x / len, 
                     other.m_y / len, 
                     other.m_z / len, 
                     other.m_w / len);
    }

    /*
      *  @brief Normalizes this vector in place.
    */
    void 
    normalizeInPlace() {
      float len = length(*this);

      m_x /= len;
      m_y /= len;
      m_z /= len;
      m_w /= len;
    }

    /*
      *  @brief Calculates the distance between two vectors.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Distance between vectors.
    */
    static float 
    distance(const Vector4& a, const Vector4& b) {
      return length(b - a);
    }

    /*
      *  @brief Linearly interpolates between two vectors.
      *  @param a Start vector.
      *  @param b End vector.
      *  @param t Interpolation factor.
      *  @return Interpolated vector.
    */
    static Vector4 
    lerp(const Vector4& a, const Vector4& b, float t) {
      return a + (b - a) * t;
    }

    /*
      *  @brief Returns a zero vector.
      *  @return Vector with all components set to zero.
    */
    static Vector4 
    Zero() {
      return Vector4(0.0f, 0.0f, 0.0f, 0.0f);
    }

    /*
      *  @brief Returns a unit vector.
      *  @return Vector with all components set to one.
    */
    static Vector4 
    Unit() {
      return Vector4(1.0f, 1.0f, 1.0f, 1.0f);
    }

    /*
      *  @brief Spherically interpolates between two vectors.
      *  @param a Start vector.
      *  @param b End vector.
      *  @param t Interpolation factor.
      *  @return Interpolated vector.
    */
    static Vector4 
    slerp(const Vector4& a, const Vector4& b, float t) {
      float dot = dotProduct(a, b);
      dot = EngineMath::clamp(dot, -1.0f, 1.0f);
      float theta = EngineMath::aCos(dot) * t;

      Vector4 relative = b - a * dot;
      float relLen = length(relative);

      if (relLen < 1e-6f) {
        return a;
      }

      relative = Vector4::normalize(relative);

      return a * EngineMath::cos(theta) + relative * EngineMath::sin(theta);
    }

    /*
      *  @brief Reflects a direction vector off a surface with the given normal.
      *  @param direction Direction vector.
      *  @param normal Surface normal vector.
      *  @return Reflected vector.
    */
    static Vector4 
    reflect(const Vector4& direction, const Vector4& normal) {
      float dot = dotProduct(direction, normal);
      
      return direction - normal * (2.0f * dot);
    }

    /*
      *  @brief Projects vector a onto vector b.
      *  @param a Vector to project.
      *  @param b Vector to project onto.
      *  @return Projected vector.
    */
    static Vector4 
    project(const Vector4& a, const Vector4& b) {
      float dot = dotProduct(a, b);
      float lenSq = squaredLength(b);
      
      return b * (dot / lenSq);
    }

    /*
      *  @brief Clamps the components of a vector between min and max vectors.
      *  @param value Vector to clamp.
      *  @param min Minimum vector.
      *  @param max Maximum vector.
      *  @return Clamped vector.
    */
    static Vector4 
    clamp(const Vector4& value, const Vector4& min, const Vector4& max) {
      float x = EngineMath::clamp(value.m_x, min.m_x, max.m_x);
      float y = EngineMath::clamp(value.m_y, min.m_y, max.m_y);
      float z = EngineMath::clamp(value.m_z, min.m_z, max.m_z);
      float w = EngineMath::clamp(value.m_w, min.m_w, max.m_w);
      
      return Vector4(x, y, z, w);
    }

    /*
      *  @brief Calculates the angle between two vectors in radians.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Angle in radians.
    */
    static float 
    angle(const Vector4& a, const Vector4& b) {
      float dot = dotProduct(a, b);
      float lenA = length(a);
      float lenB = length(b);
      
      return EngineMath::aCos(dot / (lenA * lenB));
    }
  };
}
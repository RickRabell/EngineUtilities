#pragma once
#pragma once
#include "../Utilities/EngineMath.h"

namespace EngineUtilities {
  /**
    *  @class Vector3
    *  @brief Represents a 3D vector and provides arithmetic and geometric operations.
  */
  class Vector3 {
  public:
    /**
      *  @brief x component of the vector.
    */
    float m_x;
    /**
      *  @brief y component of the vector.
    */
    float m_y;
    /**
      *  @brief z component of the vector.
    */
    float m_z;

    /**
      *  @brief Default constructor. Initializes all components to zero.
    */
    Vector3() { 
      m_x = 0; 
      m_y = 0; 
      m_z = 0; 
    }
    /**
      *  @brief Initializes all components to the given value.
      *  @param value Value to assign to x, y, and z.
    */
    Vector3(float value) { 
      m_x = value; 
      m_y = value; 
      m_z = value; 
    }
    /**
      *  @brief Initializes vector with given x, y, and z values.
      *  @param x X component.
      *  @param y Y component.
      *  @param z Z component.
    */
    Vector3(float x, float y, float z) { 
      this->m_x = x; 
      this->m_y = y; 
      this->m_z = z; 
    }

    /**
      *  @brief Adds two vectors.
      *  @param other Vector to add.
      *  @return Resulting vector.
    */
    Vector3 
    operator+(const Vector3& other) const {
      return Vector3(m_x + other.m_x, 
                     m_y + other.m_y, 
                     m_z + other.m_z);
    }

    /**
      *  @brief Subtracts two vectors.
      *  @param other Vector to subtract.
      *  @return Resulting vector.
    */
    Vector3 
    operator-(const Vector3& other) const {
      return Vector3(m_x - other.m_x, 
                     m_y - other.m_y, 
                     m_z - other.m_z);
    }

    /**
      *  @brief Multiplies vector by a scalar.
      *  @param scalar Scalar value.
      *  @return Resulting vector.
    */
    Vector3 
    operator*(float scalar) const {
      return Vector3(m_x * scalar, 
                     m_y * scalar, 
                     m_z * scalar);
    }

    /**
      *  @brief Divides vector by a scalar.
      *  @param scalar Scalar value.
      *  @return Resulting vector.
    */
    Vector3 
    operator/(float scalar) const {
      return Vector3(m_x / scalar, 
                     m_y / scalar, 
                     m_z / scalar);
    }

    /**
      *  @brief Adds another vector to this vector.
      *  @param other Vector to add.
      *  @return Reference to this vector.
    */
    Vector3& 
    operator+=(const Vector3& other) {
      m_x += other.m_x;
      m_y += other.m_y;
      m_z += other.m_z;
      
      return *this;
    }

    /**
      *  @brief Subtracts another vector from this vector.
      *  @param other Vector to subtract.
      *  @return Reference to this vector.
    */
    Vector3& 
    operator-=(const Vector3& other) {
      m_x -= other.m_x;
      m_y -= other.m_y;
      m_z -= other.m_z;
      
      return *this;
    }

    /**
      *  @brief Multiplies this vector by another vector component-wise.
      *  @param other Vector to multiply.
      *  @return Reference to this vector.
    */
    Vector3& 
    operator*=(const Vector3& other) {
      m_x *= other.m_x;
      m_y *= other.m_y;
      m_z *= other.m_z;
      
      return *this;
    }

    /**
      *  @brief Multiplies this vector by a scalar.
      *  @param scalar Scalar value.
      *  @return Reference to this vector.
    */
    Vector3& 
    operator*=(float scalar) {
      m_x *= scalar;
      m_y *= scalar;
      m_z *= scalar;
      
      return *this;
    }

    /**
      *  @brief Checks if two vectors are equal.
      *  @param other Vector to compare.
      *  @return True if equal, false otherwise.
    */
    bool 
    operator==(const Vector3& other) const {
      return (m_x == other.m_x && 
              m_y == other.m_y && 
              m_z == other.m_z);
    }

    /**
      *  @brief Checks if two vectors are not equal.
      *  @param other Vector to compare.
      *  @return True if not equal, false otherwise.
    */
    bool 
    operator!=(const Vector3& other) const {
      return (m_x != other.m_x || 
              m_y != other.m_y || 
              m_z != other.m_z);
    }

    /**
      *  @brief Accesses vector components by index.
      *  @param index Index (0 for x, 1 for y, 2 for z).
      *  @return Reference to the component.
      *  @throws "Index out of range for Vector3" if index is invalid.
    */
    float& 
    operator[](int index) {
      switch (index) {
        case 0: return m_x;
        case 1: return m_y;
        case 2: return m_z;
				default: throw "Index out of range for Vector3";
      }
    }

    /**
      *  @brief Calculates the length of a vector.
      *  @param other Vector to calculate length for.
      *  @return Length of the vector.
    */
    static float 
    length(const Vector3& other) {
      return EngineMath::sqrt((other.m_x * other.m_x) + 
                              (other.m_y * other.m_y) + 
                              (other.m_z * other.m_z));
    }

    /**
      *  @brief Calculates the squared length of a vector.
      *  @param other Vector to calculate squared length for.
      *  @return Squared length of the vector.
    */
    static float 
    squaredLength(const Vector3& other) {
      return (other.m_x * other.m_x) + 
             (other.m_y * other.m_y) + 
             (other.m_z * other.m_z);
    }

    /**
      *  @brief Calculates the dot product of two vectors.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Dot product.
    */
    static float 
    dotProduct(const Vector3& a, const Vector3& b) {
      return (a.m_x * b.m_x) + 
             (a.m_y * b.m_y) + 
             (a.m_z * b.m_z);
    }

    /**
      *  @brief Calculates the cross product of two vectors.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Cross product vector.
    */
    Vector3 
    crossProduct(const Vector3& a, const Vector3& b) const {
      return Vector3(
        a.m_y * b.m_z - a.m_z * b.m_y,
        a.m_z * b.m_x - a.m_x * b.m_z,
        a.m_x * b.m_y - a.m_y * b.m_x
      );
    }

    /**
      *  @brief Normalizes a vector.
      *  @param other Vector to normalize.
      *  @return Normalized vector.
    */
    Vector3 
    normalize(const Vector3& other) const {
      float len = other.length(other);
      
      return Vector3(other.m_x / len, 
                     other.m_y / len, 
                     other.m_z / len);
    }

    /**
      *  @brief Normalizes a vector (non-const version).
      *  @param other Vector to normalize.
      *  @return Normalized vector.
    */
    Vector3 
    normalize(Vector3& other) {
      float len = other.length(other);
      
      return Vector3(other.m_x / len, 
                     other.m_y / len, 
                     other.m_z / len);
    }

    /**
      *  @brief Normalizes this vector in place.
    */
    void 
    normalizeInPlace() {
      float len = length(*this);
      m_x /= len;
      m_y /= len;
      m_z /= len;
    }

    /**
      *  @brief Calculates the distance between two vectors.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Distance.
    */
    static float 
    distance(const Vector3& a, const Vector3& b) {
      return (b - a).length(b - a);
    }

    /**
      *  @brief Linearly interpolates between two vectors.
      *  @param a Start vector.
      *  @param b End vector.
      *  @param t Interpolation factor.
      *  @return Interpolated vector.
    */
    static Vector3 
    lerp(const Vector3& a, const Vector3& b, float t) {
      return a + (b - a) * t;
    }

    /**
      *  @brief Returns a zero vector.
      *  @return Zero vector.
    */
    static Vector3 
    zero() { 
      return Vector3(0.0f, 0.0f, 0.0f); 
    }

    /**
      *  @brief Returns a unit vector (all components are 1).
      *  @return Unit vector.
    */
    static Vector3 
    unit() { 
      return Vector3(1.0f, 1.0f, 1.0f); 
    }

    /**
      *  @brief Spherically interpolates between two vectors.
      *  @param a Start vector.
      *  @param b End vector.
      *  @param t Interpolation factor.
      *  @return Interpolated vector.
    */
    static Vector3 
    slerp(const Vector3& a, const Vector3& b, float t) {
      float dot = dotProduct(a, b);
      dot = EngineMath::clamp(dot, -1.0f, 1.0f);
      float theta = EngineMath::aCos(dot) * t;

      Vector3 relative = b - a * dot;
      float relLen = Vector3::length(relative);
      relative = relative / relLen;

      return a * EngineMath::cos(theta) + relative * EngineMath::sin(theta);
    }

    /**
      *  @brief Reflects a vector off a surface with the given normal.
      *  @param direction Incident vector.
      *  @param normal Surface normal.
      *  @return Reflected vector.
    */
    static Vector3 
    reflect(const Vector3& direction, const Vector3& normal) {
			float dot = dotProduct(direction, normal);
      
      return direction - normal * (2.0f * dot);
    }

    /**
      *  @brief Projects vector a onto vector b.
      *  @param a Vector to project.
      *  @param b Vector to project onto.
      *  @return Projected vector.
    */
    static Vector3 
    project(const Vector3& a, const Vector3& b) {
			float dot = dotProduct(a, b);
      float lenSq = squaredLength(b);
      
      return b * (dot / lenSq);
    }

    /**
      *  @brief Clamps a vector's components between min and max vectors.
      *  @param value Vector to clamp.
      *  @param min Minimum vector.
      *  @param max Maximum vector.
      *  @return Clamped vector.
    */
    static Vector3 
    clamp(const Vector3& value, const Vector3& min, const Vector3& max) {
      float cx = EngineMath::clamp(value.m_x, min.m_x, max.m_x);
      float cy = EngineMath::clamp(value.m_y, min.m_y, max.m_y);
      float cz = EngineMath::clamp(value.m_z, min.m_z, max.m_z);
      
      return Vector3(cx, cy, cz);
    }

    /**
      *  @brief Calculates the angle between two vectors.
      *  @param a First vector.
      *  @param b Second vector.
      *  @return Angle in radians.
    */
    static float 
    angle(const Vector3& a, const Vector3& b) {
			float dot = dotProduct(a, b);
      float lenA = length(a);
      float lenB = length(b);
      
      return EngineMath::aCos(dot / (lenA * lenB));
    }
  };
}
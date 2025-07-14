#pragma once

#include "../Vectors/Vector2.h"
#include "../Utilities/EngineMath.h"

/**
  * @file Matrix3x3.h
  * @brief Defines the Matrix3x3 class for 2D transformations using 3x3 matrices.
*/

namespace EngineUtilities {

  /**
    * @class Matrix3x3
    * @brief Represents a 3x3 matrix for 2D transformations (rotation, scaling, translation).
  */
  class Matrix3x3 {
  public:
    /**
      * @brief Matrix elements in row-major order.
    */
    float m_m00, m_m01, m_m02;
    float m_m10, m_m11, m_m12;
    float m_m20, m21, m22;

    /**
      * @brief Default constructor. Initializes as identity matrix.
    */
    Matrix3x3()
      : m_m00(1.0f), m_m01(0.0f), m_m02(0.0f),
        m_m10(0.0f), m_m11(1.0f), m_m12(0.0f),
        m_m20(0.0f), m21(0.0f), m22(1.0f) {
    }

    /**
      * @brief Constructor with all 9 components.
      * @param a00 Row 0, Col 0
      * @param a01 Row 0, Col 1
      * @param a02 Row 0, Col 2
      * @param a10 Row 1, Col 0
      * @param a11 Row 1, Col 1
      * @param a12 Row 1, Col 2
      * @param a20 Row 2, Col 0
      * @param a21 Row 2, Col 1
      * @param a22 Row 2, Col 2
    */
    Matrix3x3(float a00, float a01, float a02,
              float a10, float a11, float a12,
              float a20, float a21, float a22)
      : m_m00(a00), m_m01(a01), m_m02(a02),
        m_m10(a10), m_m11(a11), m_m12(a12),
        m_m20(a20), m21(a21), m22(a22) {}

    /**
      * @brief Returns the identity matrix.
      * @return Identity matrix.
    */
    static Matrix3x3 
    Identity() {
      return Matrix3x3(
        1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a translation matrix.
      * @param tx Translation in x.
      * @param ty Translation in y.
      * @return Translation matrix.
    */
    static Matrix3x3 
    Translation(float tx, float ty) {
      return Matrix3x3(
        1.0f, 0.0f, tx,
        0.0f, 1.0f, ty,
        0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a rotation matrix.
      * @param radians Angle in radians.
      * @return Rotation matrix.
    */
    static Matrix3x3 
    Rotation(float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      
      return Matrix3x3(
        c, -s, 0.0f,
        s,  c, 0.0f,
        0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a scaling matrix.
      * @param sx Scale in x.
      * @param sy Scale in y.
      * @return Scaling matrix.
    */
    static Matrix3x3 
    Scale(float sx, float sy) {
      return Matrix3x3(
        sx, 0.0f, 0.0f,
        0.0f, sy, 0.0f,
        0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Returns the transpose of the matrix.
      * @return Transposed matrix.
    */
    Matrix3x3 
    Transpose() const {
      return Matrix3x3(
        m_m00, m_m10, m_m20,
        m_m01, m_m11, m21,
        m_m02, m_m12, m22
      );
    }

    /**
      * @brief Returns the inverse of the matrix.
      * @return Inverse matrix, or identity if not invertible.
    */
    Matrix3x3 
    Inverse() const {
      float det = Determinant();
      if (EngineMath::abs(det) < EngineMath::epsilon) {
        return Matrix3x3(); // Return identity if not invertible
      }

      float invDet = 1.0f / det;

      // Calculate cofactors for 3x3 matrix
      float c00 =  m_m11 * m22 - m_m12 * m21;
      float c01 = -(m_m10 * m22 - m_m12 * m_m20);
      float c02 =  m_m10 * m21 - m_m11 * m_m20;
      float c10 = -(m_m01 * m22 - m_m02 * m21);
      float c11 =  m_m00 * m22 - m_m02 * m_m20;
      float c12 = -(m_m00 * m21 - m_m01 * m_m20);
      float c20 =  m_m01 * m_m12 - m_m02 * m_m11;
      float c21 = -(m_m00 * m_m12 - m_m02 * m_m10);
      float c22 =  m_m00 * m_m11 - m_m01 * m_m10;

      return Matrix3x3(
        c00 * invDet, c10 * invDet, c20 * invDet,
        c01 * invDet, c11 * invDet, c21 * invDet,
        c02 * invDet, c12 * invDet, c22 * invDet
      );
    }

    /**
      * @brief Calculates the determinant of the matrix.
      * @return Determinant value.
    */
    float 
    Determinant() const {
      return
        m_m00 * (m_m11 * m22 - m_m12 * m21) -
        m_m01 * (m_m10 * m22 - m_m12 * m_m20) +
        m_m02 * (m_m10 * m21 - m_m11 * m_m20);
    }

    /**
      * @brief Adds another matrix to this matrix.
      * @param other Matrix to add.
      * @return Sum of matrices.
    */
    Matrix3x3 
    operator+(const Matrix3x3& other) const {
      return Matrix3x3(
        m_m00 + other.m_m00, m_m01 + other.m_m01, m_m02 + other.m_m02,
        m_m10 + other.m_m10, m_m11 + other.m_m11, m_m12 + other.m_m12,
        m_m20 + other.m_m20, m21 + other.m21, m22 + other.m22
      );
    }

    /**
      * @brief Subtracts another matrix from this matrix.
      * @param other Matrix to subtract.
      * @return Difference of matrices.
    */
    Matrix3x3 
    operator-(const Matrix3x3& other) const {
      return Matrix3x3(
        m_m00 - other.m_m00, m_m01 - other.m_m01, m_m02 - other.m_m02,
        m_m10 - other.m_m10, m_m11 - other.m_m11, m_m12 - other.m_m12,
        m_m20 - other.m_m20, m21 - other.m21, m22 - other.m22
      );
    }

    /**
      * @brief Multiplies this matrix by another matrix.
      * @param other Matrix to multiply.
      * @return Product of matrices.
    */
    Matrix3x3 
    operator*(const Matrix3x3& other) const {
      return Matrix3x3(
        m_m00 * other.m_m00 + m_m01 * other.m_m10 + m_m02 * other.m_m20,
        m_m00 * other.m_m01 + m_m01 * other.m_m11 + m_m02 * other.m21,
        m_m00 * other.m_m02 + m_m01 * other.m_m12 + m_m02 * other.m22,

        m_m10 * other.m_m00 + m_m11 * other.m_m10 + m_m12 * other.m_m20,
        m_m10 * other.m_m01 + m_m11 * other.m_m11 + m_m12 * other.m21,
        m_m10 * other.m_m02 + m_m11 * other.m_m12 + m_m12 * other.m22,

        m_m20 * other.m_m00 + m21 * other.m_m10 + m22 * other.m_m20,
        m_m20 * other.m_m01 + m21 * other.m_m11 + m22 * other.m21,
        m_m20 * other.m_m02 + m21 * other.m_m12 + m22 * other.m22
      );
    }

    /**
      * @brief Multiplies this matrix by a Vector2 (assumes homogeneous coordinates).
      * @param v Vector2 to transform.
      * @return Transformed Vector2.
    */
    Vector2 
    operator*(const Vector2& v) const {
      float x = m_m00 * v.m_x + m_m01 * v.m_y + m_m02;
      float y = m_m10 * v.m_x + m_m11 * v.m_y + m_m12;
      
      return Vector2(x, y);
    }

    /**
      * @brief Multiplies this matrix by a scalar.
      * @param scalar Scalar value.
      * @return Scaled matrix.
    */
    Matrix3x3 
    operator*(float scalar) const {
      return Matrix3x3(
        m_m00 * scalar, m_m01 * scalar, m_m02 * scalar,
        m_m10 * scalar, m_m11 * scalar, m_m12 * scalar,
        m_m20 * scalar, m21 * scalar, m22 * scalar
      );
    }

    /**
      * @brief Compares this matrix with another for equality.
      * @param other Matrix to compare.
      * @return True if all elements are equal within epsilon.
    */
    bool 
    operator==(const Matrix3x3& other) const {
      return EngineMath::abs(m_m00 - other.m_m00) < EngineMath::epsilon &&
             EngineMath::abs(m_m01 - other.m_m01) < EngineMath::epsilon &&
             EngineMath::abs(m_m02 - other.m_m02) < EngineMath::epsilon &&
             EngineMath::abs(m_m10 - other.m_m10) < EngineMath::epsilon &&
             EngineMath::abs(m_m11 - other.m_m11) < EngineMath::epsilon &&
             EngineMath::abs(m_m12 - other.m_m12) < EngineMath::epsilon &&
             EngineMath::abs(m_m20 - other.m_m20) < EngineMath::epsilon &&
             EngineMath::abs(m21 - other.m21) < EngineMath::epsilon &&
             EngineMath::abs(m22 - other.m22) < EngineMath::epsilon;
    }

    /**
      * @brief Compares this matrix with another for inequality.
      * @param other Matrix to compare.
      * @return True if any element is different.
    */
    bool 
    operator!=(const Matrix3x3& other) const {
      return !(*this == other);
    }
  };
}
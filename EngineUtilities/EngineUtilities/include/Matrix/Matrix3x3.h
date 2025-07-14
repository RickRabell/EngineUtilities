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
    float m00, m01, m02;
    float m10, m11, m12;
    float m20, m21, m22;

    /**
     * @brief Default constructor. Initializes as identity matrix.
     */
    Matrix3x3()
      : m00(1.0f), m01(0.0f), m02(0.0f),
        m10(0.0f), m11(1.0f), m12(0.0f),
        m20(0.0f), m21(0.0f), m22(1.0f) {
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
      : m00(a00), m01(a01), m02(a02),
        m10(a10), m11(a11), m12(a12),
        m20(a20), m21(a21), m22(a22) {}

    /**
     * @brief Returns the identity matrix.
     * @return Identity matrix.
     */
    static Matrix3x3 Identity() {
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
    static Matrix3x3 Translation(float tx, float ty) {
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
    static Matrix3x3 Rotation(float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      return Matrix3x3(
        c, -s, 0.0f,
        s,  c, 0.0f,
        0.0f, 0.0f, 1.0f
      );
    }
    
    // FromEuler

    /**
     * @brief Creates a scaling matrix.
     * @param sx Scale in x.
     * @param sy Scale in y.
     * @return Scaling matrix.
     */
    static Matrix3x3 Scale(float sx, float sy) {
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
    Matrix3x3 Transpose() const {
      return Matrix3x3(
        m00, m10, m20,
        m01, m11, m21,
        m02, m12, m22
      );
    }

    /**
     * @brief Returns the inverse of the matrix.
     * @return Inverse matrix, or identity if not invertible.
     */
    Matrix3x3 Inverse() const {
      float det = Determinant();
      if (EngineMath::abs(det) < EngineMath::epsilon)
        return Matrix3x3(); // Return identity if not invertible

      float invDet = 1.0f / det;

      // Calculate cofactors for 3x3 matrix
      float c00 =  m11 * m22 - m12 * m21;
      float c01 = -(m10 * m22 - m12 * m20);
      float c02 =  m10 * m21 - m11 * m20;
      float c10 = -(m01 * m22 - m02 * m21);
      float c11 =  m00 * m22 - m02 * m20;
      float c12 = -(m00 * m21 - m01 * m20);
      float c20 =  m01 * m12 - m02 * m11;
      float c21 = -(m00 * m12 - m02 * m10);
      float c22 =  m00 * m11 - m01 * m10;

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
    float Determinant() const {
      return
        m00 * (m11 * m22 - m12 * m21) -
        m01 * (m10 * m22 - m12 * m20) +
        m02 * (m10 * m21 - m11 * m20);
    }

    /**
     * @brief Adds another matrix to this matrix.
     * @param other Matrix to add.
     * @return Sum of matrices.
     */
    Matrix3x3 operator+(const Matrix3x3& other) const {
      return Matrix3x3(
        m00 + other.m00, m01 + other.m01, m02 + other.m02,
        m10 + other.m10, m11 + other.m11, m12 + other.m12,
        m20 + other.m20, m21 + other.m21, m22 + other.m22
      );
    }

    /**
     * @brief Subtracts another matrix from this matrix.
     * @param other Matrix to subtract.
     * @return Difference of matrices.
     */
    Matrix3x3 operator-(const Matrix3x3& other) const {
      return Matrix3x3(
        m00 - other.m00, m01 - other.m01, m02 - other.m02,
        m10 - other.m10, m11 - other.m11, m12 - other.m12,
        m20 - other.m20, m21 - other.m21, m22 - other.m22
      );
    }

    /**
     * @brief Multiplies this matrix by another matrix.
     * @param other Matrix to multiply.
     * @return Product of matrices.
     */
    Matrix3x3 operator*(const Matrix3x3& other) const {
      return Matrix3x3(
        m00 * other.m00 + m01 * other.m10 + m02 * other.m20,
        m00 * other.m01 + m01 * other.m11 + m02 * other.m21,
        m00 * other.m02 + m01 * other.m12 + m02 * other.m22,

        m10 * other.m00 + m11 * other.m10 + m12 * other.m20,
        m10 * other.m01 + m11 * other.m11 + m12 * other.m21,
        m10 * other.m02 + m11 * other.m12 + m12 * other.m22,

        m20 * other.m00 + m21 * other.m10 + m22 * other.m20,
        m20 * other.m01 + m21 * other.m11 + m22 * other.m21,
        m20 * other.m02 + m21 * other.m12 + m22 * other.m22
      );
    }

    /**
     * @brief Multiplies this matrix by a Vector2 (assumes homogeneous coordinates).
     * @param v Vector2 to transform.
     * @return Transformed Vector2.
     */
    Vector2 operator*(const Vector2& v) const {
      float x = m00 * v.m_x + m01 * v.m_y + m02;
      float y = m10 * v.m_x + m11 * v.m_y + m12;
      return Vector2(x, y);
    }

    /**
     * @brief Multiplies this matrix by a scalar.
     * @param scalar Scalar value.
     * @return Scaled matrix.
     */
    Matrix3x3 operator*(float scalar) const {
      return Matrix3x3(
        m00 * scalar, m01 * scalar, m02 * scalar,
        m10 * scalar, m11 * scalar, m12 * scalar,
        m20 * scalar, m21 * scalar, m22 * scalar
      );
    }

    // operator/ para multiplicar con un escalar

    /**
     * @brief Compares this matrix with another for equality.
     * @param other Matrix to compare.
     * @return True if all elements are equal within epsilon.
     */
    bool operator==(const Matrix3x3& other) const {
      return EngineMath::abs(m00 - other.m00) < EngineMath::epsilon &&
             EngineMath::abs(m01 - other.m01) < EngineMath::epsilon &&
             EngineMath::abs(m02 - other.m02) < EngineMath::epsilon &&
             EngineMath::abs(m10 - other.m10) < EngineMath::epsilon &&
             EngineMath::abs(m11 - other.m11) < EngineMath::epsilon &&
             EngineMath::abs(m12 - other.m12) < EngineMath::epsilon &&
             EngineMath::abs(m20 - other.m20) < EngineMath::epsilon &&
             EngineMath::abs(m21 - other.m21) < EngineMath::epsilon &&
             EngineMath::abs(m22 - other.m22) < EngineMath::epsilon;
    }

    /**
     * @brief Compares this matrix with another for inequality.
     * @param other Matrix to compare.
     * @return True if any element is different.
     */
    bool operator!=(const Matrix3x3& other) const {
      return !(*this == other);
    }
  };

}
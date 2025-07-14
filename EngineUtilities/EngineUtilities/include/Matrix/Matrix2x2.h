#pragma once

#include "../Vectors/Vector2.h"
#include "../Utilities/EngineMath.h"

/*!
 * @namespace EngineUtilities
 * @brief Contains utility classes and functions for the engine.
 */
namespace EngineUtilities {

  /*!
   * @class Matrix2x2
   * @brief Represents a 2x2 matrix, used for 2D transformations like rotation and scaling.
   *
   * The matrix is stored in row-major order:
   * | m00 m01 |
   * | m10 m11 |
   */
  class Matrix2x2 {
  public:
    /// @brief The element at row 0, column 0.
    float m00;
    /// @brief The element at row 0, column 1.
    float m01;
    /// @brief The element at row 1, column 0.
    float m10;
    /// @brief The element at row 1, column 1.
    float m11;

    /*!
     * @brief Default constructor. Initializes the matrix to an identity matrix.
     */
    Matrix2x2()
      : m00(1.0f), m01(0.0f), m10(0.0f), m11(1.0f) {
    }

    /*!
     * @brief Constructs a matrix from four float components.
     * @param a00 The element at row 0, column 0.
     * @param a01 The element at row 0, column 1.
     * @param a10 The element at row 1, column 0.
     * @param a11 The element at row 1, column 1.
     */
    Matrix2x2(float a00, float a01, float a10, float a11)
      : m00(a00), m01(a01), m10(a10), m11(a11) {
    }

    /*!
     * @brief Creates and returns a 2x2 identity matrix.
     * @return The 2x2 identity matrix.
     */
    static Matrix2x2 Identity() {
      return Matrix2x2(1.0f, 0.0f, 0.0f, 1.0f);
    }

    // Traslación (no aplica en 2x2, normalmente es para 3x3 o 4x4, se omite)

    /*!
     * @brief Creates a 2D rotation matrix from an angle in radians.
     * @param radians The angle of rotation in radians.
     * @return The corresponding rotation matrix.
     */
    static Matrix2x2 Rotation(float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      return Matrix2x2(c, -s, s, c);
    }

    /*!
     * @brief Creates a 2D scaling matrix.
     * @param sx The scaling factor along the x-axis.
     * @param sy The scaling factor along the y-axis.
     * @return The corresponding scaling matrix.
     */
    static Matrix2x2 Scale(float sx, float sy) {
      return Matrix2x2(sx, 0.0f, 0.0f, sy);
    }

    /*!
     * @brief Computes the transpose of this matrix.
     * @return The transposed matrix.
     */
    Matrix2x2 Transpose() const {
      return Matrix2x2(m00, m10, m01, m11);
    }

    /*!
     * @brief Computes the inverse of this matrix.
     * @return The inverse of the matrix. If the matrix is not invertible (determinant is close to zero), it returns an identity matrix.
     */
    Matrix2x2 Inverse() const {
      float det = Determinant();
      if (EngineMath::abs(det) < EngineMath::epsilon)
        return Matrix2x2(); // Return identity if not invertible
      float invDet = 1.0f / det;
      return Matrix2x2(
        m11 * invDet, -m01 * invDet,
        -m10 * invDet, m00 * invDet
      );
    }

    /*!
     * @brief Computes the determinant of this matrix.
     * @return The determinant.
     */
    float Determinant() const {
      return m00 * m11 - m01 * m10;
    }

    /*!
     * @brief Adds another matrix to this matrix.
     * @param other The matrix to add.
     * @return The result of the addition.
     */
    Matrix2x2 operator+(const Matrix2x2& other) const {
      return Matrix2x2(
        m00 + other.m00, m01 + other.m01,
        m10 + other.m10, m11 + other.m11
      );
    }

    /*!
     * @brief Subtracts another matrix from this matrix.
     * @param other The matrix to subtract.
     * @return The result of the subtraction.
     */
    Matrix2x2 operator-(const Matrix2x2& other) const {
      return Matrix2x2(
        m00 - other.m00, m01 - other.m01,
        m10 - other.m10, m11 - other.m11
      );
    }

    /*!
     * @brief Multiplies this matrix by another matrix.
     * @param other The matrix to multiply by.
     * @return The result of the multiplication.
     */
    Matrix2x2 operator*(const Matrix2x2& other) const {
      return Matrix2x2(
        m00 * other.m00 + m01 * other.m10, m00 * other.m01 + m01 * other.m11,
        m10 * other.m00 + m11 * other.m10, m10 * other.m01 + m11 * other.m11
      );
    }

    /*!
     * @brief Multiplies this matrix by a 2D vector, applying the transformation.
     * @param v The vector to transform.
     * @return The transformed vector.
     */
    Vector2 operator*(const Vector2& v) const {
      return Vector2(
        m00 * v.m_x + m01 * v.m_y,
        m10 * v.m_x + m11 * v.m_y
      );
    }

    /*!
     * @brief Multiplies this matrix by a scalar.
     * @param scalar The scalar value to multiply with.
     * @return The result of the scalar multiplication.
     */
    Matrix2x2 operator*(float scalar) const {
      return Matrix2x2(
        m00 * scalar, m01 * scalar,
        m10 * scalar, m11 * scalar
      );
    }

    /*!
     * @brief Compares this matrix with another for equality.
     * @param other The matrix to compare with.
     * @return True if the matrices are equal within a small tolerance (epsilon), false otherwise.
     */
    bool operator==(const Matrix2x2& other) const {
      return EngineMath::abs(m00 - other.m00) < EngineMath::epsilon &&
        EngineMath::abs(m01 - other.m01) < EngineMath::epsilon &&
        EngineMath::abs(m10 - other.m10) < EngineMath::epsilon &&
        EngineMath::abs(m11 - other.m11) < EngineMath::epsilon;
    }

    /*!
     * @brief Compares this matrix with another for inequality.
     * @param other The matrix to compare with.
     * @return True if the matrices are not equal, false otherwise.
     */
    bool operator!=(const Matrix2x2& other) const {
      return !(*this == other);
    }
  };

}
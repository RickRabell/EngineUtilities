#pragma once

#include "../Vectors/Vector3.h"
#include "../Vectors/Vector4.h"
#include "../Utilities/EngineMath.h"

// Forward declaration for Quaternion to avoid circular dependency
namespace EngineUtilities {
  class Quaternion;
}

/**
  * @file Matrix4x4.h
  * @brief Defines the Matrix4x4 class for 3D transformations and projections.
*/
namespace EngineUtilities {

  /**
    * @class Matrix4x4
    * @brief Represents a 4x4 matrix for 3D transformations (rotation, scaling, translation) and projections.
  */
  class Matrix4x4 {
  public:
    /**
      * @brief Matrix elements in row-major order.
    */
    float m_m00, m_m01, m_m02, m_m03;
    float m_m10, m_m11, m_m12, m_m13;
    float m_m20, m_m21, m_m22, m_m23;
    float m_m30, m_m31, m_m32, m_m33;

    /**
      * @brief Default constructor. Initializes as identity matrix.
    */
    Matrix4x4()
      : m_m00(1.0f), m_m01(0.0f), m_m02(0.0f), m_m03(0.0f),
      m_m10(0.0f), m_m11(1.0f), m_m12(0.0f), m_m13(0.0f),
      m_m20(0.0f), m_m21(0.0f), m_m22(1.0f), m_m23(0.0f),
      m_m30(0.0f), m_m31(0.0f), m_m32(0.0f), m_m33(1.0f) {
    }

    /**
      * @brief Constructor with all 16 components.
    */
    Matrix4x4(float a00, float a01, float a02, float a03,
              float a10, float a11, float a12, float a13,
              float a20, float a21, float a22, float a23,
              float a30, float a31, float a32, float a33)
      : m_m00(a00), m_m01(a01), m_m02(a02), m_m03(a03),
      m_m10(a10), m_m11(a11), m_m12(a12), m_m13(a13),
      m_m20(a20), m_m21(a21), m_m22(a22), m_m23(a23),
      m_m30(a30), m_m31(a31), m_m32(a32), m_m33(a33) {
    }

    /**
      * @brief Returns the identity matrix.
      * @return Identity matrix.
    */
    static Matrix4x4 
    Identity() {
      return Matrix4x4(
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a translation matrix.
      * @param v Vector representing the translation.
      * @return Translation matrix.
    */
    static Matrix4x4 
    Translation(const Vector3& v) {
      return Matrix4x4(
        1.0f, 0.0f, 0.0f, v.m_x,
        0.0f, 1.0f, 0.0f, v.m_y,
        0.0f, 0.0f, 1.0f, v.m_z,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a rotation matrix around an arbitrary axis.
      * @param axis The axis to rotate around (should be normalized).
      * @param radians Angle in radians.
      * @return Rotation matrix.
    */
    static Matrix4x4 
    RotateAxis(const Vector3& axis, float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      float t = 1.0f - c;
      float x = axis.m_x;
      float y = axis.m_y;
      float z = axis.m_z;

      return Matrix4x4(
        t * x * x + c, t * x * y - s * z, t * x * z + s * y, 0.0f,
        t * x * y + s * z, t * y * y + c, t * y * z - s * x, 0.0f,
        t * x * z - s * y, t * y * z + s * x, t * z * z + c, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a rotation matrix around the X axis.
      * @param radians Angle in radians.
      * @return Rotation matrix.
    */
    static Matrix4x4 
    RotateX(float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      
      return Matrix4x4(
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, c, -s, 0.0f,
        0.0f, s, c, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a rotation matrix around the Y axis.
      * @param radians Angle in radians.
      * @return Rotation matrix.
    */
    static Matrix4x4 
    RotateY(float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      
      return Matrix4x4(
        c, 0.0f, s, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        -s, 0.0f, c, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a rotation matrix around the Z axis.
      * @param radians Angle in radians.
      * @return Rotation matrix.
    */
    static Matrix4x4 
    RotateZ(float radians) {
      float c = EngineMath::cos(radians);
      float s = EngineMath::sin(radians);
      
      return Matrix4x4(
        c, -s, 0.0f, 0.0f,
        s, c, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Creates a scaling matrix.
      * @param v Vector representing the scale on each axis.
      * @return Scaling matrix.
    */
    static Matrix4x4 
    Scale(const Vector3& v) {
      return Matrix4x4(
        v.m_x, 0.0f, 0.0f, 0.0f,
        0.0f, v.m_y, 0.0f, 0.0f,
        0.0f, 0.0f, v.m_z, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /**
      * @brief Returns the transpose of the matrix.
      * @return Transposed matrix.
    */
    Matrix4x4 
    Transpose() const {
      return Matrix4x4(
        m_m00, m_m10, m_m20, m_m30,
        m_m01, m_m11, m_m21, m_m31,
        m_m02, m_m12, m_m22, m_m32,
        m_m03, m_m13, m_m23, m_m33
      );
    }

    /**
      * @brief Calculates the determinant of the matrix.
      * @return Determinant value.
    */
    float 
    Determinant() const {
      return
        m_m00 * (m_m11 * (m_m22 * m_m33 - m_m23 * m_m32) - m_m12 * (m_m21 * m_m33 - m_m23 * m_m31) + m_m13 * (m_m21 * m_m32 - m_m22 * m_m31)) -
        m_m01 * (m_m10 * (m_m22 * m_m33 - m_m23 * m_m32) - m_m12 * (m_m20 * m_m33 - m_m23 * m_m30) + m_m13 * (m_m20 * m_m32 - m_m22 * m_m30)) +
        m_m02 * (m_m10 * (m_m21 * m_m33 - m_m23 * m_m31) - m_m11 * (m_m20 * m_m33 - m_m23 * m_m30) + m_m13 * (m_m20 * m_m31 - m_m21 * m_m30)) -
        m_m03 * (m_m10 * (m_m21 * m_m32 - m_m22 * m_m31) - m_m11 * (m_m20 * m_m32 - m_m22 * m_m30) + m_m12 * (m_m20 * m_m31 - m_m21 * m_m30));
    }

    /**
      * @brief Returns the inverse of the matrix.
      * @return Inverse matrix, or identity if not invertible.
    */
    Matrix4x4 
    Inverse() const {
        float det = Determinant();
        if (EngineMath::absf(det) < EngineMath::epsilon) {
            return Matrix4x4::Identity(); // Not invertible
        }

        float invDet = 1.0f / det;
        Matrix4x4 result;

        // Compute the matrix of cofactors, then transpose (adjugate), then multiply by 1/det
        result.m_m00 =  (m_m11 * (m_m22 * m_m33 - m_m23 * m_m32) - 
                         m_m12 * (m_m21 * m_m33 - m_m23 * m_m31) + 
                         m_m13 * (m_m21 * m_m32 - m_m22 * m_m31)) 
                         * invDet;
        result.m_m01 = -(m_m01 * (m_m22 * m_m33 - m_m23 * m_m32) - 
                         m_m02 * (m_m21 * m_m33 - m_m23 * m_m31) + 
                         m_m03 * (m_m21 * m_m32 - m_m22 * m_m31)) 
                         * invDet;
        result.m_m02 =  (m_m01 * (m_m12 * m_m33 - m_m13 * m_m32) - 
                         m_m02 * (m_m11 * m_m33 - m_m13 * m_m31) + 
                         m_m03 * (m_m11 * m_m32 - m_m12 * m_m31)) 
                         * invDet;
        result.m_m03 = -(m_m01 * (m_m12 * m_m23 - m_m13 * m_m22) - 
                         m_m02 * (m_m11 * m_m23 - m_m13 * m_m21) + 
                         m_m03 * (m_m11 * m_m22 - m_m12 * m_m21)) 
                         * invDet;

        result.m_m10 = -(m_m10 * (m_m22 * m_m33 - m_m23 * m_m32) - 
                         m_m12 * (m_m20 * m_m33 - m_m23 * m_m30) + 
                         m_m13 * (m_m20 * m_m32 - m_m22 * m_m30)) 
                         * invDet;
        result.m_m11 =  (m_m00 * (m_m22 * m_m33 - m_m23 * m_m32) - 
                         m_m02 * (m_m20 * m_m33 - m_m23 * m_m30) + 
                         m_m03 * (m_m20 * m_m32 - m_m22 * m_m30)) 
                         * invDet;
        result.m_m12 = -(m_m00 * (m_m12 * m_m33 - m_m13 * m_m32) - 
                         m_m02 * (m_m10 * m_m33 - m_m13 * m_m30) + 
                         m_m03 * (m_m10 * m_m32 - m_m12 * m_m30)) 
                         * invDet;
        result.m_m13 =  (m_m00 * (m_m12 * m_m23 - m_m13 * m_m22) - 
                         m_m02 * (m_m10 * m_m23 - m_m13 * m_m20) + 
                         m_m03 * (m_m10 * m_m22 - m_m12 * m_m20)) 
                         * invDet;

        result.m_m20 =  (m_m10 * (m_m21 * m_m33 - m_m23 * m_m31) - 
                         m_m11 * (m_m20 * m_m33 - m_m23 * m_m30) + 
                         m_m13 * (m_m20 * m_m31 - m_m21 * m_m30)) 
                         * invDet;
        result.m_m21 = -(m_m00 * (m_m21 * m_m33 - m_m23 * m_m31) - 
                         m_m01 * (m_m20 * m_m33 - m_m23 * m_m30) + 
                         m_m03 * (m_m20 * m_m31 - m_m21 * m_m30)) 
                         * invDet;
        result.m_m22 =  (m_m00 * (m_m11 * m_m33 - m_m13 * m_m31) - 
                         m_m01 * (m_m10 * m_m33 - m_m13 * m_m30) + 
                         m_m03 * (m_m10 * m_m31 - m_m11 * m_m30)) 
                         * invDet;
        result.m_m23 = -(m_m00 * (m_m11 * m_m23 - m_m13 * m_m21) - 
                         m_m01 * (m_m10 * m_m23 - m_m13 * m_m20) + 
                         m_m03 * (m_m10 * m_m21 - m_m11 * m_m20)) 
                         * invDet;

        result.m_m30 = -(m_m10 * (m_m21 * m_m32 - m_m22 * m_m31) - 
                         m_m11 * (m_m20 * m_m32 - m_m22 * m_m30) + 
                         m_m12 * (m_m20 * m_m31 - m_m21 * m_m30)) 
                         * invDet;
        result.m_m31 =  (m_m00 * (m_m21 * m_m32 - m_m22 * m_m31) - 
                         m_m01 * (m_m20 * m_m32 - m_m22 * m_m30) + 
                         m_m02 * (m_m20 * m_m31 - m_m21 * m_m30)) 
                        * invDet;
        result.m_m32 = -(m_m00 * (m_m11 * m_m32 - m_m12 * m_m31) - 
                         m_m01 * (m_m10 * m_m32 - m_m12 * m_m30) + 
                         m_m02 * (m_m10 * m_m31 - m_m11 * m_m30)) 
                         * invDet;
        result.m_m33 =  (m_m00 * (m_m11 * m_m22 - m_m12 * m_m21) - 
                         m_m01 * (m_m10 * m_m22 - m_m12 * m_m20) + 
                         m_m02 * (m_m10 * m_m21 - m_m11 * m_m20)) 
                         * invDet;

        return result;
    }

    /**
      * @brief Creates a perspective projection matrix.
      * @param fovY Field of view in the Y direction, in radians.
      * @param aspect Aspect ratio (width / height).
      * @param zNear Near clipping plane.
      * @param zFar Far clipping plane.
      * @return Perspective projection matrix.
    */
    static Matrix4x4 
    Perspective(float fovY, float aspect, float zNear, float zFar) {
      float tanHalfFovY = EngineMath::tan(fovY / 2.0f);
      
      Matrix4x4 result = Matrix4x4(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); // Zero matrix
      result.m_m00 = 1.0f / (aspect * tanHalfFovY);
      result.m_m11 = 1.0f / tanHalfFovY;
      result.m_m22 = -(zFar + zNear) / (zFar - zNear);
      result.m_m23 = -1.0f;
      result.m_m32 = -(2.0f * zFar * zNear) / (zFar - zNear);

      return result;
    }

    /**
      * @brief Creates an orthographic projection matrix.
      * @param left Left viewing plane.
      * @param right Right viewing plane.
      * @param bottom Bottom viewing plane.
      * @param top Top viewing plane.
      * @param zNear Near clipping plane.
      * @param zFar Far clipping plane.
      * @return Orthographic projection matrix.
    */
    static Matrix4x4 
    Ortho(float left, float right, float bottom, float top, float zNear, float zFar) {
      Matrix4x4 result = Identity();
      result.m_m00 = 2.0f / (right - left);
      result.m_m11 = 2.0f / (top - bottom);
      result.m_m22 = -2.0f / (zFar - zNear);
      result.m_m03 = -(right + left) / (right - left);
      result.m_m13 = -(top + bottom) / (top - bottom);
      result.m_m23 = -(zFar + zNear) / (zFar - zNear);
      
      return result;
    }

    /**
      * @brief Creates a view matrix (LookAt).
      * @param eye The position of the camera.
      * @param target The point the camera is looking at.
      * @param up The up direction of the world (usually (0,1,0)).
      * @return View matrix.
    */
    static Matrix4x4 
    LookAt(const Vector3& eye, const Vector3& target, const Vector3& up) {
      Vector3 zaxis = (target - eye);
      zaxis.normalizeInPlace();

      Vector3 xaxis = zaxis.crossProduct(zaxis, up);
      xaxis.normalizeInPlace();
      Vector3 yaxis = xaxis.crossProduct(xaxis, zaxis);

      return Matrix4x4(
        xaxis.m_x, xaxis.m_y, xaxis.m_z, -xaxis.dotProduct(xaxis, eye),
        yaxis.m_x, yaxis.m_y, yaxis.m_z, -yaxis.dotProduct(yaxis, eye),
        -zaxis.m_x, -zaxis.m_y, -zaxis.m_z, zaxis.dotProduct(zaxis, eye),
        0.0f, 0.0f, 0.0f, 1.0f
      );
    }

    /*
      *  @brief Adds two matrices element-wise.
    */
    Matrix4x4 
    operator+(const Matrix4x4& other) const {
        return Matrix4x4(
            m_m00 + other.m_m00, m_m01 + other.m_m01, m_m02 + other.m_m02, m_m03 + other.m_m03,
            m_m10 + other.m_m10, m_m11 + other.m_m11, m_m12 + other.m_m12, m_m13 + other.m_m13,
            m_m20 + other.m_m20, m_m21 + other.m_m21, m_m22 + other.m_m22, m_m23 + other.m_m23,
            m_m30 + other.m_m30, m_m31 + other.m_m31, m_m32 + other.m_m32, m_m33 + other.m_m33
        );
    }

    /*
      *  @brief Subtracts two matrices element-wise.
    */
    Matrix4x4 
    operator-(const Matrix4x4& other) const {
        return Matrix4x4(
            m_m00 - other.m_m00, m_m01 - other.m_m01, m_m02 - other.m_m02, m_m03 - other.m_m03,
            m_m10 - other.m_m10, m_m11 - other.m_m11, m_m12 - other.m_m12, m_m13 - other.m_m13,
            m_m20 - other.m_m20, m_m21 - other.m_m21, m_m22 - other.m_m22, m_m23 - other.m_m23,
            m_m30 - other.m_m30, m_m31 - other.m_m31, m_m32 - other.m_m32, m_m33 - other.m_m33
        );
    }

    /*
      *  @brief Multiplies two matrices (matrix multiplication).
    */
    Matrix4x4 
    operator*(const Matrix4x4& other) const {
        return Matrix4x4(
            m_m00 * other.m_m00 + m_m01 * other.m_m10 + m_m02 * other.m_m20 + m_m03 * other.m_m30,
            m_m00 * other.m_m01 + m_m01 * other.m_m11 + m_m02 * other.m_m21 + m_m03 * other.m_m31,
            m_m00 * other.m_m02 + m_m01 * other.m_m12 + m_m02 * other.m_m22 + m_m03 * other.m_m32,
            m_m00 * other.m_m03 + m_m01 * other.m_m13 + m_m02 * other.m_m23 + m_m03 * other.m_m33,

            m_m10 * other.m_m00 + m_m11 * other.m_m10 + m_m12 * other.m_m20 + m_m13 * other.m_m30,
            m_m10 * other.m_m01 + m_m11 * other.m_m11 + m_m12 * other.m_m21 + m_m13 * other.m_m31,
            m_m10 * other.m_m02 + m_m11 * other.m_m12 + m_m12 * other.m_m22 + m_m13 * other.m_m32,
            m_m10 * other.m_m03 + m_m11 * other.m_m13 + m_m12 * other.m_m23 + m_m13 * other.m_m33,

            m_m20 * other.m_m00 + m_m21 * other.m_m10 + m_m22 * other.m_m20 + m_m23 * other.m_m30,
            m_m20 * other.m_m01 + m_m21 * other.m_m11 + m_m22 * other.m_m21 + m_m23 * other.m_m31,
            m_m20 * other.m_m02 + m_m21 * other.m_m12 + m_m22 * other.m_m22 + m_m23 * other.m_m32,
            m_m20 * other.m_m03 + m_m21 * other.m_m13 + m_m22 * other.m_m23 + m_m23 * other.m_m33,

            m_m30 * other.m_m00 + m_m31 * other.m_m10 + m_m32 * other.m_m20 + m_m33 * other.m_m30,
            m_m30 * other.m_m01 + m_m31 * other.m_m11 + m_m32 * other.m_m21 + m_m33 * other.m_m31,
            m_m30 * other.m_m02 + m_m31 * other.m_m12 + m_m32 * other.m_m22 + m_m33 * other.m_m32,
            m_m30 * other.m_m03 + m_m31 * other.m_m13 + m_m32 * other.m_m23 + m_m33 * other.m_m33
        );
    }

    /*
      *  @brief Multiplies the matrix by a 4D vector.
    */
    Vector4 
    operator*(const Vector4& v) const {
        return Vector4(
            m_m00 * v.m_x + m_m01 * v.m_y + m_m02 * v.m_z + m_m03 * v.m_w,
            m_m10 * v.m_x + m_m11 * v.m_y + m_m12 * v.m_z + m_m13 * v.m_w,
            m_m20 * v.m_x + m_m21 * v.m_y + m_m22 * v.m_z + m_m23 * v.m_w,
            m_m30 * v.m_x + m_m31 * v.m_y + m_m32 * v.m_z + m_m33 * v.m_w
        );
    }

    /*
      *  @brief Multiplies the matrix by a scalar.
    */
    Matrix4x4 
    operator*(float scalar) const {
        return Matrix4x4(
            m_m00 * scalar, m_m01 * scalar, m_m02 * scalar, m_m03 * scalar,
            m_m10 * scalar, m_m11 * scalar, m_m12 * scalar, m_m13 * scalar,
            m_m20 * scalar, m_m21 * scalar, m_m22 * scalar, m_m23 * scalar,
            m_m30 * scalar, m_m31 * scalar, m_m32 * scalar, m_m33 * scalar
        );
    }

    /*
      *  @brief Divides the matrix by a scalar.
    */
    Matrix4x4 
    operator/(float scalar) const {
        return Matrix4x4(
            m_m00 / scalar, m_m01 / scalar, m_m02 / scalar, m_m03 / scalar,
            m_m10 / scalar, m_m11 / scalar, m_m12 / scalar, m_m13 / scalar,
            m_m20 / scalar, m_m21 / scalar, m_m22 / scalar, m_m23 / scalar,
            m_m30 / scalar, m_m31 / scalar, m_m32 / scalar, m_m33 / scalar
        );
    }

    /*
      *  @brief Checks if two matrices are equal (element-wise comparison).
    */
    bool 
    operator==(const Matrix4x4& other) const {
        return
            m_m00 == other.m_m00 && m_m01 == other.m_m01 && m_m02 == other.m_m02 && m_m03 == other.m_m03 &&
            m_m10 == other.m_m10 && m_m11 == other.m_m11 && m_m12 == other.m_m12 && m_m13 == other.m_m13 &&
            m_m20 == other.m_m20 && m_m21 == other.m_m21 && m_m22 == other.m_m22 && m_m23 == other.m_m23 &&
            m_m30 == other.m_m30 && m_m31 == other.m_m31 && m_m32 == other.m_m32 && m_m33 == other.m_m33;
    }

    /*
      *  @brief Checks if two matrices are not equal (element-wise comparison).
    */
    bool 
    operator!=(const Matrix4x4& other) const {
        return !(*this == other);
    }

    /**
      * @brief Transforms a 3D point (w=1 is assumed).
      * @param pt The point to transform.
      * @return The transformed point.
    */
    Vector3 
    transformPoint(const Vector3& pt) const {
      float x = pt.m_x * m_m00 + pt.m_y * m_m01 + pt.m_z * m_m02 + m_m03;
      float y = pt.m_x * m_m10 + pt.m_y * m_m11 + pt.m_z * m_m12 + m_m13;
      float z = pt.m_x * m_m20 + pt.m_y * m_m21 + pt.m_z * m_m22 + m_m23;
      float w = pt.m_x * m_m30 + pt.m_y * m_m31 + pt.m_z * m_m32 + m_m33;
      
      return Vector3(x / w, y / w, z / w);
    }

    /**
      * @brief Transforms a 3D vector/direction (w=0 is assumed).
      * @param vec The vector to transform.
      * @return The transformed vector.
    */
    Vector3 
    transformVector(const Vector3& vec) const {
      return Vector3(
        vec.m_x * m_m00 + vec.m_y * m_m01 + vec.m_z * m_m02,
        vec.m_x * m_m10 + vec.m_y * m_m11 + vec.m_z * m_m12,
        vec.m_x * m_m20 + vec.m_y * m_m21 + vec.m_z * m_m22
      );
    }

    /**
      * @brief Transforms a 3D normal.
      * @param n The normal to transform.
      * @return The transformed normal.
    */
    Vector3 
    transformNormal(const Vector3& n) const {
      Matrix4x4 inv = this->Inverse().Transpose();
      
      return inv.transformVector(n);
    }
  };
}
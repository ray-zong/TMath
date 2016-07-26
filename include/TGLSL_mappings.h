#ifndef TGLSL_MAPPINGS_H
#define TGLSL_MAPPINGS_H


#include "TMatrix.h"
#include "TMatrix_4x4.h"
#include "TQuaternion.h"
#include "TVector.h"
#include "TVector_2.h"
#include "TVector_4.h"
#include "TVector_3.h"

/// @file TGLSL_mappings.h
/// @brief 兼容GLSL(OpenGL着色语言)数据类型.
/// @addtogroup TGLSL
///
/// 为了便于利用TMath模板类，使编写代码与
/// <a href="http://www.opengl.org/documentation/glsl/">GLSL</a>相似,
/// TMath提供了一系列与GLSL的向量和矩阵相似的数据类型.

/// @brief TMath命名空间.
TMATH_BEGIN_NAMESPACE

/// @addtogroup TGLSL
/// @{

/// 2-dimensional <code>float</code> Vector.
typedef TVector<float, 2> vec2;
/// 3-dimensional <code>float</code> Vector.
typedef TVector<float, 3> vec3;
/// 4-dimensional <code>float</code> Vector.
typedef TVector<float, 4> vec4;

/// 2-dimensional <code>int</code> Vector.
typedef TVector<int, 2> vec2i;
/// 3-dimensional <code>int</code> Vector.
typedef TVector<int, 3> vec3i;
/// 4-dimensional <code>int</code> Vector.
typedef TVector<int, 4> vec4i;

/// 2x2 <code>float</code> Matrix.
typedef TMatrix<float, 2, 2> mat2;
/// 3x3 <code>float</code> Matrix.
typedef TMatrix<float, 3, 3> mat3;
/// 3x3 <code>float</code> Matrix.
typedef TMatrix<float, 4, 4> mat4;

/// 2-dimensional <code>float</code> packed Vector (VectorPacked).
typedef TVectorPacked<float, 2> vec2_packed;
/// 3-dimensional <code>float</code> packed Vector (VectorPacked).
typedef TVectorPacked<float, 3> vec3_packed;
/// 4-dimensional <code>float</code> packed Vector (VectorPacked).
typedef TVectorPacked<float, 4> vec4_packed;

/// 2-dimensional <code>int</code> packed Vector (VectorPacked).
typedef TVectorPacked<int, 2> vec2i_packed;
/// 3-dimensional <code>int</code> packed Vector (VectorPacked).
typedef TVectorPacked<int, 3> vec3i_packed;
/// 4-dimensional <code>int</code> packed Vector (VectorPacked).
typedef TVectorPacked<int, 4> vec4i_packed;

/// Float-based quaternion.  Note that this is not technically
/// a GLES type, but is included for convenience.
typedef Ray::TQuaternion<float> quat;

/// @brief 计算两个3维向量叉乘.
///
/// @param v1 Vector to multiply
/// @param v2 Vector to multiply
/// @return 3-dimensional vector that contains the result.
template<class T>
inline TVector<T, 3> cross(const TVector<T, 3>& v1, const TVector<T, 3>& v2) {
  return TVector<T, 3>::CrossProduct(v1,v2);
}

/// @brief 计算任意类型的两个n维向量的点乘.
///
/// @param v1 Vector to multiply
/// @param v2 Vector to multiply
/// @return Scalar dot product result.
template<class TV>
inline typename TV::Scalar dot(const TV& v1, const TV& v2) {
  return TV::DotProduct(v1,v2);
}

/// @brief 单位化一个任意类型的n维向量.
///
/// @param v1 Vector to normalize.
/// @return 单位化向量.
template<class TV>
inline TV normalize(const TV& v1) {
  return v1.Normalized();
}

/// @}

TMATH_END_NAMESPACE

#endif // TGLSL_MAPPINGS_H

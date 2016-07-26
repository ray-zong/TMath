#ifndef TVECTOR_H
#define TVECTOR_H

#include <string>
#include <cmath>

#include "TGlobal.h"
#include "TUtilities.h"

/// @file TVector.h Vector
/// @brief Vector class and functions.
/// @addtogroup TMath_vector

// Disable spurious warnings generated by TMATH_VECTOR_OPERATION().
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4127)  // conditional expression is constant
#if _MSC_VER >= 1900             // MSVC 2015
#pragma warning(disable : 4456)  // allow shadowing in unrolled loops
#endif                           // _MSC_VER >= 1900
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warray-bounds"
#endif

/// @cond 内部
#define TMATH_VECTOR_OPERATION(OP) TMATH_UNROLLED_LOOP(i, d, OP)
/// @endcond

/// @cond 内部
#define TMATH_VECTOR_OPERATOR(OP)           \
  {                                          \
    TVector<T, d> result;                     \
    TMATH_VECTOR_OPERATION(result[i] = OP); \
    return result;                           \
  }
/// @endcond

/// @cond 内部
#define TMATH_VECTOR_SELF_OPERATOR(OP) \
  {                                     \
    TMATH_VECTOR_OPERATION(OP);        \
    return *this;                       \
  }
/// @endcond

TMATH_BEGIN_NAMESPACE

template <class T, int d>
class TVector;

/// @cond 内部使用
template <class T, int d>
static inline T DotProductHelper(const TVector<T, d>& v1,
                                 const TVector<T, d>& v2);
template <class T>
static inline T DotProductHelper(const TVector<T, 2>& v1,
                                 const TVector<T, 2>& v2);
template <class T>
static inline T DotProductHelper(const TVector<T, 3>& v1,
                                 const TVector<T, 3>& v2);
template <class T>
static inline T DotProductHelper(const TVector<T, 4>& v1,
                                 const TVector<T, 4>& v2);
/// @endcond

/// @addtogroup 向量
/// @{

/// @class TVectorPacked "ector.h"
/// @brief 封装n维向量.
///
/// Some Vector classes are padded so that it's possible to use the data
/// structures with SIMD 指令.  This structure can be used in
/// conjunction with unpacked Vector classes to pack data
/// into flat arrays suitable for sending to a GPU (e.g vertex buffers).
///
/// <p>
/// For example, to pack (store) an unpacked to packed vector:<br>
/// <blockquote><code><pre>
/// VectorPacked<float, 3> packed;
/// Vector<float, 3> vector(3, 2, 1);
/// vector.Pack(&packed);
/// </pre></code></blockquote>
/// or<br>
/// <blockquote><code><pre>
/// Vector<float, 3> vector(3, 2, 1);
/// VectorPacked<float, 3> packed = vector;
/// </pre></code></blockquote>
/// </p>
///
/// <p>
/// To initialize a vector from a packed vector:<br>
/// <blockquote><code><pre>
/// VectorPacked<float, 3> packed = { 3, 2, 1 };
/// Vector<float, 3> vector(packed);
/// </pre></code></blockquote>
///
/// @tparam T type of VectorPacked elements.
/// @tparam d dimensions (number of elements) in the VectorPacked structure.
template <class T, int d>
struct TVectorPacked {
  /// Create an uninitialized VectorPacked.
  TVectorPacked() {}

  /// Create a VectorPacked from a Vector.
  ///
  /// Both VectorPacked and Vector must have the same number of dimensions.
  /// @param vector Vector to create the VectorPacked from.
  explicit TVectorPacked(const TVector<T, d>& vector) { vector.Pack(this); }

  /// Copy a Vector to a VectorPacked.
  ///
  /// Both VectorPacked and Vector must have the same number of dimensions.
  /// @param vector Vector to copy to the VectorPacked.
  /// @returns A reference to this VectorPacked.
  TVectorPacked& operator=(const TVector<T, d>& vector) {
    vector.Pack(this);
    return *this;
  }

  /// Elements of the packed vector one per dimension.
  T data[d];
};
/// @}

/// @addtogroup 向量
/// @{
/// @class TVector "TVector.h"
/// @brief TVector of d elements with type T
///
/// Vector stores <b>d</b> elements of type <b>T</b> and provides a set
/// functions to perform operations on the set of elements.
///
/// @tparam T type of Vector elements.
/// @tparam d dimensions (number of elements) in the Vector structure.
template <class T, int d>
class TVector {
 public:
  /// @brief Element type to enable reference by other classes.
  typedef T Scalar;

  /// @brief Create an uninitialized Vector.
  inline TVector() {
      TMATH_VECTOR_OPERATION(data_[i] = 0);
  }

  /// @brief Create a vector from another vector copying each element.
  ///
  /// @param v Vector that the data will be copied from.
  inline TVector(const TVector<T, d>& v) {
    TMATH_VECTOR_OPERATION(data_[i] = v.data_[i]);
  }

  /// @brief Create a vector from another vector of a different type.
  ///
  /// This copies each element of a Vector which makes it possible to between
  /// vectors of different types, for example
  /// <code>float/double/int</code> vectors.
  /// @param v Vector that the data will be copied from.
  /// @tparam U type of Vector elements to copy.
  template <typename U>
  explicit inline TVector(const TVector<U, d>& v) {
    TMATH_VECTOR_OPERATION(data_[i] = static_cast<T>(v[i]));
  }

  /// @brief Create a vector from a single float.
  ///
  /// Each elements is set to be equal to the value given.
  /// @param s Scalar value that the vector will be initialized to.
  explicit inline TVector(const T& s) { TMATH_VECTOR_OPERATION(data_[i] = s); }

  /// @brief Create a vector form the first d elements of an array.
  ///
  /// @param a Array of values that the vector will be iniitlized to.
  explicit inline TVector(const T* a) {
    TMATH_VECTOR_OPERATION(data_[i] = a[i]);
  }

  /// @brief Create a vector from two values.
  ///
  /// @note This method only works when the vector is of size two.
  ///
  /// @param s1 Scalar value for the first element of the vector.
  /// @param s2 Scalar value for the second element of the vector.
  inline TVector(const T& s1, const T& s2) {
    TMATH_STATIC_ASSERT(d >= 2);
    data_[0] = s1;
    data_[1] = s2;
  }

  /// @brief Create a vector from three values.
  ///
  /// @note This method only works when the vector is of size three.
  ///
  /// @param s1 Scalar value for the first element of the vector.
  /// @param s2 Scalar value for the second element of the vector.
  /// @param s3 Scalar value for the third element of the vector.
  inline TVector(const T& s1, const T& s2, const T& s3) {
    TMATH_STATIC_ASSERT(d >= 3);
    data_[0] = s1;
    data_[1] = s2;
    data_[2] = s3;
  }

  /// @brief Create a vector from a 2 component vector and a third value.
  ///
  /// @note This method only works when the vector is of size three.
  ///
  /// @param v12 Vector containing the first 2 values.
  /// @param s3 Scalar value for the third element of the vector.
  inline TVector(const TVector<T, 2>& v12, const T& s3) {
    TMATH_STATIC_ASSERT(d >= 3);
    data_[0] = v12.x();
    data_[1] = v12.y();
    data_[2] = s3;
  }

  /// @brief Create a vector from four values.
  ///
  /// @note This method only works when the vector is of size four.
  ///
  /// @param s1 Scalar value for the first element of the vector.
  /// @param s2 Scalar value for the second element of the vector.
  /// @param s3 Scalar value for the third element of the vector.
  /// @param s4 Scalar value for the forth element of the vector.
  inline TVector(const T& s1, const T& s2, const T& s3, const T& s4) {
    TMATH_STATIC_ASSERT(d >= 4);
    data_[0] = s1;
    data_[1] = s2;
    data_[2] = s3;
    data_[3] = s4;
  }

  /// @brief Create a 4-dimensional vector from a Vector<T, 3>.
  ///
  /// The last element is initialized to the specified value.
  /// @note This method only works with 4 element vectors.
  ///
  /// @param vector3 Vector used to initialize the first 3 elements.
  /// @param value Value used to set the last element of the vector.
  inline TVector(const TVector<T, 3>& vector3, const T& value) {
    TMATH_STATIC_ASSERT(d == 4);
    data_[0] = vector3[0];
    data_[1] = vector3[1];
    data_[2] = vector3[2];
    data_[3] = value;
  }

  /// @brief Create a vector from two 2 component vectors.
  ///
  /// @note This method only works when the vector is of size four.
  ///
  /// @param v12 Vector containing the first 2 values.
  /// @param v34 Vector containing the last 2 values.
  inline TVector(const TVector<T, 2>& v12, const TVector<T, 2>& v34) {
    TMATH_STATIC_ASSERT(d == 4);
    data_[0] = v12.x();
    data_[1] = v12.y();
    data_[2] = v34.x();
    data_[3] = v34.y();
  }

  /// @brief Create a vector from packed vector (VectorPacked).
  ///
  /// @param vector Packed vector used to initialize an unpacked.
  explicit inline TVector(const TVectorPacked<T, d>& vector) {
    TMATH_VECTOR_OPERATION(data_[i] = vector.data[i]);
  }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data that can be modified by the
  /// caller.
  inline T& operator()(const int i) { return data_[i]; }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data.
  inline const T& operator()(const int i) const { return data_[i]; }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data that can be modified by the
  /// caller.
  inline T& operator[](const int i) { return data_[i]; }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A const reference to the accessed.
  inline const T& operator[](const int i) const { return data_[i]; }

  /// Get the first element (X axis) of the Vector.
  inline T& x() {
    TMATH_STATIC_ASSERT(d > 0);
    return data_[0];
  }
  /// Get the second element (Y axis) of the Vector.
  inline T& y() {
    TMATH_STATIC_ASSERT(d > 1);
    return data_[1];
  }
  /// Get the third element (Z axis) of the Vector.
  inline T& z() {
    TMATH_STATIC_ASSERT(d > 2);
    return data_[2];
  }
  /// Get the fourth element (W axis) of the Vector.
  inline T& w() {
    TMATH_STATIC_ASSERT(d > 3);
    return data_[3];
  }

  /// Get the first element (X axis) of the Vector.
  inline const T& x() const {
    TMATH_STATIC_ASSERT(d > 0);
    return data_[0];
  }
  /// Get the second element (Y axis) of the Vector.
  inline const T& y() const {
    TMATH_STATIC_ASSERT(d > 1);
    return data_[1];
  }
  /// Get the third element (Z axis) of the Vector.
  inline const T& z() const {
    TMATH_STATIC_ASSERT(d > 2);
    return data_[2];
  }
  /// Get the fourth element (W axis) of the Vector.
  inline const T& w() const {
    TMATH_STATIC_ASSERT(d > 3);
    return data_[3];
  }

  /// @brief GLSL style 3 element accessor.
  ///
  /// This only works with vectors that contain more than 3 elements.
  /// @returns A 3-dimensional Vector containing the first 3 elements of
  // this Vector.
  inline TVector<T, 3> xyz() {
    TMATH_STATIC_ASSERT(d > 3);
    return TVector<T, 3>(x(), y(), z());
  }

  /// @brief GLSL style 3 element accessor.
  ///
  /// This only works with vectors that contain more than 3 elements.
  /// @returns A 3-dimensional Vector containing the first 3 elements of
  // this Vector.
  inline const TVector<T, 3> xyz() const {
    TMATH_STATIC_ASSERT(d > 3);
    return TVector<T, 3>(x(), y(), z());
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain more than 2 elements.
  /// @returns A 2-dimensional Vector with the first 2 elements of this Vector.
  inline TVector<T, 2> xy() {
    TMATH_STATIC_ASSERT(d > 2);
    return TVector<T, 2>(x(), y());
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain more than 2 elements.
  /// @returns A 2-dimensional Vector with the first 2 elements of this Vector.
  inline const TVector<T, 2> xy() const {
    TMATH_STATIC_ASSERT(d > 2);
    return TVector<T, 2>(x(), y());
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain 4 elements.
  /// @returns A 2-dimensional Vector with the last 2 elements of this Vector.
  inline TVector<T, 2> zw() {
    TMATH_STATIC_ASSERT(d == 4);
    return TVector<T, 2>(z(), w());
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain 4 elements.
  /// @returns A 2-dimensional Vector with the last 2 elements of this Vector.
  inline const TVector<T, 2> zw() const {
    TMATH_STATIC_ASSERT(d == 4);
    return TVector<T, 2>(z(), w());
  }

  /// @brief Pack a Vector to a packed "d" element vector structure.
  ///
  /// @param vector Packed "d" element vector to write to.
  inline void pack(TVectorPacked<T, d>* const vector) const {
    TMATH_VECTOR_OPERATION(vector->data[i] = data_[i]);
  }

  /// @brief Negate all elements of the Vector.
  ///
  /// @return A new Vector containing the result.
  inline TVector<T, d> operator-() const { TMATH_VECTOR_OPERATOR(-data_[i]); }

  /// @brief Multiply this Vector by another Vector.
  ///
  /// In line with GLSL, this performs component-wise multiplication.
  /// @param v A Vector to mulitply this Vector with.
  /// @return A new Vector containing the result.
  inline TVector<T, d> operator*(const TVector<T, d>& v) const {
    return hadamardProduct(*this, v);
  }

  /// @brief Divide this Vector by another Vector.
  ///
  /// In line with GLSL, this performs component-wise division.
  /// @param v A Vector to divide this Vector by.
  /// @return A new Vector containing the result.
  inline TVector<T, d> operator/(const TVector<T, d>& v) const {
    TMATH_VECTOR_OPERATOR(data_[i] / v[i]);
  }

  /// @brief Add this Vector with another Vector.
  ///
  /// @param v A vector to add this vector with.
  /// @return A new vector containing the result.
  inline TVector<T, d> operator+(const TVector<T, d>& v) const {
    TMATH_VECTOR_OPERATOR(data_[i] + v[i]);
  }

  /// @brief Add this Vector with another Vector.
  ///
  /// @param v A vector to subtract from this vector.
  /// @return A new vector containing the result.
  inline TVector<T, d> operator-(const TVector<T, d>& v) const {
    TMATH_VECTOR_OPERATOR(data_[i] - v[i]);
  }

  /// @brief Multiply this Vector with a scalar.
  ///
  /// @param s A scalar to multiply this vector with.
  /// @return A new vector containing the result.
  inline TVector<T, d> operator*(const T& s) const {
    TMATH_VECTOR_OPERATOR(data_[i] * s);
  }

  /// @brief Divide this Vector by a scalar.
  ///
  /// @param s A scalar to divide this vector with.
  /// @return A new vector containing the result.
  inline TVector<T, d> operator/(const T& s) const {
    TMATH_VECTOR_OPERATOR(data_[i] / s);
  }

  /// @brief Add a scalar to all elements of this Vector.
  ///
  /// @param s A scalar to add to this vector.
  /// @return A new vector containing the result.
  inline TVector<T, d> operator+(const T& s) const {
    TMATH_VECTOR_OPERATOR(data_[i] + s);
  }

  /// @brief Subtract a scalar from all elements of this vector.
  ///
  /// @param s A scalar to subtract from this vector.
  /// @return A new vector that stores the result.
  inline TVector<T, d> operator-(const T& s) const {
    TMATH_VECTOR_OPERATOR(data_[i] - s);
  }

  /// @brief 赋值运算符重载.
  ///
  /// @return A new Vector containing the result.
  inline TVector<T, d>& operator=(const TVector<T, d>& v){ 
      TMATH_VECTOR_SELF_OPERATOR(data_[i] = v[i]);
  }

  /// @brief Multiply (in-place) this Vector with another Vector.
  ///
  /// In line with GLSL, this performs component-wise multiplication.
  /// @param v A vector to multiply this vector with.
  /// @return A reference to this class.
  inline TVector<T, d>& operator*=(const TVector<T, d>& v) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] *= v[i]);
  }

  /// @brief Divide (in-place) this Vector by another Vector.
  ///
  /// In line with GLSL, this performs component-wise division.
  /// @param v A vector to divide this vector by.
  /// @return A reference to this class.
  inline TVector<T, d>& operator/=(const TVector<T, d>& v) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] /= v[i]);
  }

  /// @brief Add (in-place) this Vector with another Vector.
  ///
  /// @param v A vector to add this vector with.
  /// @return A reference to this class.
  inline TVector<T, d>& operator+=(const TVector<T, d>& v) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] += v[i]);
  }

  /// @brief Subtract (in-place) another Vector from this Vector.
  ///
  /// @param v A vector to subtract this vector by.
  /// @return A reference to this class.
  inline TVector<T, d>& operator-=(const TVector<T, d>& v) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] -= v[i]);
  }

  /// @brief Multiply (in-place) each element of this Vector with a scalar.
  ///
  /// @param s A scalar to mulitply this vector with.
  /// @return A reference to this class.
  inline TVector<T, d>& operator*=(const T& s) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] *= s);
  }

  /// @brief Divide (in-place) each element of this Vector by a scalar.
  ///
  /// @param s A scalar to divide this vector by.
  /// @return A reference to this class.
  inline TVector<T, d>& operator/=(const T& s) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] /= s);
  }

  /// @brief Add (in-place) a scalar to each element of this Vector.
  ///
  /// @param s A scalar to add this vector to.
  /// @return A reference to this class.
  inline TVector<T, d>& operator+=(const T& s) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] += s);
  }

  /// @brief Subtract (in-place) a scalar from each element of this Vector.
  ///
  /// @param s A scalar to subtract from this vector.
  /// @return A reference to this class.
  inline TVector<T, d>& operator-=(const T& s) {
    TMATH_VECTOR_SELF_OPERATOR(data_[i] -= s);
  }

  /// @brief Calculate the squared length of this vector.
  ///
  /// @return The length of this vector squared.
  inline T lengthSquared() const { return dotProduct(*this, *this); }

  /// @brief Calculate the length of this vector.
  ///
  /// @return The length of this vector.
  inline T length() const { return sqrt(lengthSquared()); }

  /// @brief Normalize this vector in-place.
  ///
  /// @return The length of this vector.
  inline T normalize() {
    const T len = length();
    *this = *this * (1 / len);
    return len;
  }

  /// @brief Calculate the normalized version of this vector.
  ///
  /// @return The normalized vector.
  inline TVector<T, d> normalized() const { return *this * (1 / length()); }

  /// @brief Calculate the dot product of two vectors.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return The dot product of v1 and v2.
  static inline T dotProduct(const TVector<T, d>& v1, const TVector<T, d>& v2) {
    return DotProductHelper(v1, v2);
  }

  /// @brief Calculate the hadamard or componentwise product of two vectors.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return The hadamard product of v1 and v2.
  static inline TVector<T, d> hadamardProduct(const TVector<T, d>& v1,
                                             const TVector<T, d>& v2) {
    TMATH_VECTOR_OPERATOR(v1[i] * v2[i]);
  }

  /// @brief Calculate the cross product of two vectors.
  ///
  /// Note that this function is only defined for 3-dimensional Vectors.
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return The cross product of v1 and v2.
  static inline TVector<T, 3> crossProduct(const TVector<T, 3>& v1,
                                          const TVector<T, 3>& v2) {
    return TVector<T, 3>(v1[1] * v2[2] - v1[2] * v2[1],
                        v1[2] * v2[0] - v1[0] * v2[2],
                        v1[0] * v2[1] - v1[1] * v2[0]);
  }

  /// @brief Linearly interpolate two vectors.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @param percent Percentage from v1 to v2 in range 0.0...1.0.
  /// @return The hadamard product of v1 and v2.
  static inline TVector<T, d> lerp(const TVector<T, d>& v1,
                                  const TVector<T, d>& v2, const T percent) {
    const T one_minus_percent = static_cast<T>(1.0) - percent;
    TMATH_VECTOR_OPERATOR(one_minus_percent * v1[i] + percent * v2[i]);
  }

  /// @brief Generates a random vector.
  ///
  /// The range of each component is bounded by min and max.
  /// @param min Minimum value of the vector.
  /// @param max Maximum value of the vector.
  static inline TVector<T, d> randomInRange(const TVector<T, d>& min,
                                           const TVector<T, d>& max) {
    TVector<T, d> result;
    TMATH_VECTOR_OPERATION(result[i] =
                                RandomInRange<T>(min[i], max[i]));
    return result;
  }

  /// @brief Compare each component and returns max values.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return Max value of v1 and v2.
  static inline TVector<T, d> max(const TVector<T, d>& v1,
                                 const TVector<T, d>& v2) {
    TVector<T, d> result;
    TMATH_VECTOR_OPERATION(result[i] = std::max(v1[i], v2[i]));
    return result;
  }

  /// @brief Compare each component and returns min values.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return Min value of v1 and v2.
  static inline TVector<T, d> min(const TVector<T, d>& v1,
                                 const TVector<T, d>& v2) {
    TVector<T, d> result;
    TMATH_VECTOR_OPERATION(result[i] = std::min(v1[i], v2[i]));
    return result;
  }

 private:
  /// Elements of the vector.
  T data_[d];
};
/// @}

/// @addtogroup TMath_vector
/// @{

/// @brief Multiply a Vector by a scalar.
///
/// Multiplies each component of the specified Vector with a scalar.
///
/// @param s scalar to multiply.
/// @param v Vector to multiply.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline TVector<T, d> operator*(const T& s, const TVector<T, d>& v) {
  return v * s;
}

/// @brief Divide a Vector by a scalar.
///
/// Divides each component of the specified Vector by a scalar.
///
/// @param v Vector to be divided.
/// @param s scalar to divide the vector by.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline TVector<T, d> operator/(const TVector<T, d>& v, const T& s) {
  return v / s;
}

/// @brief Add a scalar to each element of a Vector.
///
/// @param s scalar to add to each element of a Vector.
/// @param v Vector to add the scalar to.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline TVector<T, d> operator+(const T& s, const TVector<T, d>& v) {
  return v + s;
}

/// @brief Subtract a scalar from each element of a Vector.
///
/// @param s scalar to subtract from each element of a Vector.
/// @param v Vector to subtract the scalar from.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline TVector<T, d> operator-(const T& s, const TVector<T, d>& v) {
  return v - s;
}

/// @brief Check if val is within [range_start..range_end), denoting a
/// rectangular area.
///
/// @param val 2D vector to be tested.
/// @param range_start Starting point of the range (inclusive).
/// @param range_end Ending point of the range (non-inclusive).
/// @return Bool indicating success.
///
/// @tparam T Type of vector components to test.
template <class T>
bool InRange2D(const TVector<T, 2>& val,
               const TVector<T, 2>& range_start,
               const TVector<T, 2>& range_end) {
  return InRange(val.x(), range_start.x(), range_end.x()) &&
         InRange(val.y(), range_start.y(), range_end.y());
}

/// @cond 内部
/// @brief Calculate the dot product of two vectors.
///
/// @param v1 First vector.
/// @param v2 Second vector.
/// @return The dot product of v1 and v2.
/// @related Vector
template <class T, int d>
static inline T DotProductHelper(const TVector<T, d>& v1,
                                 const TVector<T, d>& v2) {
  T result = 0;
  TMATH_VECTOR_OPERATION(result += v1[i] * v2[i]);
  return result;
}
/// @endcond

/// @cond 内部
template <class T>
static inline T DotProductHelper(const TVector<T, 2>& v1,
                                 const TVector<T, 2>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1];
}
/// @endcond

/// @cond 内部
template <class T>
static inline T DotProductHelper(const TVector<T, 3>& v1,
                                 const TVector<T, 3>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
/// @endcond

/// @cond 内部
template <class T>
static inline T DotProductHelper(const TVector<T, 4>& v1,
                                 const TVector<T, 4>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}
/// @endcond

/// @}

TMATH_END_NAMESPACE

#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

#endif // TVECTOR_H

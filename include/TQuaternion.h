#ifndef TQUATERNION_H
#define TQUATERNION_H


#ifdef _WIN32
#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES  // For M_PI.
#endif                     // !defined(_USE_MATH_DEFINES)
#endif                     // _WIN32

#include <cmath>

#include "TMatrix.h"


/// @file TQuaternion.h
/// @brief Quaternion class and functions.
/// @addtogroup TQuaternion
///
/// TMath provides a Quaternion class that utilizes SIMD optimized
/// Matrix and Vector classes.

TMATH_BEGIN_NAMESPACE

/// @addtogroup TQuaternion
/// @{
/// @class Quaternion
///
/// @brief Stores a Quaternion of type T and provides a set of utility
/// operations on each Quaternion.
/// @tparam T Type of each element in the Quaternion.
template <class T>
class TQuaternion {
 public:
  /// @brief Construct an uninitialized Quaternion.
  inline TQuaternion() {}

  /// @brief Construct a Quaternion from a copy.
  /// @param q Quaternion to copy.
  inline TQuaternion(const TQuaternion<T>& q) {
    s_ = q.s_;
    v_ = q.v_;
  }

  /// @brief Construct a Quaternion using scalar values to initialize each
  /// element.
  ///
  /// @param s1 Scalar component.
  /// @param s2 First element of the Vector component.
  /// @param s3 Second element of the Vector component.
  /// @param s4 Third element of the Vector component.
  inline TQuaternion(const T& s1, const T& s2, const T& s3, const T& s4) {
    s_ = s1;
    v_ = TVector<T, 3>(s2, s3, s4);
  }

  /// @brief Construct a quaternion from a scalar and 3-dimensional Vector.
  ///
  /// @param s1 Scalar component.
  /// @param v1 Vector component.
  inline TQuaternion(const T& s1, const TVector<T, 3>& v1) {
    s_ = s1;
    v_ = v1;
  }

  /// @brief Return the scalar component of the quaternion.
  ///
  /// @return The scalar component
  inline T& scalar() { return s_; }

  /// @brief Return the scalar component of the quaternion.
  ///
  /// @return The scalar component
  inline const T& scalar() const { return s_; }

  /// @brief Set the scalar component of the quaternion.
  ///
  /// @param s Scalar component.
  inline void set_scalar(const T& s) { s_ = s; }

  /// @brief Return the vector component of the quaternion.
  ///
  /// @return The scalar component
  inline TVector<T, 3>& vector() { return v_; }

  /// @brief Return the vector component of the quaternion.
  ///
  /// @return The scalar component
  inline const TVector<T, 3>& vector() const { return v_; }

  /// @brief Set the vector component of the quaternion.
  ///
  /// @param v Vector component.
  inline void set_vector(const TVector<T, 3>& v) { v_ = v; }

  /// @brief Calculate the inverse Quaternion.
  ///
  /// This calculates the inverse such that <code>(q * q).Inverse()</code>
  /// is the identity.
  ///
  /// @return Quaternion containing the result.
  inline TQuaternion<T> Inverse() const { return TQuaternion<T>(s_, -v_); }

  /// @brief Multiply this Quaternion with another Quaternion.
  ///
  /// @note This is equivalent to
  /// <code>FromMatrix(ToMatrix() * q.ToMatrix()).</code>
  /// @param q Quaternion to multiply with.
  /// @return Quaternion containing the result.
  inline TQuaternion<T> operator*(const TQuaternion<T>& q) const {
    return TQuaternion<T>(
        s_ * q.s_ - Vector<T, 3>::DotProduct(v_, q.v_),
        s_ * q.v_ + q.s_ * v_ + Vector<T, 3>::CrossProduct(v_, q.v_));
  }

  /// @brief Multiply this Quaternion by a scalar.
  ///
  /// This multiplies the angle of the rotation by a scalar factor.
  /// @param s1 Scalar to multiply with.
  /// @return Quaternion containing the result.
  inline TQuaternion<T> operator*(const T& s1) const {
    T angle;
    TVector<T, 3> axis;
    ToAngleAxis(&angle, &axis);
    angle *= s1;
    return Quaternion<T>(cos(0.5f * angle),
                         axis.Normalized() * sin(0.5f * angle));
  }

  /// @brief Multiply a Vector by this Quaternion.
  ///
  /// This will rotate the specified vector by the rotation specified by this
  /// Quaternion.
  ///
  /// @param v1 Vector to multiply by this Quaternion.
  /// @return Rotated Vector.
  inline TVector<T, 3> operator*(const TVector<T, 3>& v1) const {
    T ss = s_ + s_;
    return ss * TVector<T, 3>::CrossProduct(v_, v1) + (ss * s_ - 1) * v1 +
           2 * TVector<T, 3>::DotProduct(v_, v1) * v_;
  }

  /// @brief Normalize this quaterion (in-place).
  ///
  /// @return Length of the quaternion.
  inline T Normalize() {
    T length = sqrt(s_ * s_ + TVector<T, 3>::DotProduct(v_, v_));
    T scale = (1 / length);
    s_ *= scale;
    v_ *= scale;
    return length;
  }

  /// @brief Calculate the normalized version of this quaternion.
  ///
  /// @return The normalized quaternion.
  inline TQuaternion<T> Normalized() {
    TQuaternion<T> q(*this);
    q.Normalize();
    return q;
  }

  /// @brief Convert this Quaternion to an Angle and axis.
  ///
  /// The returned  angle is the size of the rotation in radians about the
  /// axis represented by this Quaternion.
  ///
  /// @param angle Receives the angle.
  /// @param axis Receives the normalized axis.
  inline void ToAngleAxis(T* angle, TVector<T, 3>* axis) const {
    *axis = s_ > 0 ? v_ : -v_;
    *angle = 2 * atan2(axis->Normalize(), s_ > 0 ? s_ : -s_);
  }

  /// @brief Convert this Quaternion to 3 Euler Angles.
  ///
  /// @return 3-dimensional Vector where each element is a angle of rotation
  /// (in radians) around the x, y, and z axes.
  inline TVector<T, 3> ToEulerAngles() const {
    TMatrix<T, 3> m(ToMatrix());
    T cos2 = m[0] * m[0] + m[1] * m[1];
    if (cos2 < 1e-6f) {
      return TVector<T, 3>(0, m[2] < 0 ? static_cast<T>(0.5 * M_PI)
                                      : static_cast<T>(-0.5 * M_PI),
                          -atan2(m[3], m[4]));
    } else {
      return TVector<T, 3>(atan2(m[5], m[8]), atan2(-m[2], sqrt(cos2)),
                          atan2(m[1], m[0]));
    }
  }

  /// @brief Convert to a 3x3 Matrix.
  ///
  /// @return 3x3 rotation Matrix.
  inline TMatrix<T, 3> ToMatrix() const {
    const T x2 = v_[0] * v_[0], y2 = v_[1] * v_[1], z2 = v_[2] * v_[2];
    const T sx = s_ * v_[0], sy = s_ * v_[1], sz = s_ * v_[2];
    const T xz = v_[0] * v_[2], yz = v_[1] * v_[2], xy = v_[0] * v_[1];
    return TMatrix<T, 3>(1 - 2 * (y2 + z2), 2 * (xy + sz), 2 * (xz - sy),
                        2 * (xy - sz), 1 - 2 * (x2 + z2), 2 * (sx + yz),
                        2 * (sy + xz), 2 * (yz - sx), 1 - 2 * (x2 + y2));
  }

  /// @brief Convert to a 4x4 Matrix.
  ///
  /// @return 4x4 transform Matrix.
  inline TMatrix<T, 4> ToMatrix4() const {
    const T x2 = v_[0] * v_[0], y2 = v_[1] * v_[1], z2 = v_[2] * v_[2];
    const T sx = s_ * v_[0], sy = s_ * v_[1], sz = s_ * v_[2];
    const T xz = v_[0] * v_[2], yz = v_[1] * v_[2], xy = v_[0] * v_[1];
    return TMatrix<T, 4>(1 - 2 * (y2 + z2), 2 * (xy + sz), 2 * (xz - sy), 0.0f,
                        2 * (xy - sz), 1 - 2 * (x2 + z2), 2 * (sx + yz), 0.0f,
                        2 * (sy + xz), 2 * (yz - sx), 1 - 2 * (x2 + y2), 0.0f,
                        0.0f, 0.0f, 0.0f, 1.0f);
  }

  /// @brief Create a Quaternion from an angle and axis.
  ///
  /// @param angle Angle in radians to rotate by.
  /// @param axis Axis in 3D space to rotate around.
  /// @return Quaternion containing the result.
  static TQuaternion<T> FromAngleAxis(const T& angle, const TVector<T, 3>& axis) {
    const T halfAngle = static_cast<T>(0.5) * angle;
    TVector<T, 3> localAxis(axis);
    return Quaternion<T>(cos(halfAngle),
                         localAxis.Normalized() * sin(halfAngle));
  }

  /// @brief Create a quaternion from 3 euler angles.
  ///
  /// @param angles 3-dimensional Vector where each element contains an
  /// angle in radius to rotate by about the x, y and z axes.
  /// @return Quaternion containing the result.
  static TQuaternion<T> FromEulerAngles(const TVector<T, 3>& angles) {
    const TVector<T, 3> halfAngles(static_cast<T>(0.5) * angles[0],
                                  static_cast<T>(0.5) * angles[1],
                                  static_cast<T>(0.5) * angles[2]);
    const T sinx = sin(halfAngles[0]), cosx = cos(halfAngles[0]);
    const T siny = sin(halfAngles[1]), cosy = cos(halfAngles[1]);
    const T sinz = sin(halfAngles[2]), cosz = cos(halfAngles[2]);
    return TQuaternion<T>(cosx * cosy * cosz + sinx * siny * sinz,
                         sinx * cosy * cosz - cosx * siny * sinz,
                         cosx * siny * cosz + sinx * cosy * sinz,
                         cosx * cosy * sinz - sinx * siny * cosz);
  }

  /// @brief Create a quaternion from a rotation Matrix.
  ///
  /// @param m 3x3 rotation Matrix.
  /// @return Quaternion containing the result.
  static TQuaternion<T> FromMatrix(const TMatrix<T, 3>& m) {
    const T trace = m(0, 0) + m(1, 1) + m(2, 2);
    if (trace > 0) {
      const T s = sqrt(trace + 1) * 2;
      const T oneOverS = 1 / s;
      return TQuaternion<T>(static_cast<T>(0.25) * s, (m[5] - m[7]) * oneOverS,
                           (m[6] - m[2]) * oneOverS, (m[1] - m[3]) * oneOverS);
    } else if (m[0] > m[4] && m[0] > m[8]) {
      const T s = sqrt(m[0] - m[4] - m[8] + 1) * 2;
      const T oneOverS = 1 / s;
      return TQuaternion<T>((m[5] - m[7]) * oneOverS, static_cast<T>(0.25) * s,
                           (m[3] + m[1]) * oneOverS, (m[6] + m[2]) * oneOverS);
    } else if (m[4] > m[8]) {
      const T s = sqrt(m[4] - m[0] - m[8] + 1) * 2;
      const T oneOverS = 1 / s;
      return TQuaternion<T>((m[6] - m[2]) * oneOverS, (m[3] + m[1]) * oneOverS,
                           static_cast<T>(0.25) * s, (m[5] + m[7]) * oneOverS);
    } else {
      const T s = sqrt(m[8] - m[0] - m[4] + 1) * 2;
      const T oneOverS = 1 / s;
      return TQuaternion<T>((m[1] - m[3]) * oneOverS, (m[6] + m[2]) * oneOverS,
                           (m[5] + m[7]) * oneOverS, static_cast<T>(0.25) * s);
    }
  }

  /// @brief Calculate the spherical linear interpolation between two
  /// Quaternions.
  ///
  /// @param q1 Start Quaternion.
  /// @param q2 End Quaternion.
  /// @param s1 The scalar value determining how far from q1 and q2 the
  /// resulting quaternion should be.  A value of 0 corresponds to q1 and a
  /// value of 1 corresponds to q2.
  /// @result Quaternion containing the result.
  static inline TQuaternion<T> Slerp(const TQuaternion<T>& q1,
                                    const TQuaternion<T>& q2, const T& s1) {
    if (q1.s_ * q2.s_ + Vector<T, 3>::DotProduct(q1.v_, q2.v_) > 0.999999f)
      return TQuaternion<T>(q1.s_ * (1 - s1) + q2.s_ * s1,
                           q1.v_ * (1 - s1) + q2.v_ * s1);
    return q1 * ((q1.Inverse() * q2) * s1);
  }

  /// @brief Access an element of the quaternion.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data that can be modified by the
  /// caller.
  inline T& operator[](const int i) {
    if (i == 0) return s_;
    return v_[i - 1];
  }

  /// @brief Access an element of the quaternion.
  ///
  /// @param i Index of the element to access.
  /// @return A const reference to the accessed.
  inline const T& operator[](const int i) const {
    return const_cast<TQuaternion<T>*>(this)->operator[](i);
  }

  /// @brief Returns a vector that is perpendicular to the supplied vector.
  ///
  /// @param v1 An arbitrary vector
  /// @return A vector perpendicular to v1.  Normally this will just be
  /// the cross product of v1, v2.  If they are parallel or opposite though,
  /// the routine will attempt to pick a vector.
  static inline TVector<T, 3> PerpendicularVector(const TVector<T, 3>& v) {
    // We start out by taking the cross product of the vector and the x-axis to
    // find something parallel to the input vectors.  If that cross product
    // turns out to be length 0 (i. e. the vectors already lie along the x axis)
    // then we use the y-axis instead.
    TVector<T, 3> axis = TVector<T, 3>::CrossProduct(
        TVector<T, 3>(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0)),
        v);
    // We use a fairly high epsilon here because we know that if this number
    // is too small, the axis we'll get from a cross product with the y axis
    // will be much better and more numerically stable.
    if (axis.LengthSquared() < static_cast<T>(0.05)) {
      axis = TVector<T, 3>::CrossProduct(
          TVector<T, 3>(static_cast<T>(0), static_cast<T>(1), static_cast<T>(0)),
          v);
    }
    return axis;
  }

  /// @brief Returns the a Quaternion that rotates from start to end.
  ///
  /// @param v1 The starting vector
  /// @param v2 The vector to rotate to
  /// @param preferred_axis the axis to use, if v1 and v2 are parallel.
  /// @return A Quaternion describing the rotation from v1 to v2
  /// See the comment on RotateFromToWithAxis for an explanation of the math.
  static inline TQuaternion<T> RotateFromToWithAxis(
      const TVector<T, 3>& v1, const TVector<T, 3>& v2,
      const TVector<T, 3>& preferred_axis) {
    TVector<T, 3> start = v1.Normalized();
    TVector<T, 3> end = v2.Normalized();

    T dot_product = TVector<T, 3>::DotProduct(start, end);
    // Any rotation < 0.1 degrees is treated as no rotation
    // in order to avoid division by zero errors.
    // So we early-out in cases where it's less then 0.1 degrees.
    // cos( 0.1 degrees) = 0.99999847691
    if (dot_product >= static_cast<T>(0.99999847691)) {
      return TQuaternion<T>::identity;
    }
    // If the vectors point in opposite directions, return a 180 degree
    // rotation, on the axis that they asked for.
    if (dot_product <= static_cast<T>(-0.99999847691)) {
      return TQuaternion<T>(static_cast<T>(0), preferred_axis);
    }
    // Degenerate cases have been handled, so if we're here, we have to
    // actually compute the angle we want:
    TVector<T, 3> cross_product = TVector<T, 3>::CrossProduct(start, end);

    return TQuaternion<T>(static_cast<T>(1.0) + dot_product, cross_product)
        .Normalized();
  }

  /// @brief Returns the a Quaternion that rotates from start to end.
  ///
  /// @param v1 The starting vector
  /// @param v2 The vector to rotate to
  /// @return A Quaternion describing the rotation from v1 to v2.  In the case
  /// where the vectors are parallel, it returns the identity.  In the case
  /// where
  /// they point in opposite directions, it picks an arbitrary axis.  (Since
  /// there
  /// are technically infinite possible quaternions to represent a 180 degree
  /// rotation.)
  ///
  /// The final equation used here is fairly elegant, but its derivation is
  /// not obvious:  We want to find the quaternion that represents the angle
  /// between Start and End.
  ///
  /// The angle can be expressed as a quaternion with the values:
  ///  angle: ArcCos(dotproduct(start, end) / (|start|*|end|)
  ///  axis: crossproduct(start, end).normalized * sin(angle/2)
  ///
  /// or written as:
  ///  quaternion(cos(angle/2), axis * sin(angle/2))
  ///
  /// Using the trig identity:
  ///  sin(angle * 2) = 2 * sin(angle) * cos*angle)
  /// Via substitution, we can turn this into:
  ///  sin(angle/2) = 0.5 * sin(angle)/cos(angle/2)
  ///
  /// Using this substitution, we get:
  ///  quaternion( cos(angle/2),
  ///             0.5 * crossproduct(start, end).normalized
  ///                 * sin(angle) / cos(angle/2))
  ///
  /// If we scale the whole thing up by 2 * cos(angle/2) then we get:
  ///  quaternion(2 * cos(angle/2) * cos(angle/2),
  ///             crossproduct(start, end).normalized * sin(angle))
  ///
  /// (Note that the quaternion is no longer normalized after this scaling)
  ///
  /// Another trig identity:
  ///   cos(angle/2) = sqrt((1 + cos(angle) / 2)
  ///
  /// Substituting this in, we can simplify the quaternion scalar:
  ///
  ///  quaternion(1 + cos(angle),
  ///             crossproduct(start, end).normalized * sin(angle))
  ///
  /// Because cross(start, end) has a magnitude of |start|*|end|*sin(angle),
  ///  crossproduct(start,end).normalized
  /// is equivalent to
  ///  crossproduct(start,end) / |start| * |end| * sin(angle)
  /// So after that substitution:
  ///
  ///  quaternion(1 + cos(angle),
  ///             crossproduct(start, end) / (|start| * |end|))
  ///
  /// dotproduct(start, end) has the value of |start| * |end| * cos(angle),
  /// so by algebra,
  ///  cos(angle) = dotproduct(start, end) / (|start| * |end|)
  /// we can replace our quaternion scalar here also:
  ///
  ///  quaternion(1 + dotproduct(start, end) / (|start| * |end|),
  ///             crossproduct(start, end) / (|start| * |end|))
  ///
  /// If start and end are normalized, then |start| * |end| = 1, giving us a
  /// final quaternion of:
  ///
  /// quaternion(1 + dotproduct(start, end), crossproduct(start, end))
  static inline TQuaternion<T> RotateFromTo(const TVector<T, 3>& v1,
                                           const TVector<T, 3>& v2) {
    TVector<T, 3> start = v1.Normalized();
    TVector<T, 3> end = v2.Normalized();

    T dot_product = TVector<T, 3>::DotProduct(start, end);
    // Any rotation < 0.1 degrees is treated as no rotation
    // in order to avoid division by zero errors.
    // So we early-out in cases where it's less then 0.1 degrees.
    // cos( 0.1 degrees) = 0.99999847691
    if (dot_product >= static_cast<T>(0.99999847691)) {
      return TQuaternion<T>::identity;
    }
    // If the vectors point in opposite directions, return a 180 degree
    // rotation, on an arbitrary axis.
    if (dot_product <= static_cast<T>(-0.99999847691)) {
      return TQuaternion<T>(0, PerpendicularVector(start));
    }
    // Degenerate cases have been handled, so if we're here, we have to
    // actually compute the angle we want:
    TVector<T, 3> cross_product = TVector<T, 3>::CrossProduct(start, end);

    return TQuaternion<T>(static_cast<T>(1.0) + dot_product, cross_product)
        .Normalized();
  }

  /// @brief Contains a quaternion doing the identity transform.
  static TQuaternion<T> identity;

 private:
  T s_;
  TVector<T, 3> v_;
};

template <typename T>
TQuaternion<T> TQuaternion<T>::identity = TQuaternion<T>(1, 0, 0, 0);
/// @}

/// @addtogroup TQuaternion
/// @{

/// @brief Multiply a Quaternion by a scalar.
///
/// This multiplies the angle of the rotation of the specified Quaternion
/// by a scalar factor.
/// @param s Scalar to multiply with.
/// @param q Quaternion to scale.
/// @return Quaternion containing the result.
///
/// @related Quaternion
template <class T>
inline TQuaternion<T> operator*(const T& s, const TQuaternion<T>& q) {
  return q * s;
}
/// @}

TMATH_END_NAMESPACE

#endif // TQUATERNION_H

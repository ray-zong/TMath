#ifndef TUTILITIES_H
#define TUTILITIES_H


#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <string>

#include "TGlobal.h"

/// @file TUtilities.h Utilities
/// @brief Utility macros and functions.
///
///

/// @cond 内部
template <bool>
struct static_assert_util;
template <>
struct static_assert_util<true> {};
/// @endcond

/// @addtogroup TMath_utilities
/// @{
/// @def TMATH_STATIC_ASSERT
/// @brief Compile time assert for pre-C++11 compilers.
///
/// For example:
/// <blockquote><code>
/// TMATH_STATIC_ASSERT(0 == 1);
/// </code></blockquote> will result in a compile error.
#define TMATH_STATIC_ASSERT(x) static_assert_util<(x)>()
/// @}

/// @cond 内部
/// Unroll an loop up to 4 iterations, where iterator is the identifier
/// used in each operation (e.g "i"), number_of_iterations is a constant which
/// specifies the number of times to perform the operation and "operation" is
/// the statement to execute for each iteration of the loop (e.g data[i] = v).
#define TMATH_UNROLLED_LOOP(iterator, number_of_iterations, operation) \
  {                                                                     \
    const int iterator = 0;                                             \
    { operation; }                                                      \
    if ((number_of_iterations) > 1) {                                   \
      const int iterator = 1;                                           \
      { operation; }                                                    \
      if ((number_of_iterations) > 2) {                                 \
        const int iterator = 2;                                         \
        { operation; }                                                  \
        if ((number_of_iterations) > 3) {                               \
          const int iterator = 3;                                       \
          { operation; }                                                \
          if ((number_of_iterations) > 4) {                             \
            for (int iterator = 4; iterator < (number_of_iterations);   \
                 ++iterator) {                                          \
              operation;                                                \
            }                                                           \
          }                                                             \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }
/// @endcond

TMATH_BEGIN_NAMESPACE

/// @addtogroup TMath_utilities
/// @{

/// @brief Clamp x within [lower, upper].
/// @anchor TMath_Clamp
///
/// @note Results are undefined if lower > upper.
///
/// @param x Value to clamp.
/// @param lower Lower value of the range.
/// @param upper Upper value of the range.
/// @returns Clamped value.
template <class T>
T Clamp(const T &x, const T &lower, const T &upper) {
  return std::max<T>(lower, std::min<T>(x, upper));
}

/// @brief Linearly interpolate between range_start and range_end, based on
/// percent.
/// @anchor TMath_Lerp
///
/// @param range_start Start of the range.
/// @param range_end End of the range.
/// @param percent Value between 0.0 and 1.0 used to interpolate between
/// range_start and range_end.  Where a value of 0.0 results in a return
/// value of range_start and 1.0 results in a return value of range_end.
/// @return Value between range_start and range_end.
///
/// @tparam T Type of the range to interpolate over.
/// @tparam T2 Type of the value used to perform interpolation
///         (e.g float or double).
template <class T, class T2>
T Lerp(const T &range_start, const T &range_end, const T2 &percent) {
  const T2 one_minus_percent = static_cast<T2>(1.0) - percent;
  return range_start * one_minus_percent + range_end * percent;
}

/// @brief Linearly interpolate between range_start and range_end, based on
/// percent.
/// @anchor TMath_Lerp2
///
/// @param range_start Start of the range.
/// @param range_end End of the range.
/// @param percent Value between 0.0 and 1.0 used to interpolate between
/// range_start and range_end.  Where a value of 0.0 results in a return
/// value of range_start and 1.0 results in a return value of range_end.
/// @return Value between range_start and range_end.
///
/// @tparam T Type of the range to interpolate over.
template <class T>
T Lerp(const T &range_start, const T &range_end, const T &percent) {
  return Lerp<T, T>(range_start, range_end, percent);
}

/// @brief Check if val is within [range_start..range_end).
/// @anchor TMath_InRange
///
/// @param val Value to be tested.
/// @param range_start Starting point of the range (inclusive).
/// @param range_end Ending point of the range (non-inclusive).
/// @return Bool indicating success.
///
/// @tparam T Type of values to test.
template <class T>
bool InRange(T val, T range_start, T range_end) {
  return val >= range_start && val < range_end;
}

/// @brief  Generate a random value of type T.
/// @anchor TMath_Random
///
/// This method generates a random value of type T, greater than or equal to
/// 0.0 and less than 1.0.
///
/// This function uses the standard C library function rand() from math.h to
/// generate the random number.
///
/// @returns Random number greater than or equal to 0.0 and less than 1.0.
///
/// @see RandomRange()
/// @see RandomInRange()
template <class T>
inline T Random() {
  return static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
}

/// @cond 内部
template <>
inline float Random() {
  return static_cast<float>(rand() >> 8) /
         (static_cast<float>((RAND_MAX >> 8) + 1));
}
/// @endcond

/// @cond 内部
template <>
inline double Random() {
  return static_cast<double>(rand()) / (static_cast<double>(RAND_MAX + 1LL));
}
/// @endcond

/// @brief Generate a random value of type T in the range -range...+range
/// @anchor TMath_RandomRange
///
/// This function uses the standard C library function rand() from math.h to
/// generate the random number.
///
/// @param range Range of the random value to generate.
/// @return Random value in the range -range to +range
///
/// @see Random()
template <class T>
inline T RandomRange(T range) {
  return (Random<T>() * range * 2) - range;
}

/// @brief Generate a random number between [range_start, range_end]
/// @anchor TMath_RandomInRange
///
/// This function uses the standard C library function rand() from math.h to
/// generate the random number.
///
/// @param range_start Minimum value.
/// @param range_end Maximum value.
/// @return Random value in the range [range_start, range_end].
///
/// @see Lerp()
/// @see Random()
template <class T>
inline T RandomInRange(T range_start, T range_end) {
  return Lerp(range_start, range_end, Random<T>());
}

/// @cond 内部
template <>
inline int RandomInRange<int>(int range_start, int range_end) {
  return static_cast<int>(RandomInRange<float>(static_cast<float>(range_start),
                                               static_cast<float>(range_end)));
}
/// @endcond

/// @brief Round a value up to the nearest power of 2.
///
/// @param x Value to round up.
/// @returns Value rounded up to the nearest power of 2.
template <class T>
T RoundUpToPowerOf2(T x) {
  return static_cast<T>(
      pow(static_cast<T>(2), ceil(log(x) / log(static_cast<T>(2)))));
}

/// @brief 字符串转为数字.
///
/// @param str 字符串.
/// @param val 数字.
/// @returns 状态量,转换是否成功.
template <typename T>
bool FromString(const std::string& str, T& val)
{
    std::stringstream ss(str);
    return !(ss >> val).fail();
}

/// @brief 求余运算(恒为正)
/// Mod(-3,4)= 1   fmod(-3,4)= -3
///
/// @param x 
/// @param y 
/// @return 运算结果(恒为正)
inline TMATH_EXPORT double Mod(const double x, const double y);

inline TMATH_EXPORT double WrapNegPosPI(const double a);

inline TMATH_EXPORT double WrapTwoPI(const double a);

inline TMATH_EXPORT double WrapNegPos180(const double a);

inline TMATH_EXPORT double Wrap360(const double a);

inline TMATH_EXPORT double DegreesToRadians(const double degrees);

inline TMATH_EXPORT double RadiansToDegrees(const double radians);

inline TMATH_EXPORT double AcTan(const double sinx, const double cosx);

TMATH_EXPORT void TrimLeft(std::string &s);

TMATH_EXPORT void TrimRight(std::string &s);

TMATH_EXPORT void Trim(std::string &s);

/// @}

#if defined(_MSC_VER)
#if _MSC_VER <= 1800  // MSVC 2013
#if !defined(noexcept)
#define noexcept
#endif  // !defined(noexcept)
#endif  // _MSC_VER <= 1800
#endif  //  defined(_MSC_VER)

/// @}

TMATH_END_NAMESPACE


#endif // TUTILITIES_H

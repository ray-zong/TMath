#ifndef TMATH_GLOBAL_H
#define TMATH_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(TMATH_LIBRARY)
#  define TMATHSHARED_EXPORT Q_DECL_EXPORT
#else
#  define TMATHSHARED_EXPORT Q_DECL_IMPORT
#endif

#define TMATH_VERSION_MAJOR 1
#define TMATH_VERSION_MINOR 0
#define TMATH_VERSION_REVISION 0

#define RAY_BEGIN_NAMESPACE namespace Ray {
#define RAY_END_NAMESPACE }


template <class T>
inline T randomInRange(T range_start,
                       T range_end)
{
  return Lerp(range_start, range_end, Random<T>());
}

template <>
inline int randomInRange<int>(int range_start,
                              int range_end)
{
  return static_cast<int>(randomInRange<float>(static_cast<float>(range_start),
                                               static_cast<float>(range_end)));
}

#endif // TMATH_GLOBAL_H

#ifndef TMATH_GLOBAL_H
#define TMATH_GLOBAL_H

#if defined(TMATH_LIBRARY)
  #define TMATH_EXPORT __declspec(dllexport)
#else
  #define TMATH_EXPORT __declspec(dllimport)
#endif

#define TMATH_BEGIN_NAMESPACE namespace TMath {
#define TMATH_END_NAMESPACE }

/// @addtogroup TMath_version
/// @{

/// @def TMATH_VERSION_MAJOR
/// @brief Major version number of the library.
/// @see kMathVersionString
#define TMATH_VERSION_MAJOR 1
/// @def MATH_VERSION_MINOR
/// @brief Minor version number of the library.
/// @see kMathVersionString
#define TMATH_VERSION_MINOR 0
/// @def TMATH_VERSION_REVISION
/// @brief Revision number of the library.
/// @see kMathVersionString
#define TMATH_VERSION_REVISION 0

/// @}

/// @cond 内部
#define TMATH_STRING_EXPAND(X) #X
#define TMATH_STRING(X) TMATH_STRING_EXPAND(X)
/// @endcond


// Weak linkage is culled by VS & doesn't work on cygwin.
#if !defined(_WIN32) && !defined(__CYGWIN__)

extern volatile __attribute__((weak)) const char *TMathVersionString;
/// @addtogroup TMath_version
/// @{

/// @var kTMathVersionString
/// @brief String which identifies the current version of TMath.
///
/// @ref kTMathVersionString is used by Google developers to identify which
/// applications uploaded to Google Play are using this library.  This allows
/// the development team at Google to determine the popularity of the library.
/// How it works: Applications that are uploaded to the Google Play Store are
/// scanned for this version string.  We track which applications are using it
/// to measure popularity.  You are free to remove it (of course) but we would
/// appreciate if you left it in.
///
/// @see TMATH_VERSION_MAJOR
/// @see TMATH_VERSION_MINOR
/// @see TMATH_VERSION_REVISION
volatile __attribute__((weak)) const char *TMathVersionString =
    "TMath " TMATH_STRING(TMATH_VERSION_MAJOR) "." TMATH_STRING(
        TMATH_VERSION_MINOR) "." TMATH_STRING(TMATH_VERSION_REVISION);
/// @}

#endif  // !defined(_WIN32) && !defined(__CYGWIN__)

#endif // TMATH_GLOBAL_H

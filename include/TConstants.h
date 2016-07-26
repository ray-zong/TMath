#ifndef TCONSTANTS_H
#define TCONSTANTS_H

#include "TQuaternion.h"
#include "TVector.h"

TMATH_BEGIN_NAMESPACE

/// @file TConstants.h
/// @brief 特定维度的向量常量.
/// @addtogroup TConstants
///
/// 向量常量比自定义（调用构造函数）效率更高。在加载到内存时，构造函数大多数情况下会更慢。
/// <p>
/// 例如, 这个:<br>
/// <code>
/// lookat = mat4::LookAt(target, position, ray::kAxisY3f);
/// </code>
/// <br>比:<br>
/// <code>
/// lookat = mat4::LookAt(target, position,
///                       ray::Vector<float, 3>(0.0f, 1.0f, 0.0f));
/// </code>
/// <br> 更加高效简洁.
/// </p>
///
/// Depending on your linker's sophistication and settings, these constants may
/// be duplicated in every compilation unit in which they're used. However,
/// most linkers should be able to detect and eliminate this duplication.

/// @addtogroup TConstants
/// @{

/// 2-dimensional <code>float</code> Vector of zeros.
static const TVector<float, 2> kZeros2f(0.0f, 0.0f);
/// 2-dimensional <code>float</code> Vector of ones.
static const TVector<float, 2> kOnes2f(1.0f, 1.0f);
/// 2-dimensional <code>float</code> unit Vector pointing along the X axis.
static const TVector<float, 2> kAxisX2f(1.0f, 0.0f);
/// 2-dimensional <code>float</code> unit Vector pointing along the Y axis.
static const TVector<float, 2> kAxisY2f(0.0f, 1.0f);

/// 3-dimensional <code>float</code> Vector of zeros.
static const TVector<float, 3> kZeros3f(0.0f, 0.0f, 0.0f);
/// 3-dimensional <code>float</code> Vector of ones.
static const TVector<float, 3> kOnes3f(1.0f, 1.0f, 1.0f);
/// 3-dimensional <code>float</code> unit Vector pointing along the X axis.
static const TVector<float, 3> kAxisX3f(1.0f, 0.0f, 0.0f);
/// 3-dimensional <code>float</code> unit Vector pointing along the Y axis.
static const TVector<float, 3> kAxisY3f(0.0f, 1.0f, 0.0f);
/// 3-dimensional <code>float</code> unit Vector pointing along the Z axis.
static const TVector<float, 3> kAxisZ3f(0.0f, 0.0f, 1.0f);

/// 4-dimensional <code>float</code> Vector of zeros.
static const TVector<float, 4> kZeros4f(0.0f, 0.0f, 0.0f, 0.0f);
/// 4-dimensional <code>float</code> Vector of ones.
static const TVector<float, 4> kOnes4f(1.0f, 1.0f, 1.0f, 1.0f);
/// 4-dimensional <code>float</code> unit Vector pointing along the X axis.
static const TVector<float, 4> kAxisX4f(1.0f, 0.0f, 0.0f, 0.0f);
/// 4-dimensional <code>float</code> unit Vector pointing along the Y axis.
static const TVector<float, 4> kAxisY4f(0.0f, 1.0f, 0.0f, 0.0f);
/// 4-dimensional <code>float</code> unit Vector pointing along the Z axis.
static const TVector<float, 4> kAxisZ4f(0.0f, 0.0f, 1.0f, 0.0f);
/// 4-dimensional <code>float</code> unit Vector pointing along the W axis.
static const TVector<float, 4> kAxisW4f(0.0f, 0.0f, 0.0f, 1.0f);

/// 2-dimensional <code>double</code> Vector of zeros.
static const TVector<double, 2> kZeros2d(0.0, 0.0);
/// 2-dimensional <code>double</code> Vector of ones.
static const TVector<double, 2> kOnes2d(1.0, 1.0);
/// 2-dimensional <code>double</code> unit Vector pointing along the X axis.
static const TVector<double, 2> kAxisX2d(1.0, 0.0);
/// 2-dimensional <code>double</code> unit Vector pointing along the Y axis.
static const TVector<double, 2> kAxisY2d(0.0, 1.0);

/// 3-dimensional <code>double</code> Vector of zeros.
static const TVector<double, 3> kZeros3d(0.0, 0.0, 0.0);
/// 3-dimensional <code>double</code> Vector of ones.
static const TVector<double, 3> kOnes3d(1.0, 1.0, 1.0);
/// 3-dimensional <code>double</code> unit Vector pointing along the X axis.
static const TVector<double, 3> kAxisX3d(1.0, 0.0, 0.0);
/// 3-dimensional <code>double</code> unit Vector pointing along the Y axis.
static const TVector<double, 3> kAxisY3d(0.0, 1.0, 0.0);
/// 3-dimensional <code>double</code> unit Vector pointing along the Z axis.
static const TVector<double, 3> kAxisZ3d(0.0, 0.0, 1.0);

/// 4-dimensional <code>double</code> Vector of zeros.
static const TVector<double, 4> kZeros4d(0.0, 0.0, 0.0, 0.0);
/// 4-dimensional <code>double</code> Vector of ones.
static const TVector<double, 4> kOnes4d(1.0, 1.0, 1.0, 1.0);
/// 4-dimensional <code>double</code> unit Vector pointing along the X axis.
static const TVector<double, 4> kAxisX4d(1.0, 0.0, 0.0, 0.0);
/// 4-dimensional <code>double</code> unit Vector pointing along the Y axis.
static const TVector<double, 4> kAxisY4d(0.0, 1.0, 0.0, 0.0);
/// 4-dimensional <code>double</code> unit Vector pointing along the Z axis.
static const TVector<double, 4> kAxisZ4d(0.0, 0.0, 1.0, 0.0);
/// 4-dimensional <code>double</code> unit Vector pointing along the W axis.
static const TVector<double, 4> kAxisW4d(0.0, 0.0, 0.0, 1.0);

/// 2-dimensional <code>int</code> Vector of zeros.
static const TVector<int, 2> kOnes2i(1, 1);
/// 2-dimensional <code>int</code> Vector of ones.
static const TVector<int, 2> kZeros2i(0, 0);
/// 2-dimensional <code>int</code> unit Vector pointing along the X axis.
static const TVector<int, 2> kAxisX2i(1, 0);
/// 2-dimensional <code>int</code> unit Vector pointing along the Y axis.
static const TVector<int, 2> kAxisY2i(0, 1);

/// 3-dimensional <code>int</code> Vector of zeros.
static const TVector<int, 3> kZeros3i(0, 0, 0);
/// 3-dimensional <code>int</code> Vector of ones.
static const TVector<int, 3> kOnes3i(1, 1, 1);
/// 3-dimensional <code>int</code> unit Vector pointing along the X axis.
static const TVector<int, 3> kAxisX3i(1, 0, 0);
/// 3-dimensional <code>int</code> unit Vector pointing along the Y axis.
static const TVector<int, 3> kAxisY3i(0, 1, 0);
/// 3-dimensional <code>int</code> unit Vector pointing along the Z axis.
static const TVector<int, 3> kAxisZ3i(0, 0, 1);

/// 4-dimensional <code>int</code> Vector of zeros.
static const TVector<int, 4> kZeros4i(0, 0, 0, 0);
/// 4-dimensional <code>int</code> Vector of ones.
static const TVector<int, 4> kOnes4i(1, 1, 1 ,1);
/// 4-dimensional <code>int</code> unit Vector pointing along the X axis.
static const TVector<int, 4> kAxisX4i(1, 0, 0, 0);
/// 4-dimensional <code>int</code> unit Vector pointing along the Z axis.
static const TVector<int, 4> kAxisY4i(0, 1, 0, 0);
/// 4-dimensional <code>int</code> unit Vector pointing along the Y axis.
static const TVector<int, 4> kAxisZ4i(0, 0, 1, 0);
/// 4-dimensional <code>int</code> unit Vector pointing along the W axis.
static const TVector<int, 4> kAxisW4i(0, 0, 0, 1);

/// Quaternion Identity
static const TQuaternion<float> kQuatIdentityf(0.0f, 0.0f, 0.0f, 1.0f);
/// Quaternion Identity
static const TQuaternion<double> kQuatIdentityd(0.0, 0.0, 0.0, 1.0);

// An AffineTransform versoin of the mat4 Identity matrix.
static const TAffineTransform kAffineIdentity(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                             0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                             0.0f);

const double kAE = 1.0;
const double kQ0 = 120.0;
const double kS0 = 78.0;
const double kG = 6.67408e-11;                //万有引力常数(+/-0.00031)
const double kEarthMass = 5.97219e24;         //地球质量
const double kGM = 398603e9;                  //地心引力常数
const double kXKMPER = 6378.135;
const double kXJ2 = 1.082616e-3;
const double kXJ3 = -2.53881e-6;
const double kXJ4 = -1.65597e-6;

/*
 * alternative XKE
 * affects final results
 * aiaa-2006-6573
 * const double kXKE = 60.0 / sqrt(kXKMPER * kXKMPER * kXKMPER / kMU);
 * dundee
 * const double kXKE = 7.43669161331734132e-2;
 */
const double kXKE = 60.0 / sqrt(kXKMPER * kXKMPER * kXKMPER / kGM);
const double kCK2 = 0.5 * kXJ2 * kAE * kAE;
const double kCK4 = -0.375 * kXJ4 * kAE * kAE * kAE * kAE;

/*
 * alternative QOMS2T
 * affects final results
 * aiaa-2006-6573
 * #define QOMS2T   (pow(((Q0 - S0) / XKMPER), 4.0))
 * dundee
 * #define QOMS2T   (1.880279159015270643865e-9)
 */
const double kQOMS2T = pow(((kQ0 - kS0) / kXKMPER), 4.0);

const double kS = kAE * (1.0 + kS0 / kXKMPER);
const double kPI = 3.14159265358979323846264338327950288419716939937510582;
const double kTWOPI = 2.0 * kPI;
const double kTWOTHIRD = 2.0 / 3.0;
const double kTHDT = 4.37526908801129966e-3;
/*
 * earth flattening
 */
const double kF = 1.0 / 298.26;
/*
 * earth rotation per sideral day
 */
const double kOMEGA_E = 1.00273790934;
const double kAU = 1.49597870691e8;

const double kSECONDS_PER_DAY = 86400.0;
const double kMINUTES_PER_DAY = 1440.0;
const double kHOURS_PER_DAY = 24.0;

// Jan 1.0 1900 = Jan 1 1900 00h UTC
const double kEPOCH_JAN1_00H_1900 = 2415019.5;

// Jan 1.5 1900 = Jan 1 1900 12h UTC
const double kEPOCH_JAN1_12H_1900 = 2415020.0;

// Jan 1.5 2000 = Jan 1 2000 12h UTC
const double kEPOCH_JAN1_12H_2000 = 2451545.0;

/// @}

TMATH_END_NAMESPACE


#endif // TCONSTANTS_H

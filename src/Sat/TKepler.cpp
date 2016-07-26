#include "TKepler.h"

#include <cmath>
#include <iostream>

#include "TECI.h"

TMATH_BEGIN_NAMESPACE

TKepler::TKepler()
{

}

TKepler::~TKepler()
{

}

void TKepler::elemToCoor( const TVector<double, 6> &elements, const double &time, 
    TVector<double, 3> &rVector, TVector<double, 3> &vVector, const double GM)
{
    ///开普勒六根数
    double a = elements[0];
    double e = elements[1];
    double i = elements[2];
    double Omega = elements[3];
    double omega = elements[4];
    double M0 = elements[5];

    //计算平近点角位置
    double M = M0 + sqrt(GM / (a * a *a)) * time;

    //计算偏近点角
    double E = eccAnom(M, e);

    double cosE = cos(E);
    double sinE = sin(E);

    //轨道平面坐标系
    double fac = sqrt((1.0-e) * (1.0 +e));

    double R = a * (1.0 - e * cosE);  //距离
    double V = sqrt(GM * a) / R;      //速度

    TVector<double, 3> r(a * (cosE - e), a * fac * sinE, 0.0);
    TVector<double, 3> v(-V * sinE, V * fac * cosE, 0.0);

    //轨道面到惯性系旋转矩阵
    TMatrix<double, 3, 3> MR = TMatrix<double, 3, 3>::RotationZ(-Omega) * TMatrix<double, 3, 3>::RotationX(-i);
    MR *=  TMatrix<double, 3, 3>::RotationZ(-omega);

    rVector = MR * r;
    vVector = MR * v;
}

bool TKepler::coorToElem(const TVector<double, 3> &rVector, const TVector<double, 3> &vVector,
    TVector<double, 6> &elements, const double GM)
{

    TVector<double, 3> hVector = TVector<double, 3>::crossProduct(rVector, vVector);
    double H = hVector.length();

    double Omega = atan2 ( hVector[0], -hVector[1] );                     /// 升交点赤经
    Omega = fmod(Omega, kTWOPI);
    double i     = atan2 ( sqrt(hVector[0] * hVector[0] + hVector[1] * hVector[1]), hVector[2] ); /// 轨道倾角
    double u     = atan2 ( rVector[2] * H, - rVector[0] * hVector[1] + rVector[1] * hVector[0] ); /// Arg. of latitude

    double R  = rVector.length();                                                            /// 地心距离

    double a = R / (2.0 - R * TVector<double, 3>::dotProduct(vVector, vVector) / GM);        /// 长半轴

    double eCosE = 1.0 - R / a;                                                              /// e*cos(E)
    double eSinE = TVector<double, 3>::dotProduct(rVector, vVector) / sqrt( GM * a);         /// e*sin(E)

    double e2 = eCosE * eCosE + eSinE * eSinE;
    double e  = sqrt(e2);                                      /// 偏心率
    double E  = atan2(eSinE, eCosE);                           /// 偏近点角

    double M  = Mod( E - eSinE, kTWOPI);                      /// 平近点角

    double nu = atan2(sqrt(1.0 - e2) * eSinE, eCosE - e2);     /// 真近点角

    double omega = Mod( u - nu, kTWOPI);                      /// 近地点幅角

    /// 返回开普勒六根数
    elements[0] = a;
    elements[1] = e;
    elements[2] = i;
    elements[3] = Omega;
    elements[4] = omega;
    elements[5] = M;

    return true;
}

double TKepler::eccAnom( double M, double e )
{
    // Constants
    const int maxit = 15;

    // Variables

    int    i=0;
    double E, f;

    // Starting value

    M = fmod(M, kTWOPI);
    if (e<0.8)
    {
        E = M;
    }
    else
    {
        E = kPI;
    }

    /// 牛顿迭代
    do
    {
        f = E - e * sin(E) - M;
        E = E - f / ( 1.0 - e * cos(E) );
        ++i;
        if (i == maxit)
        {
            std::cerr << " convergence problems in EccAnom" << std::endl;
            break;
        }
    }
    while (fabs(f) > 1.0e-15);

    return(E);
}

TMATH_END_NAMESPACE

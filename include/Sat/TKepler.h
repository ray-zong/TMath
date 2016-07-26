#ifndef TKEPLER_H
#define TKEPLER_H

#include "TGlobal.h"
#include "TVector.h"
#include "TConstants.h"
#include "TECI.h"

TMATH_BEGIN_NAMESPACE
/**
 * @brief 
 */
class TMATH_EXPORT TKepler
{
public:
    TKepler();
    ~TKepler();

    /** 
     * @brief 解开普勒方程
     * @param elements 开普勒六根数    (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                                  [m, ,rad,rad,rad,rad]
     * @param time  相对于历元的时间       [s]
     * @param GM  引力常量,默认为地球引力常量
     * @return eci   ECI坐标系下位置速度
     */
    static void elemToCoor(const TVector<double, 6> &elements, const double &time,
                      TVector<double, 3> &rVector, TVector<double, 3> &vVector, const double GM = kGM);

    /** 
     * @brief 由位置、方向矢量求轨道根数
     * @param rVector  位置矢量
     * @param vVector  方向矢量
     * @param GM  引力常量,默认为地球引力常量
     * @return elements 开普勒六根数    (半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角)
     *                                  [m, ,rad,rad,rad,rad]
     */
    static bool coorToElem(const TVector<double, 3> &rVector, const TVector<double, 3> &vVector,
        TVector<double, 6> &elements, const double GM = kGM);

private:
    /** 
     * @brief 计算偏近点角
     * @param M 平近点角
     * @param e 离心率
     * @return double 偏近点角
     */
    static double eccAnom(double M, double e);
};

TMATH_END_NAMESPACE

#endif

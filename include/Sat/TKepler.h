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
     * @brief �⿪���շ���
     * @param elements ������������    (�볤��,������,������,������ྭ,���ص����,ƽ�����)
     *                                  [m, ,rad,rad,rad,rad]
     * @param time  �������Ԫ��ʱ��       [s]
     * @param GM  ��������,Ĭ��Ϊ������������
     * @return eci   ECI����ϵ��λ���ٶ�
     */
    static void elemToCoor(const TVector<double, 6> &elements, const double &time,
                      TVector<double, 3> &rVector, TVector<double, 3> &vVector, const double GM = kGM);

    /** 
     * @brief ��λ�á�����ʸ����������
     * @param rVector  λ��ʸ��
     * @param vVector  ����ʸ��
     * @param GM  ��������,Ĭ��Ϊ������������
     * @return elements ������������    (�볤��,������,������,������ྭ,���ص����,ƽ�����)
     *                                  [m, ,rad,rad,rad,rad]
     */
    static bool coorToElem(const TVector<double, 3> &rVector, const TVector<double, 3> &vVector,
        TVector<double, 6> &elements, const double GM = kGM);

private:
    /** 
     * @brief ����ƫ�����
     * @param M ƽ�����
     * @param e ������
     * @return double ƫ�����
     */
    static double eccAnom(double M, double e);
};

TMATH_END_NAMESPACE

#endif

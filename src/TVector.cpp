#include "TVector.h"

#include <cassert>

RAY_BEGIN_NAMESPACE

template<class T, int d>
        TVector<T, d>::TVector()
{

}

template<class T, int d>
TVector<T, d>::TVector(const TVector<T, d> &v)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] = v.m_data[i];
    }
}

template<class T, int d>
TVector<T, d>::TVector(const T &s)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] = s;
    }
}

template<class T, int d>
TVector<T, d>::TVector(const T *a)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] = a[i];
    }
}

template<class T, int d>
TVector<T, d>::TVector(const T &s1, const T &s2, const T &s3)
{
    assert(d == 3);
    m_data[0] = s1;
    m_data[1] = s2;
    m_data[2] = s3;
}

template<class T, int d>
TVector<T, d>::TVector(const T &s1, const T &s2, const T &s3, const T &s4)
{
    assert(d == 4);
    m_data[0] = s1;
    m_data[1] = s2;
    m_data[2] = s3;
    m_data[3] = s4;
}

template<class T, int d>
T &TVector<T, d>::operator()(const int i)
{
    return m_data[i];
}

template<class T, int d>
const T &TVector<T, d>::operator()(const int i) const
{
    return m_data[i];
}

template<class T, int d>
T &TVector<T, d>::operator[](const int i)
{
    return m_data[i];
}

template<class T, int d>
const T &TVector<T, d>::operator[](const int i) const
{
    return m_data[i];
}

template<class T, int d>
T &TVector<T, d>::x()
{
    return m_data[0];
}

template<class T, int d>
T &TVector<T, d>::y()
{
    return m_data[1];
}

template<class T, int d>
T &TVector<T, d>::z()
{
    return m_data[2];
}

template<class T, int d>
T &TVector<T, d>::w()
{
    return m_data[3];
}

template<class T, int d>
const T &TVector<T, d>::x() const
{
    return m_data[0];
}

template<class T, int d>
const T &TVector<T, d>::y() const
{
    return m_data[1];
}

template<class T, int d>
const T &TVector<T, d>::z() const
{
    return m_data[2];
}

template<class T, int d>
const T &TVector<T, d>::w() const
{
    return m_data[3];
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator-() const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i];
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator*(const TVector<T, d> &v) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] * v.m_data[i];
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator/(const TVector<T, d> &v) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] / v.m_data[i];
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator+(const TVector<T, d> &v) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] + v.m_data[i];
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator-(const TVector<T, d> &v) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] - v.m_data[i];
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator*(const T &s) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] * s;
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator/(const T &s) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] / s;
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator+(const T &s) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] + s;
    }
    return result;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::operator-(const T &s) const
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = m_data[i] - s;
    }
    return result;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator*=(const TVector<T, d> &v)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] *= v.m_data[i];
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator/=(const TVector<T, d> &v)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] /= v.m_data[i];
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator+=(const TVector<T, d> &v)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] += v.m_data[i];
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator-=(const TVector<T, d> &v)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] -= v.m_data[i];
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator*=(const T &s)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] *= s;
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator/=(const T &s)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] /= s;
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator+=(const T &s)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] += s;
    }
    return *this;
}

template<class T, int d>
TVector<T, d> &TVector<T, d>::operator-=(const T &s)
{
    for(int i = 0; i < d; ++i)
    {
        m_data[i] -= s;
    }
    return *this;
}

template<class T, int d>
T TVector<T, d>::lengthSquared() const
{
    return dotProduct(*this, *this);
}

template<class T, int d>
T TVector<T, d>::length() const
{
    return sqrt(lengthSquared());
}

template<class T, int d>
T TVector<T, d>::normalize()
{
    const T length = length();
    *this = *this * (1 / length);
    return length;
}

template<class T, int d>
TVector<T, d> TVector<T, d>::normalized() const
{
    return *this * (1 / length());
}

template<class T, int d>
T TVector<T, d>::dotProduct(const TVector<T, d> &v1,
                            const TVector<T, d> &v2)
{
    dotProductHelper(v1, v2);
}

template<class T, int d>
TVector<T, d> TVector<T, d>::hadamardProduct(const TVector<T, d> &v1,
                                             const TVector<T, d> &v2)
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = v1[i] * v2[i];
    }
    return result;
}

TVector<T, 3> TVector<T, d>::crossProduct(const TVector<T, 3> &v1,
                                          const TVector<T, 3> &v2)
{
    return TVector<T, 3>(v1[1] * v2[2] - v1[2] * v2[1],
                         v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]);
}

TVector<T, d> TVector<T, d>::lerp(const TVector<T, d> &v1,
                                  const TVector<T, d> &v2,
                                  const T percent)
{
    const T one_minus_percent = static_cast<T>(1.0) - percent;

    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = one_minus_percent * v1[i]
                + percent * v2[i];
    }
    return result;
}

TVector<T, d> TVector<T, d>::randomInRange(const TVector<T, d> &min,
                                           const TVector<T, d> &max)
{
    TVector<T, d> result;
    for(int i = 0; i < d; ++i)
    {
        result[i] = one_minus_percent * v1[i]
                + percent * v2[i];
    }
    return result;
}

template<class T, int d>
T TVector<T, d>::dotProductHelper(const TVector<T, d> &v1,
                                  const TVector<T, d> &v2)
{
    T result = 0;
    for(int i = 0; i < d; ++i)
    {
        result += v1[i] * v2[i];
    }
    return result;
}

template<class T, int d>
T TVector<T, d>::dotProductHelper(const TVector<T, 2> &v1,
                                  const TVector<T, 2> &v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

template<class T, int d>
T TVector<T, d>::dotProductHelper(const TVector<T, 3> &v1,
                                  const TVector<T, 3> &v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

template<class T, int d>
T TVector<T, d>::dotProductHelper(const TVector<T, 4> &v1,
                                  const TVector<T, 4> &v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1]
            + v1[2] * v2[2] + v1[3] * v2[3];
}
//template<class T, int d, typename U>
//TVector<T, d>::TVector(const TVector<U, d> &v)
//{
//    for(int i = 0; i < d; ++i)
//    {
//        m_data[i] = static_cast<T>(v.m_data[i]);
//    }
//}

RAY_END_NAMESPACE

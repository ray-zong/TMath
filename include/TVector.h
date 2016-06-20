#ifndef TVECTOR_H
#define TVECTOR_H

#include "TGlobal.h"

RAY_BEGIN_NAMESPACE

template<class T, int d>
class TVector
{
public:
    inline TVector();

    inline TVector(const TVector<T, d> &v);

    //template <typename U>
    //explicit inline TVector(const TVector<U, d> &v);

    explicit inline TVector(const T &s);

    explicit inline TVector(const T *a);

    inline TVector(const T &s1, const T &s2, const T &s3);

    inline TVector(const T &s1, const T &s2, const T &s3, const T &s4);

    inline T &operator()(const int i);

    inline const T &operator()(const int i) const;

    inline T &operator[](const int i);

    inline const T &operator[](const int i) const;

    inline T &x();

    inline T &y();

    inline T &z();

    inline T &w();

    inline const T &x() const;

    inline const T &y() const;

    inline const T &z() const;

    inline const T &w() const;

    inline TVector<T, d> operator-() const;

    inline TVector<T, d> operator*(const TVector<T, d> &v) const;

    inline TVector<T, d> operator/(const TVector<T, d> &v) const;

    inline TVector<T, d> operator+(const TVector<T, d> &v) const;

    inline TVector<T, d> operator-(const TVector<T, d> &v) const;

    inline TVector<T, d> operator*(const T &s) const;

    inline TVector<T, d> operator/(const T &s) const;

    inline TVector<T, d> operator+(const T &s) const;

    inline TVector<T, d> operator-(const T &s) const;

    inline TVector<T, d> &operator*=(const TVector<T, d> &v);

    inline TVector<T, d> &operator/=(const TVector<T, d> &v);

    inline TVector<T, d> &operator+=(const TVector<T, d> &v);

    inline TVector<T, d> &operator-=(const TVector<T, d> &v);

    inline TVector<T, d> &operator*=(const T &s);

    inline TVector<T, d> &operator/=(const T &s);

    inline TVector<T, d> &operator+=(const T &s);

    inline TVector<T, d> &operator-=(const T &s);

    inline T lengthSquared() const;

    inline T length() const;

    inline T normalize();

    inline TVector<T, d> normalized() const;

    static inline T dotProduct(const TVector<T, d> &v1,
                               const TVector<T, d> &v2);

    static inline TVector<T, d> hadamardProduct(const TVector<T, d> &v1,
                                                const TVector<T, d> &v2);

    static inline TVector<T, 3> crossProduct(const TVector<T, 3> &v1,
                                             const TVector<T, 3> &v2);

    static inline TVector<T, d> lerp(const TVector<T, d> &v1,
                                     const TVector<T, d> &v2,
                                     const T percent);

    static inline TVector<T, d> randomInRange(const TVector<T, d> &min,
                                              const TVector<T, d> &max);


    static inline T dotProductHelper(const TVector<T, d> &v1,
                                     const TVector<T, d> &v2);

    static inline T dotProductHelper(const TVector<T, 2> &v1,
                                     const TVector<T, 2> &v2);

    static inline T dotProductHelper(const TVector<T, 3> &v1,
                                     const TVector<T, 3> &v2);

    static inline T dotProductHelper(const TVector<T, 4> &v1,
                                     const TVector<T, 4> &v2);
private:
    T m_data[d];
};

RAY_END_NAMESPACE

#endif // TVECTOR_H

#ifndef TMATRIX_H
#define TMATRIX_H

#include "TGlobal.h"

RAY_BEGIN_NAMESPACE

template <typename T, int rows, int columns>
class TMatrix
{
public:
    inline TMatrix();
    inline TMatrix(const TMatrix<T, rows, columns> &m);
    explicit TMatrix(const T &values);

private:
    T m_data[rows][columns];
};

RAY_END_NAMESPACE

#endif // TMATRIX_H

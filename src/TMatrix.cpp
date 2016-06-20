#include "TMatrix.h"

RAY_BEGIN_NAMESPACE

template <typename T, int rows, int columns>
TMatrix<T, rows, columns>::TMatrix()
{

}

template <typename T, int rows, int columns>
TMatrix<T, rows, columns>::TMatrix(const TMatrix<Ray::T, Ray::rows, Ray::columns> &m)
{
}

template <typename T, int rows, int columns>
TMatrix<T, rows, columns>::TMatrix(const Ray::T &values)
{
}

RAY_END_NAMESPACE

#ifndef TSATELLITEEXCEPTION_H
#define TSATELLITEEXCEPTION_H

#include <stdexcept>
#include <string>

#include "TGlobal.h"

TMATH_BEGIN_NAMESPACE
/**
 * @brief The exception that the SGP4 class throws upon an error.
 */
class TSatelliteException : public std::runtime_error
{
public:
    TSatelliteException(const char* message)
        : runtime_error(message)
    {
    }

    virtual ~TSatelliteException(){}
};

TMATH_END_NAMESPACE

#endif



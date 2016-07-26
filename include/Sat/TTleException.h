#ifndef TTLEEXCEPTION_H
#define TTLEEXCEPTION_H

#include <stdexcept>
#include <string>

#include "TGlobal.h"

TMATH_BEGIN_NAMESPACE
/**
 * @brief The exception that the Tle class throws on an error.
 *
 * The exception that the Tle decoder will throw on an error.
 */
class TTleException : public std::runtime_error
{
public:
    /**
     * Constructor
     * @param message Exception message
     */
    TTleException(const char* message)
        : runtime_error(message)
    {
    }

    //TTleException(const TTleException&) = default;

    virtual ~TTleException(){}
};

TMATH_END_NAMESPACE

#endif
